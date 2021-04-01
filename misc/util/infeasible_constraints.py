# -*- coding: utf-8 -*-
#  ___________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright 2017 National Technology and Engineering Solutions of Sandia, LLC
#  Under the terms of Contract DE-NA0003525 with National Technology and
#  Engineering Solutions of Sandia, LLC, the U.S. Government retains certain
#  rights in this software.
#  This software is distributed under the 3-clause BSD License.
#  ___________________________________________________________________________

"""Module with diagnostic utilities for infeasible models."""
from pyomo.environ import ConcreteModel, Constraint, Var, inequality
from pyomo.common.log import LoggingIntercept
import pyutilib.th as unittest
from six import StringIO
from pyomo.core import Constraint, Var, value
from math import fabs
import logging

from pyomo.common import deprecated
from pyomo.core.expr.visitor import identify_variables
from pyomo.util.blockutil import log_model_constraints

logger = logging.getLogger(__name__)

def log_infeasible_constraints(
        m, tol=1E-6, logger=logger,
        log_expression=False, log_variables=False
):
    """Print the infeasible constraints in the model.

    Uses the current model state. Uses pyomo.util.infeasible logger unless one
    is provided.

    Args:
        m (Block): Pyomo block or model to check
        tol (float): feasibility tolerance
        log_expression (bool): If true, prints the constraint expression
        log_variables (bool): If true, prints the constraint variable names and values

    """
    # Iterate through all active constraints on the model
    for constr in m.component_data_objects(
            ctype=Constraint, active=True, descend_into=True):
        constr_body_value = value(constr.body, exception=False)
        constr_lb_value = value(constr.lower, exception=False)
        constr_ub_value = value(constr.upper, exception=False)

        constr_undefined = False
        equality_violated = False
        lb_violated = False
        ub_violated = False

        if constr_body_value is None:
            # Undefined constraint body value due to missing variable value
            constr_undefined = True
            pass
        else:
            # Check for infeasibilities
            if constr.equality:
                if fabs(constr_lb_value - constr_body_value) >= tol:
                    equality_violated = True
            else:
                if constr.has_lb() and constr_lb_value - constr_body_value >= tol:
                    lb_violated = True
                if constr.has_ub() and constr_body_value - constr_ub_value >= tol:
                    ub_violated = True

        if not any((constr_undefined, equality_violated, lb_violated, ub_violated)):
            # constraint is fine. skip to next constraint
            continue

        output_dict = dict(name=constr.name)

        log_template = "CONSTR {name}: {lb_value}{lb_operator}{body_value}{ub_operator}{ub_value}"
        if log_expression:
            log_template += "\n  - EXPR: {lb_expr}{lb_operator}{body_expr}{ub_operator}{ub_expr}"
        if log_variables:
            vars_template = "\n  - VAR {name}: {value}"
            log_template += "{var_printout}"
            constraint_vars = identify_variables(constr.body, include_fixed=True)
            output_dict['var_printout'] = ''.join(
                vars_template.format(name=v.name, value=v.value) for v in constraint_vars)

        output_dict['body_value'] = "missing variable value" if constr_undefined else constr_body_value
        output_dict['body_expr'] = constr.body
        if constr.equality:
            output_dict['lb_value'] = output_dict['lb_expr'] = output_dict['lb_operator'] = ""
            output_dict['ub_value'] = constr_ub_value
            output_dict['ub_expr'] = constr.upper
            if equality_violated:
                output_dict['ub_operator'] = " =/= "
            elif constr_undefined:
                output_dict['ub_operator'] = " =?= "
        else:
            if constr.has_lb():
                output_dict['lb_value'] = constr_lb_value
                output_dict['lb_expr'] = constr.lower
                if lb_violated:
                    output_dict['lb_operator'] = " </= "
                elif constr_undefined:
                    output_dict['lb_operator'] = " <?= "
                else:
                    output_dict['lb_operator'] = " <= "
            else:
                output_dict['lb_value'] = output_dict['lb_expr'] = output_dict['lb_operator'] = ""

            if constr.has_ub():
                output_dict['ub_value'] = constr_ub_value
                output_dict['ub_expr'] = constr.upper
                if ub_violated:
                    output_dict['ub_operator'] = " </= "
                elif constr_undefined:
                    output_dict['ub_operator'] = " <?= "
                else:
                    output_dict['ub_operator'] = " <= "
            else:
                output_dict['ub_value'] = output_dict['ub_expr'] = output_dict['ub_operator'] = ""

        logger.info(log_template.format(**output_dict))


def log_infeasible_bounds(m, tol=1E-6, logger=logger):
    """Print the infeasible variable bounds in the model.

    Args:
        m (Block): Pyomo block or model to check
        tol (float): feasibility tolerance

    """
    for var in m.component_data_objects(
            ctype=Var, descend_into=True):
        var_value = var.value
        if var_value is None:
            logger.debug("Skipping VAR {} with no assigned value.")
            continue
        if var.has_lb() and value(var.lb - var) >= tol:
            logger.info('VAR {}: {} >/= LB {}'.format(
                var.name, value(var), value(var.lb)))
        if var.has_ub() and value(var - var.ub) >= tol:
            logger.info('VAR {}: {} </= UB {}'.format(
                var.name, value(var), value(var.ub)))


def log_close_to_bounds(m, tol=1E-6, logger=logger):
    """Print the variables and constraints that are near their bounds.

    Fixed variables and equality constraints are excluded from this analysis.

    Args:
        m (Block): Pyomo block or model to check
        tol (float): bound tolerance
    """
    for var in m.component_data_objects(
            ctype=Var, descend_into=True):
        if var.fixed:
            continue
        var_value = var.value
        if var_value is None:
            logger.debug("Skipping VAR {} with no assigned value.")
            continue
        if (var.has_lb() and var.has_ub() and
                fabs(value(var.ub - var.lb)) <= 2 * tol):
            continue  # if the bounds are too close, skip.
        if var.has_lb() and fabs(value(var.lb - var)) <= tol:
            logger.info('{} near LB of {}'.format(var.name, value(var.lb)))
        elif var.has_ub() and fabs(value(var.ub - var)) <= tol:
            logger.info('{} near UB of {}'.format(var.name, value(var.ub)))

    for constr in m.component_data_objects(
            ctype=Constraint, descend_into=True, active=True):
        if not constr.equality:
            # skip equality constraints, because they should always be close to
            # bounds if enforced.
            body_value = value(constr.body, exception=False)
            if body_value is None:
                logger.info("Skipping CONSTR {}: missing variable value.".format(constr.name))
                continue
            if (constr.has_ub() and
                    fabs(value(body_value - constr.upper)) <= tol):
                logger.info('{} near UB'.format(constr.name))
            if (constr.has_lb() and
                    fabs(value(body_value - constr.lower)) <= tol):
                logger.info('{} near LB'.format(constr.name))

"""Tests infeasible model debugging utilities."""


class TestInfeasible(unittest.TestCase):
    """Tests infeasible model debugging utilities."""

    def build_model(self):
        m = ConcreteModel()
        m.x = Var(initialize=1)
        m.c1 = Constraint(expr=m.x >= 2)
        m.c2 = Constraint(expr=m.x == 4)
        m.c3 = Constraint(expr=m.x <= 0)
        m.y = Var(bounds=(0, 2), initialize=1.9999999)
        m.c4 = Constraint(expr=m.y >= m.y.value)
        m.z = Var(bounds=(0, 6))
        m.c5 = Constraint(expr=inequality(5, m.z, 10), doc="Range infeasible")
        m.c6 = Constraint(expr=m.x + m.y <= 6, doc="Feasible")
        m.c7 = Constraint(expr=m.z == 6, doc="Equality infeasible")
        m.c8 = Constraint(expr=inequality(3, m.x, 6),
                          doc="Range lb infeasible")
        m.c9 = Constraint(expr=inequality(0, m.x, 0.5),
                          doc="Range ub infeasible")
        m.c10 = Constraint(expr=m.y >= 3, doc="Inactive")
        m.c10.deactivate()
        m.c11 = Constraint(expr=m.y <= m.y.value)
        m.yy = Var(bounds=(0, 1), initialize=1E-7, doc="Close to lower bound")
        m.y3 = Var(bounds=(0, 1E-7), initialize=0, doc="Bounds too close")
        m.y4 = Var(bounds=(0, 1), initialize=2, doc="Fixed out of bounds.")
        m.y4.fix()
        return m

    def test_log_infeasible_constraints(self):
        """Test for logging of infeasible constraints."""
        m = self.build_model()
        output = StringIO()
        with LoggingIntercept(output, 'pyomo.util.infeasible', logging.INFO):
            log_infeasible_constraints(m)
        expected_output = [
            "CONSTR c1: 2.0 </= 1",
            "CONSTR c2: 1 =/= 4.0",
            "CONSTR c3: 1 </= 0.0",
            "CONSTR c5: 5.0 <?= missing variable value <?= 10.0",
            "CONSTR c7: missing variable value =?= 6.0",
            "CONSTR c8: 3.0 </= 1 <= 6.0",
            "CONSTR c9: 0.0 <= 1 </= 0.5",
        ]
        self.assertEqual(expected_output, output.getvalue().splitlines())

    def test_log_infeasible_bounds(self):
        """Test for logging of infeasible variable bounds."""
        m = self.build_model()
        m.x.setlb(2)
        m.x.setub(0)
        output = StringIO()
        with LoggingIntercept(output, 'pyomo.util', logging.INFO):
            log_infeasible_bounds(m)
        expected_output = [
            "VAR x: 1 >/= LB 2", "VAR x: 1 </= UB 0", "VAR y4: 2 </= UB 1",
        ]
        self.assertEqual(expected_output, output.getvalue().splitlines())

    def test_log_active_constraints(self):
        """Test for logging of active constraints."""
        m = self.build_model()
        depr = StringIO()
        output = StringIO()
        with LoggingIntercept(depr, 'pyomo'):
            with LoggingIntercept(output, 'pyomo.util', logging.INFO):
                log_active_constraints(m)
        self.assertIn("log_active_constraints is deprecated.", depr.getvalue())
        expected_output = [
            "c1 active", "c2 active", "c3 active", "c4 active",
            "c5 active", "c6 active", "c7 active", "c8 active",
            "c9 active", "c11 active"
        ]
        self.assertEqual(expected_output, output.getvalue().splitlines())

    def test_log_close_to_bounds(self):
        """Test logging of variables and constraints near bounds."""
        m = self.build_model()
        output = StringIO()
        with LoggingIntercept(output, 'pyomo.util.infeasible', logging.INFO):
            log_close_to_bounds(m)
        expected_output = [
            "y near UB of 2", "yy near LB of 0", "c4 near LB",
            "Skipping CONSTR c5: missing variable value.",
            "c11 near UB",
        ]
        self.assertEqual(expected_output, output.getvalue().splitlines())

    def test_log_infeasible_constraints_verbose_expressions(self):
        """Test for logging of infeasible constraints."""
        m = self.build_model()
        output = StringIO()
        with LoggingIntercept(output, 'pyomo.util.infeasible', logging.INFO):
            log_infeasible_constraints(m, log_expression=True)
        expected_output = [
            "CONSTR c1: 2.0 </= 1", "  - EXPR: 2.0 </= x",
            "CONSTR c2: 1 =/= 4.0", "  - EXPR: x =/= 4.0",
            "CONSTR c3: 1 </= 0.0", "  - EXPR: x </= 0.0",
            "CONSTR c5: 5.0 <?= missing variable value <?= 10.0", "  - EXPR: 5.0 <?= z <?= 10.0",
            "CONSTR c7: missing variable value =?= 6.0", "  - EXPR: z =?= 6.0",
            "CONSTR c8: 3.0 </= 1 <= 6.0", "  - EXPR: 3.0 </= x <= 6.0",
            "CONSTR c9: 0.0 <= 1 </= 0.5", "  - EXPR: 0.0 <= x </= 0.5",
        ]
        self.assertEqual(expected_output, output.getvalue().splitlines())

    def test_log_infeasible_constraints_verbose_variables(self):
        """Test for logging of infeasible constraints."""
        m = self.build_model()
        output = StringIO()
        with LoggingIntercept(output, 'pyomo.util.infeasible', logging.INFO):
            log_infeasible_constraints(m, log_variables=True)
        expected_output = [
            "CONSTR c1: 2.0 </= 1", "  - VAR x: 1",
            "CONSTR c2: 1 =/= 4.0", "  - VAR x: 1",
            "CONSTR c3: 1 </= 0.0", "  - VAR x: 1",
            "CONSTR c5: 5.0 <?= missing variable value <?= 10.0", "  - VAR z: None",
            "CONSTR c7: missing variable value =?= 6.0", "  - VAR z: None",
            "CONSTR c8: 3.0 </= 1 <= 6.0", "  - VAR x: 1",
            "CONSTR c9: 0.0 <= 1 </= 0.5", "  - VAR x: 1",
        ]
        self.assertEqual(expected_output, output.getvalue().splitlines())


if __name__ == '__main__':
    unittest.main()
