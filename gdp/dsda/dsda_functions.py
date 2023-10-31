"""
dsda_functions.py
The function performs the Discrete-Steepes Descent Algorithm (D-SDA) for a given GDP model. 
The function allows the user to provide additional logical constraints and can optionally convert the model to its equivalent MINLP form using specified transformations. 
The function also allows the user to provide a MINLP solver to solve the subproblem model. The function returns the solution of the GDP model, the solution of the MINLP model, and the solution of the subproblem model.
The function also returns the D-SDA convergence plot.

1. The code contains the get_external_information function, which is used to obtain information from the model to perform the reformulation with external variables.
2. The code have the external_ref function, which is used to reformulate a given GDP model by taking into account external variables. 
3. dsda_functions.py has extvars_gdp_to_mip function used to transformed a GDP into MINLP model. The function maps the external variables defined for the GDP model into corresponding binary variables suitable for the MINLP formulation.  
4. With preprocess_problem and solve_subproblem functions, the code checks feasibility and solves the subproblem model.
5. The code has solve_with_minlp and solve_with_gdpopt functions to solve GDP model using MINLP and GDPopt respectively.
6. Neighborhood search is implemented in the neighborhood_k_eq_2 function and neighborhood_k_eq_inf function.
7. The model is initialized in the initialize_model and generate_initialize function.
8. The code has the dsda function, which is used to solve the GDP model using D-SDA.
9. The funtcion finds the neighborhood of the current solution and solves the subproblem model and evaluates with the neighbors.
10. The function do a line search after the neighborhood search.
11. The  solve_with_dsda function solves the GDP model using D-SDA. All the mentioned functions are used in this function.
12. visualize_dsda function draws the D-SDA convergence plot.
13. solve_complete_external_enumeration function computes complete enumeration using the external variable reformulation.

References:
[1] Bernal, David E., et al. "Process Superstructure Optimization through Discrete Steepest Descent Optimization: a GDP Analysis and Applications in Process Intensification." Computer Aided Chemical Engineering. Vol. 49. Elsevier, 2022. 1279-1284.
[2] Linan, David A., et al. "Optimal design of superstructures for placing units and streams with multiple and ordered available locations. Part I: A new mathematical framework." Computers & Chemical Engineering 137, (2020): 106794.


"""
import copy
import csv
import itertools as it
import os
import time
from math import isnan

import matplotlib.pyplot as plt
import numpy as np
import pyomo.environ as pe
from gdp.dsda.model_serializer import StoreSpec, from_json, to_json
from pyomo.common.errors import InfeasibleConstraintException
from pyomo.contrib.fbbt.fbbt import fbbt
from pyomo.contrib.gdpopt.data_class import MasterProblemResult
from pyomo.core.base.misc import display
from pyomo.core.plugins.transform.logical_to_linear import (
    update_boolean_vars_from_binary,
)
from pyomo.gdp import Disjunct, Disjunction
from pyomo.opt import SolutionStatus, SolverResults
from pyomo.opt import TerminationCondition as tc
from pyomo.opt.base.solvers import SolverFactory


def get_external_information(m: pe.ConcreteModel(), ext_ref, tee: bool = False):
    """
    Function that obtains information from the model to perform the reformulation with external variables.
    The model must be a GDP problem with exactly one "Exactly(k_j, [Y_j1,Y_j2,Y_j3,...])" constraint for each list of variables
    [Y_j1,Y_j2,Y_j3,...] that is going to be reformulated over set j.
    Args:
        m (pyomo.ConcreteModel): GDP model that is going to be reformulated
        ext_ref (dict): Dictionary with Boolean variables to be reformulated (keys) and their corresponding ordered sets (values). Both keys and values are pyomo objects.
        tee (bool): Display reformulation
    Returns:
        reformulation_dict (dict): A dictionary of dictionaries that looks as follows:
            {1:{'exactly_number':Number of external variables for this type,
                'Boolean_vars_names':list with names of the ordered Boolean variables to be reformulated,
                'Boolean_vars_ordered_index': Indexes where the external reformulation is applied,
                'Ext_var_lower_bound': Lower bound for this type of external variable,
                'Ext_var_upper_bound': Upper bound for this type of external variable },
             2:{...},...}

            The first key (positive integer) represent a type of external variable identified in the model. For this type of external variable
            a dictionary is created.
        number_of_external_variables: Number of external variables
        lower_bounds (dict): Dictionary with positive integer keys identifying the external variable, and its lower bound as value
        upper_bounds (dict): Dictionary with positive integer keys identifying the external variable, and its upper bound as value

    """

    # If Boolean variables that are going to be reformulated are defined over multiple sets try:
    try:
        # index of the set where reformultion can be applied for a given boolean variable
        ref_index = {}
        # index of the sets where the reformulation cannot be applied for a given boolean variable
        no_ref_index = {}
        for i in ext_ref:
            ref_index[i] = []
            no_ref_index[i] = []
            for index_set in range(len(i.index_set()._sets)):
                if i.index_set()._sets[index_set].name == ext_ref[i].name:
                    ref_index[i].append(index_set)
                else:
                    no_ref_index[i].append(index_set)
    # If boolean variables that are going to be reformulated are defined over a single set except:
    except:
        # index of the set where reformultion can be applied for a given boolean variable
        ref_index = {}
        # index of the sets where the reformulation cannot be applied for a given boolean variable
        no_ref_index = {}
        for i in ext_ref:
            ref_index[i] = []
            no_ref_index[i] = []
            if i.index_set().name == ext_ref[i].name:
                ref_index[i].append(0)
            else:
                no_ref_index[i].append(0)

    # Identify the variables that can be reformulated by performing a loop over logical constraints
    count = 1
    # dict of dicts: it contains information from the exactly variables that can be reformulated into external variables.
    reformulation_dict = {}
    for c in m.component_data_objects(pe.LogicalConstraint, descend_into=True):
        if c.body.getname() == 'exactly':
            exactly_number = c.body.args[0]
            for possible_Boolean in ext_ref:
                # expected boolean variable where the reformulation is going to be applied
                expected_Boolean = possible_Boolean.name
                Boolean_name_list = []
                Boolean_name_list = Boolean_name_list + [
                    c.body.args[1:][k]._component()._name
                    for k in range(len(c.body.args[1:]))
                ]
                if all(x == expected_Boolean for x in Boolean_name_list):
                    # expected ordered set index where the reformulation is going to be applied
                    expected_ordered_set_index = ref_index[possible_Boolean]
                    # index of sets where the reformulation is not applied
                    index_of_other_sets = no_ref_index[possible_Boolean]
                    if len(index_of_other_sets) >= 1:  # If there are other indexes
                        Other_Sets_listOFlists = []
                        verification_Other_Sets_listOFlists = []
                        for j in index_of_other_sets:
                            Other_Sets_listOFlists.append(
                                [
                                    c.body.args[1:][k].index()[j]
                                    for k in range(len(c.body.args[1:]))
                                ]
                            )
                            if all(
                                c.body.args[1:][x].index()[j]
                                == c.body.args[1:][0].index()[j]
                                for x in range(len(c.body.args[1:]))
                            ):
                                verification_Other_Sets_listOFlists.append(True)
                            else:
                                verification_Other_Sets_listOFlists.append(False)
                        # If we get to this point and it is true, it means that we can apply the reformulation for this combination of Boolean var and Exactly-type constraint
                        if all(verification_Other_Sets_listOFlists):
                            reformulation_dict[count] = {}
                            reformulation_dict[count]['exactly_number'] = exactly_number
                            # rearange boolean vars in constraint
                            sorted_args = sorted(
                                c.body.args[1:],
                                key=lambda x: x.index()[expected_ordered_set_index[0]],
                            )
                            # Now work with the ordered version sorted_args instead of c.body.args[1:]
                            reformulation_dict[count]['Boolean_vars_names'] = [
                                sorted_args[k].name for k in range(len(sorted_args))
                            ]
                            reformulation_dict[count]['Boolean_vars_ordered_index'] = [
                                sorted_args[k].index()[expected_ordered_set_index[0]]
                                for k in range(len(sorted_args))
                            ]
                            reformulation_dict[count]['Ext_var_lower_bound'] = 1
                            reformulation_dict[count]['Ext_var_upper_bound'] = len(
                                sorted_args
                            )

                            count = count + 1
                    # If there is only one index, then we can apply the reformulation at this point
                    else:
                        reformulation_dict[count] = {}
                        reformulation_dict[count]['exactly_number'] = exactly_number
                        # rearange boolean vars in constraint
                        sorted_args = sorted(c.body.args[1:], key=lambda x: x.index())
                        # Now work with the ordered version sorted_args instead of c.body.args[1:]
                        reformulation_dict[count]['Boolean_vars_names'] = [
                            sorted_args[k].name for k in range(len(sorted_args))
                        ]
                        reformulation_dict[count]['Boolean_vars_ordered_index'] = [
                            sorted_args[k].index() for k in range(len(sorted_args))
                        ]
                        reformulation_dict[count]['Ext_var_lower_bound'] = 1
                        reformulation_dict[count]['Ext_var_upper_bound'] = len(
                            sorted_args
                        )

                        count = count + 1

    number_of_external_variables = sum(
        reformulation_dict[j]['exactly_number'] for j in reformulation_dict
    )

    lower_bounds = {}
    upper_bounds = {}

    exvar_num = 1
    for i in reformulation_dict:
        for j in range(reformulation_dict[i]['exactly_number']):
            lower_bounds[exvar_num] = reformulation_dict[i]['Ext_var_lower_bound']
            upper_bounds[exvar_num] = reformulation_dict[i]['Ext_var_upper_bound']
        exvar_num = exvar_num + 1

    if tee:
        print(
            '\nReformulation Summary\n--------------------------------------------------------------------------'
        )
        exvar_num = 0
        for i in reformulation_dict:
            for j in range(reformulation_dict[i]['exactly_number']):
                print(
                    'External variable x['
                    + str(exvar_num)
                    + '] '
                    + ' is associated to '
                    + str(reformulation_dict[i]['Boolean_vars_names'])
                    + ' and it must be within '
                    + str(reformulation_dict[i]['Ext_var_lower_bound'])
                    + ' and '
                    + str(reformulation_dict[i]['Ext_var_upper_bound'])
                    + '.'
                )
                exvar_num = exvar_num + 1

        print(
            '\nThere are '
            + str(number_of_external_variables)
            + ' external variables in total'
        )

    return reformulation_dict, number_of_external_variables, lower_bounds, upper_bounds


def external_ref(
    m: pe.ConcreteModel(),
    x,
    extra_logic_function,
    dict_extvar: dict = {},
    mip_ref: bool = False,
    transformation: str = 'bigm',
    tee: bool = False,
):
    """
    Function that reformulates a given GDP model by taking into account external variables. This is particularly useful when translating GDP models to MINLP representations. The function allows the user to provide additional logical constraints and can optionally convert the model to its equivalent MINLP form using specified transformations.

    Args:
        m (pyomo.ConcrteteModel): GDP model that is going to be reformulated
        x (list): List with current value of the external variables
        extra_logic_function (function): Function that returns a list of lists of the form [a,b], where a is an expressions of the reformulated Boolean variables and b is an equivalent Boolean or indicator variable (b<->a)
        dict_extvar (dict): A dictionary of dictionaries that looks as follows:
            {1:{'exactly_number':Number of external variables for this type,
                'Boolean_vars_names':list with names of the ordered Boolean variables to be reformulated,
                'Boolean_vars_ordered_index': Indexes where the external reformulation is applied,
                'Binary_vars_names':list with names of the ordered Binary variables to be reformulated, [Potentially]
                'Binary_vars_ordered_index': Indexes where the external reformulation is applied, [Potentially]
                'Ext_var_lower_bound': Lower bound for this type of external variable,
                'Ext_var_upper_bound': Upper bound for this type of external variable },
             2:{...},...}

            The first key (positive integer) represent a type of external variable identified in the model. For this type of external variable
            a dictionary is created.
        mip_ref (bool): whether the reformulation will consider binary variables besides Booleans coming from a GDP->MIP reformulation
        transformation (str): GDP to MINLP transformation to be used
        tee (bool): Display reformulation
    Returns:
        m (pyomo.ConcreteModel): A model where the independent Boolean variables that were reformulated are fixed and Boolean/indicator variables that are calculated in
        terms of the independent Boolean variables are fixed too (depending on the extra_logic_function provided by the user)

    """
    # This part of code is required due to the deep copy issue: we have to compare Boolean variables by name
    for i in dict_extvar:
        dict_extvar[i]['Boolean_vars'] = []
        for j in dict_extvar[i]['Boolean_vars_names']:
            for boolean in m.component_data_objects(pe.BooleanVar, descend_into=True):
                if boolean.name == j:
                    dict_extvar[i]['Boolean_vars'] = dict_extvar[i]['Boolean_vars'] + [
                        boolean
                    ]
        if mip_ref:
            # This part of code is required due to the deep copy issue: we have to compare binary variables by name
            # By uncommenting in previous function extvars_gdp_to_mip we would pass directly dict_extvar[i]['Binary_vars']
            dict_extvar[i]['Binary_vars'] = []
            for j in dict_extvar[i]['Binary_vars_names']:
                for binary in m.component_data_objects(pe.Var, descend_into=True):
                    if binary.name == j:
                        dict_extvar[i]['Binary_vars'] = dict_extvar[i][
                            'Binary_vars'
                        ] + [binary]

    # The function would start here if there were no problems with deep copy.
    ext_var_position = 0
    for i in dict_extvar:
        for j in range(dict_extvar[i]['exactly_number']):
            for k in range(1, len(dict_extvar[i]['Boolean_vars']) + 1):
                if x[ext_var_position] == k:
                    if not mip_ref:
                        # fix True variables: depending on the current value of the external variables, some Independent Boolean variables can be fixed
                        dict_extvar[i]['Boolean_vars'][k - 1].fix(True)
                    else:
                        # fix 0 variables: depending on the current value of the external variables, some Independent Binary variables can be fixed
                        dict_extvar[i]['Binary_vars'][k - 1].fix(1)
                        dict_extvar[i]['Boolean_vars'][k - 1].set_value(True)
            ext_var_position = ext_var_position + 1
        # Double loop required from fact that exactly_number >= 1. TODO Is there a better way to do this?
        for j in range(dict_extvar[i]['exactly_number']):
            for k in range(1, len(dict_extvar[i]['Boolean_vars']) + 1):
                if not mip_ref:
                    # fix False variables: If the independent Boolean variable is not fixed at "True", then it is fixed at "False".
                    if not dict_extvar[i]['Boolean_vars'][k - 1].is_fixed():
                        dict_extvar[i]['Boolean_vars'][k - 1].fix(False)
                else:
                    # fix 0 variables: If the independent Boolean variable is not fixed at "1", then it is fixed at "0".
                    if not dict_extvar[i]['Binary_vars'][k - 1].is_fixed():
                        dict_extvar[i]['Binary_vars'][k - 1].fix(0)
                        dict_extvar[i]['Boolean_vars'][k - 1].set_value(False)

    # Other Boolean and Indicator variables are fixed depending on the information provided by the user
    logic_expr = extra_logic_function(m)
    for i in logic_expr:
        if not mip_ref:
            i[1].fix(pe.value(i[0]))
        else:
            i[1].set_value(pe.value(i[0]))

    pe.TransformationFactory('core.logical_to_linear').apply_to(m)
    if mip_ref:  # Transform problem to MINLP
        transformation_string = 'gdp.' + transformation
        pe.TransformationFactory(transformation_string).apply_to(m)
    else:  # Deactivate disjunction's constraints in the case of pure GDP
        pe.TransformationFactory('gdp.fix_disjuncts').apply_to(m)

    pe.TransformationFactory('contrib.deactivate_trivial_constraints').apply_to(
        m, tmp=False, ignore_infeasible=True
    )

    if tee:
        print('\nFixed variables at current iteration:\n')
        print('\n Independent Boolean variables\n')
        for i in dict_extvar:
            for k in range(1, len(dict_extvar[i]['Boolean_vars']) + 1):
                print(
                    dict_extvar[i]['Boolean_vars_names'][k - 1]
                    + '='
                    + str(dict_extvar[i]['Boolean_vars'][k - 1].value)
                )

        print('\n Dependent Boolean variables and disjunctions\n')
        for i in logic_expr:
            print(i[1].name + '=' + str(i[1].value))

        if mip_ref:
            print('\n Independent binary variables\n')
            for i in dict_extvar:
                for k in range(1, len(dict_extvar[i]['Binary_vars']) + 1):
                    print(
                        dict_extvar[i]['Binary_vars_names'][k - 1]
                        + '='
                        + str(dict_extvar[i]['Binary_vars'][k - 1].value)
                    )

    return m


def extvars_gdp_to_mip(
    m: pe.ConcreteModel(), gdp_dict_extvar: dict = {}, transformation: str = 'bigm'
):
    """
    Function that transforms a given Generalized Disjunctive Programming (GDP) model into its equivalent Mixed-Integer Nonlinear Programming (MINLP) representation. It additionally maps the external variables defined for the GDP model into corresponding binary variables suitable for the MINLP formulation. The function allows for different transformation techniques, depending on user preferences.

    Args:
        m (Pyomo.ConcreteModel): GDP model that is going to be reformulated
        gdp_dict_extvar (dict): A dictionary of dictionaries that looks as follows:
            {1:{'exactly_number':Number of external variables for this type,
                'Boolean_vars_names':list with names of the ordered Boolean variables to be reformulated,
                'Boolean_vars_ordered_index': Indexes where the external reformulation is applied,
                'Ext_var_lower_bound': Lower bound for this type of external variable,
                'Ext_var_upper_bound': Upper bound for this type of external variable },
             2:{...},...}
        transformation (str): GDP to MINLP transformation to be used

            The first key (positive integer) represent a type of external variable identified in the model. For this type of external variable
            a dictionary is created.
        tee (bool): Display reformulation
    Returns:
        m (Pyomo.ConcreteModel): A MIP model transformed from the original GDP m model via the 'transformation' argument
        mip_dict_extvar (dict): A dictionary of dictionaries that looks as follows:
            {1:{'exactly_number':Number of external variables for this type,
                'Boolean_vars_names':list with names of the ordered Boolean variables to be reformulated,
                'Boolean_vars_ordered_index': Indexes where the external reformulation is applied,
                'Binary_vars_names':list with names of the ordered Binary variables to be reformulated,
                'Binary_vars_ordered_index': Indexes where the external reformulation is applied,
                'Ext_var_lower_bound': Lower bound for this type of external variable,
                'Ext_var_upper_bound': Upper bound for this type of external variable },
             2:{...},...}

    """

    # Transformation step
    pe.TransformationFactory('core.logical_to_linear').apply_to(m)
    transformation_string = 'gdp.' + transformation
    pe.TransformationFactory(transformation_string).apply_to(m)

    mip_dict_extvar = copy.deepcopy(gdp_dict_extvar)

    # This part of code is required due to the deep copy issue: we have to compare Boolean variables by name
    for i in mip_dict_extvar.keys():
        mip_dict_extvar[i]['Boolean_vars'] = []
        for j in mip_dict_extvar[i]['Boolean_vars_names']:
            for boolean in m.component_data_objects(pe.BooleanVar, descend_into=True):
                if boolean.name == j:
                    mip_dict_extvar[i]['Boolean_vars'] = mip_dict_extvar[i][
                        'Boolean_vars'
                    ] + [boolean]
        # Add extra terms to the dictionary to be relevant for binary variables
        mip_dict_extvar[i]['Binary_vars_names'] = [
            boolean.get_associated_binary().name
            for boolean in mip_dict_extvar[i]['Boolean_vars']
        ]
        # Uncomment the next line in case that deepcopy works
        # mip_dict_extvar[key]['Binary_vars'] = [
        #     boolean.get_associated_binary() for boolean in gdp_dict_extvar[key]['Boolean_vars']]
        mip_dict_extvar[i]['Binary_vars_ordered_index'] = mip_dict_extvar[i][
            'Boolean_vars_ordered_index'
        ]

    return m, mip_dict_extvar


def preprocess_problem(m, simple: bool = True):
    """
    Function that applies certain tranformations to the model to first verify that it is not trivially
    infeasible (via FBBT) and second, remove extra constraints to help NLP solvers
    Args:
        m (pyomo.ConcreteModel): (MI)NLP model that is going to be preprocessed
        simple (bool): Boolean variable to carry on a simple preprocessing (only FBBT) or a more complete one, prone to fail
    Returns:

    """
    if not simple:
        pe.TransformationFactory('contrib.detect_fixed_vars').apply_to(m)
        pe.TransformationFactory('contrib.propagate_fixed_vars').apply_to(m)
        pe.TransformationFactory('contrib.remove_zero_terms').apply_to(m)
        pe.TransformationFactory('contrib.propagate_zero_sum').apply_to(m)
        pe.TransformationFactory('contrib.constraints_to_var_bounds').apply_to(m)
        pe.TransformationFactory('contrib.detect_fixed_vars').apply_to(m)
        pe.TransformationFactory('contrib.propagate_zero_sum').apply_to(m)
        pe.TransformationFactory('contrib.deactivate_trivial_constraints').apply_to(
            m, tmp=False, ignore_infeasible=True
        )
    fbbt(m)


def solve_subproblem(
    m: pe.ConcreteModel(),
    subproblem_solver: str = 'knitro',
    subproblem_solver_options: dict = {},
    timelimit: float = 10,
    gams_output: bool = False,
    tee: bool = False,
    rel_tol: float = 1e-3,
) -> pe.ConcreteModel():
    """
    Function that checks feasibility and subproblem model.
    Note integer variables have to be previously fixed in the external reformulation
    Args:
        m (pyomo.ConcreteModel): Fixed subproblem model that is to be solved
        subproblem_solver (str): MINLP or NLP solver algorithm
        subproblem_solver_options (dict): MINLP or NLP solver algorithm options
        timelimit (float): time limit in seconds for the solve statement
        gams_output (bool): Determine keeping or not GAMS files
        tee (bool): Display iteration output
        rel_tol (float): Relative optimality tolerance
    Returns:
        m (pyomo.ConcreteModel): Solved subproblem model
    """
    # Initialize D-SDA status
    m.dsda_status = 'Initialized'
    m.dsda_usertime = 0

    try:
        # Feasibility and preprocessing checks
        preprocess_problem(m, simple=True)

    except InfeasibleConstraintException:
        m.dsda_status = 'FBBT_Infeasible'
        return m

    output_options = {}

    # Output report
    if gams_output:
        dir_path = os.path.dirname(os.path.abspath(__file__))
        gams_path = os.path.join(dir_path, "gamsfiles/")
        if not (os.path.exists(gams_path)):
            print(
                'Directory for automatically generated files '
                + gams_path
                + ' does not exist. We will create it'
            )
            os.makedirs(gams_path)
        output_options = {
            'keepfiles': True,
            'tmpdir': gams_path,
            'symbolic_solver_labels': True,
        }

    subproblem_solver_options['add_options'] = subproblem_solver_options.get(
        'add_options', []
    )
    subproblem_solver_options['add_options'].append('option reslim=%s;' % timelimit)
    subproblem_solver_options['add_options'].append('option optcr=%s;' % rel_tol)

    # Solve
    solvername = 'gams'
    opt = SolverFactory(solvername, solver=subproblem_solver)
    m.results = opt.solve(
        m,
        tee=tee,
        **output_options,
        **subproblem_solver_options,
        skip_trivial_constraints=True,
    )

    m.dsda_usertime = m.results.solver.user_time

    # Assign D-SDA status
    if m.results.solver.termination_condition == 'infeasible':
        m.dsda_status = 'Evaluated_Infeasible'
    else:  # Considering locallyOptimal, optimal, globallyOptimal, and maxtime TODO: Fix this
        m.dsda_status = 'Optimal'
    # if m.results.solver.termination_condition == 'locallyOptimal' or m.results.solver.termination_condition == 'optimal' or m.results.solver.termination_condition == 'globallyOptimal':
    #     m.dsda_status = 'Optimal'

    return m


def solve_with_minlp(
    m: pe.ConcreteModel(),
    transformation: str = 'bigm',
    minlp: str = 'baron',
    minlp_options: dict = {},
    timelimit: float = 10,
    gams_output: bool = False,
    tee: bool = False,
    rel_tol: float = 0.001,
) -> pe.ConcreteModel():
    """
    Function that transforms a GDP model and solves it as a mixed-integer nonlinear
    programming (MINLP) model.
    Args:
        m (pyomo.ConcreteModel): Pyomo GDP model that is to be solved using MINLP
        transformation (str): GDP to MINLP transformation to be used
        minlp (str): MINLP solver algorithm
        minlp_options (dict): MINLP solver algorithm options
        timelimit (float): time limit in seconds for the solve statement
        gams_output (bool): Determine keeping or not GAMS files
        tee (bool): Display iterations
        rel_tol (float): Relative optimality tolerance
    Returns:
        m (pyomo.ConcreteModel): Solved MINLP model
    """

    # Transformation step from GDP to MINLP
    pe.TransformationFactory('core.logical_to_linear').apply_to(m)
    transformation_string = 'gdp.' + transformation
    pe.TransformationFactory(transformation_string).apply_to(m)

    # Output report
    output_options = {}
    if gams_output:
        # Find the directory where this file is located
        dir_path = os.path.dirname(os.path.abspath(__file__))
        # Create a directory for the GAMS files
        gams_path = os.path.join(dir_path, "gamsfiles/")
        # If the directory does not exist, create it
        if not (os.path.exists(gams_path)):
            print(
                'Directory for automatically generated files '
                + gams_path
                + ' does not exist. We will create it'
            )
            os.makedirs(gams_path)
        # Set the output options for GAMS to create the output files with naming convention
        output_options = {
            'keepfiles': True,
            'tmpdir': gams_path,
            'symbolic_solver_labels': True,
        }

    # Set the options for the MINLP solver
    # Setting options from the input
    minlp_options['add_options'] = minlp_options.get('add_options', [])
    # Setting the time limit
    minlp_options['add_options'].append('option reslim=%s;' % timelimit)
    # Setting the relative optimality tolerance
    minlp_options['add_options'].append('option optcr=%s;' % rel_tol)

    # Solve
    solvername = 'gams'
    opt = SolverFactory(solvername, solver=minlp)
    m.results = opt.solve(m, tee=tee, **output_options, **minlp_options)
    # update_boolean_vars_from_binary(m)
    return m


def solve_with_gdpopt(
    m: pe.ConcreteModel(),
    mip: str = 'cplex',
    mip_options: dict = {},
    nlp: str = 'knitro',
    nlp_options: dict = {},
    minlp: str = 'baron',
    minlp_options: dict = {},
    timelimit: float = 10,
    strategy: str = 'LOA',
    mip_output: bool = False,
    nlp_output: bool = False,
    minlp_output: bool = False,
    tee: bool = False,
    rel_tol: float = 1e-3,
) -> pe.ConcreteModel():
    """
    Function that solves GDP model using GDPopt
    Args:
        m (pyomo.ConcreteModel): GDP model that is to be solved
        mip (str): MIP solver algorithm
        mip_options (dict): MIP solver algorithm options
        nlp (str): NLP solver algorithm
        nlp_options (dict): NLP solver algorithm options
        minlp (str): MINLP solver algorithm
        minlp_options (dict): MINLP solver algorithm options
        timelimit (float): time limit in seconds for the solve statement
        strategy (str): GDPopt strategy
        mip_output (bool): Determine keeping or not GAMS files of the MIP model
        nlp_output (bool): Determine keeping or not GAMS files of the NLP model
        minlp_output (bool): Determine keeping or not GAMS files of the MINLP model
        tee (bool): Display iterations
        rel_tol (float): Relative optimality tolerance for subproblems and GDPOpt itself
    Returns:
        m (pyomo.ConcreteModel): Solved GDP model
    """

    # Transformation step
    pe.TransformationFactory('core.logical_to_linear').apply_to(m)

    # Output report

    mip_options['add_options'] = mip_options.get('add_options', [])
    mip_options['add_options'].append('option optcr=0.0;')

    nlp_options['add_options'] = nlp_options.get('add_options', [])
    nlp_options['add_options'].append('option optcr=%s;' % rel_tol)

    minlp_options['add_options'] = minlp_options.get('add_options', [])
    minlp_options['add_options'].append('option optcr=%s;' % rel_tol)

    if mip_output:
        dir_path = os.path.dirname(os.path.abspath(__file__))
        gams_path = os.path.join(dir_path, "gamsfiles/")
        if not (os.path.exists(gams_path)):
            print(
                'Directory for automatically generated files '
                + gams_path
                + ' does not exist. We will create it'
            )
            os.makedirs(gams_path)
        mip_options['keepfiles'] = True
        mip_options['tmpdir'] = gams_path
        mip_options['symbolic_solver_labels'] = True

    if nlp_output:
        dir_path = os.path.dirname(os.path.abspath(__file__))
        gams_path = os.path.join(dir_path, "gamsfiles/")
        if not (os.path.exists(gams_path)):
            print(
                'Directory for automatically generated files '
                + gams_path
                + ' does not exist. We will create it'
            )
            os.makedirs(gams_path)
        nlp_options['keepfiles'] = True
        nlp_options['tmpdir'] = gams_path
        nlp_options['symbolic_solver_labels'] = True

    if minlp_output:
        dir_path = os.path.dirname(os.path.abspath(__file__))
        gams_path = os.path.join(dir_path, "gamsfiles/")
        if not (os.path.exists(gams_path)):
            print(
                'Directory for automatically generated files '
                + gams_path
                + ' does not exist. We will create it'
            )
            os.makedirs(gams_path)
        minlp_options['keepfiles'] = True
        minlp_options['tmpdir'] = gams_path
        minlp_options['symbolic_solver_labels'] = True

    # Solve
    solvername = 'gdpopt'
    opt = SolverFactory(solvername)
    m.results = opt.solve(
        m,
        tee=tee,
        strategy=strategy,
        time_limit=timelimit,
        mip_solver='gams',
        mip_solver_args=dict(solver=mip, warmstart=True, **mip_options),
        nlp_solver='gams',
        nlp_solver_args=dict(solver=nlp, warmstart=True, tee=tee, **nlp_options),
        minlp_solver='gams',
        minlp_solver_args=dict(solver=minlp, warmstart=True, tee=tee, **minlp_options),
        #   mip_presolve=True,
        init_strategy='fix_disjuncts',
        #   set_cover_iterlim=0,
        iterlim=1000,
        force_subproblem_nlp=True,
        subproblem_presolve=False
        #   bound_tolerance=rel_tol
        #   calc_disjunctive_bounds=True
    )
    # update_boolean_vars_from_binary(m)
    return m


def neighborhood_k_eq_2(dimension: int = 2) -> dict:
    """
    Function creates a k=2 neighborhood of the given dimension
    Args:
        dimension (int): Dimension of the neighborhood
    Returns:
        directions (dict): Dictionary contaning in each item a list with a direction within the neighborhood
    """

    num_neigh = 2 * dimension
    neighbors = np.concatenate(
        (np.eye(dimension, dtype=int), -np.eye(dimension, dtype=int)), axis=1
    )
    directions = {}
    for i in range(num_neigh):
        direct = []
        directions[i + 1] = direct
        for j in range(dimension):
            direct.append(neighbors[j, i])
    return directions


def neighborhood_k_eq_inf(dimension: int = 2) -> dict:
    """
    Function creates a k=Infinity neighborhood of the given dimension
    Args:
        dimension (int): Dimension of the neighborhood
    Returns:
        temp (dict): Dictionary contaning in each item a list with a direction within the neighborhood
    Note:
        temp is temporary variable.
    """

    neighbors = list(it.product([-1, 0, 1], repeat=dimension))
    directions = {}
    for i in range(len(neighbors)):
        directions[i + 1] = list(neighbors[i])
    temp = directions.copy()
    for i in directions.keys():
        if temp[i] == [0] * dimension:
            temp.pop(i, None)
    return temp


def initialize_model(
    m: pe.ConcreteModel(),
    json_path=None,
    from_feasible: bool = False,
    feasible_model: str = '',
) -> pe.ConcreteModel():
    """
    Function that return an initialized model from an existing json file
    Args:
        m (pyomo.ConcreteModel): Pyomo model that is to be initialized
        from_feasible (bool): If initialization is made from an external file
        feasible_model (str): Feasible initialization path or example
    Returns:
        m (pyomo.ConcreteModel): Initialized Pyomo model
    """

    wts = StoreSpec.value()

    if json_path is None:
        os.path.join(os.path.curdir)

        dir_path = os.path.dirname(os.path.abspath(__file__))

        if from_feasible:
            json_path = os.path.join(dir_path, feasible_model + '_initialization.json')
        else:
            json_path = os.path.join(dir_path, 'dsda_initialization.json')

    from_json(m, fname=json_path, wts=wts)
    return m


def generate_initialization(
    m: pe.ConcreteModel(),
    starting_initialization: bool = False,
    model_name: str = '',
    human_read: bool = True,
    wts=StoreSpec.value(),
):
    """
    Function that creates a json file for initialization based on a model m
    Args:
        m (pyomo.ConcreteModel): Base Pyomo model for initializtion
        starting_intialization (bool): Use to create "dsda_starting_initialization.json" file with a known feasible initialized model m
        model_name (str): Name of the model for the initialization
        human_read (bool): Make the json file readable by a human
        wts (StoreSpec.value): What to save, initially the values, but we might want something different. Check model_serializer tests for examples
    Returns:
        json_path (os.dir): Path where json file is stored
    """

    dir_path = os.path.dirname(os.path.abspath(__file__))

    if starting_initialization:
        json_path = os.path.join(dir_path, model_name + '_initialization.json')
    else:
        if model_name != '':
            json_path = os.path.join(dir_path, model_name + '_initialization.json')
        else:
            json_path = os.path.join(dir_path, 'dsda_initialization.json')

    to_json(m, fname=json_path, human_read=human_read, wts=wts)

    return json_path


def find_actual_neighbors(
    start: list, neighborhood: dict, min_allowed: dict = {}, max_allowed: dict = {}
) -> dict:
    """
    Function that creates all neighbors of a given point. Neighbor 0 is the starting point
    Args:
        start (list): Point of which neighbors want to be created
        neighborhood (dict): Neighborhood (output of a k-Neighborhood function)
        min_allowed (dict): In keys contains external variables and in items their respective lower bounds
        max_allowed (dict): In keys contains external variables and in items their respective upper bounds
    Returns:
        new_neighbors (dict): Contains neighbors of the actual point
    """

    neighbors = {0: start}
    for i in neighborhood.keys():  # Calculate neighbors
        neighbors[i] = list(map(sum, zip(start, list(neighborhood[i]))))

    new_neighbors = {}
    num_vars = len(neighbors[0])
    for i in neighbors.keys():
        checked = 0
        for j in range(num_vars):  # Check if within bounds
            if (
                neighbors[i][j] >= min_allowed[j + 1]
                and neighbors[i][j] <= max_allowed[j + 1]
            ):
                checked += 1
        if checked == num_vars:  # Add neighbor if all variables are within bounds
            new_neighbors[i] = neighbors[i]

    return new_neighbors


def evaluate_neighbors(
    ext_vars: dict,
    fmin: float,
    model_function,
    model_args: dict,
    ext_dict: dict,
    ext_logic,
    mip_transformation: bool = False,
    transformation: str = 'bigm',
    subproblem_solver: str = 'knitro',
    subproblem_solver_options: dict = {},
    iter_timelimit: float = 10,
    current_time: float = 0,
    timelimit: float = 3600,
    gams_output: bool = False,
    tee: bool = False,
    global_tee: bool = True,
    rel_tol: float = 1e-3,
    global_evaluated: list = [],
    init_path=None,
):
    """
    Function that evaluates a group of given points and returns the best
    Args:
        ext_vars (dict): dict with neighbors where neighbor 0 is actual point
        fmin (float): Objective at actual point
        model_function (function): function that returns GDP model to be solved
        model_args (dict): Contains the argument values needed for model_function
        ext_dict (dict): Dictionary with Boolean variables to be reformulated (keys) and their corresponding ordered sets (values)
        ext_logic (function): Function that returns a list of lists of the form [a,b], where a is an expressions of the reformulated Boolean variables and b is an equivalent Boolean or indicator variable (b<->a)
        mip_transformation (bool): Whether to solve the enumeration using the external variables applied to the MIP problem insed of the GDP
        transformation (str): Which transformation to apply to the GDP
        subproblem_solver (str): MINLP or NLP solver algorithm
        subproblem_solver_options (dict): MINLP or NLP solver algorithm options
        iter_timelimit (float): time limit in seconds for the solve statement for each iteration
        current_time (float): Current time in global algorithm
        timelimit (float): time limit in seconds for the algorithm
        gams_output (bool): Determine keeping or not GAMS files
        tee (bool): Display iteration output
        global_tee (bool): display D-SDA iteration output
        rel_tol (float): Relative optimality tolerance
        global_evaluated (list): list with points already evaluated
        init_path (os.path): path to initialization file

    Returns:
        fmin (float): Gives the best neighbor's objective
        best_var (list): Type list and gives the best neighbor
        best_dir (int): Type int and is the steepest direction (key in neighborhood)
        improve (bool): Type bool and shows if an improvement was made while looking for neighbors
        evaluation_time (float): Total solver-statement time only
        ns_evaluated (list) : evaluations in neighbor search
        best_path (os.path): path to json with best solution found

    """
    # Global Tolerance parameters
    epsilon = 1e-10
    abs_tol = 1e-5
    min_improve = 1e-5
    min_improve_rel = 1e-3

    # Initialize
    ns_evaluated = []
    evaluation_time = 0
    improve = False
    best_var = ext_vars[0]
    here = ext_vars[0]
    best_dir = 0  # Position in dictionary
    best_dist = 0
    best_path = init_path
    temp = ext_vars  # TODO change name to something more saying. Points? Combinations?
    temp.pop(0, None)

    if global_tee:
        print()
        print('Neighbor search around:', best_var)

    for i in temp.keys():  # Solve all models
        if temp[i] not in global_evaluated:
            m = model_function(**model_args)
            m_init = initialize_model(m, json_path=init_path)
            if (
                mip_transformation
            ):  # If you want a MIP reformulation, go ahead and use it'
                m_init, ext_dict = extvars_gdp_to_mip(
                    m=m, gdp_dict_extvar=ext_dict, transformation=transformation
                )

            m_fixed = external_ref(
                m=m_init,
                x=temp[i],
                extra_logic_function=ext_logic,
                dict_extvar=ext_dict,
                mip_ref=mip_transformation,
                tee=False,
            )
            t_remaining = min(
                iter_timelimit, timelimit - (time.perf_counter() - current_time)
            )
            if t_remaining < 0:  # No time reamining for optimization
                break
            m_solved = solve_subproblem(
                m=m_fixed,
                subproblem_solver=subproblem_solver,
                subproblem_solver_options=subproblem_solver_options,
                timelimit=t_remaining,
                gams_output=gams_output,
                tee=tee,
            )
            evaluation_time += m_solved.dsda_usertime
            ns_evaluated.append(temp[i])
            t_end = time.perf_counter()

            if m_solved.dsda_status == 'Optimal':  # Check if D-SDA status is optimal
                if global_tee:
                    print(
                        'Evaluated:',
                        temp[i],
                        '   |   Objective:',
                        round(pe.value(m_solved.obj), 5),
                        '   |   Global Time:',
                        round(t_end - current_time, 2),
                    )
                dist = sum((x - y) ** 2 for x, y in zip(temp[i], here))
                act_obj = pe.value(m_solved.obj)

                # Assuming minimization problem
                # Implements heuristic of largest move
                if not improve:
                    # We want a minimum improvement in the first found solution
                    if (fmin - act_obj) > min_improve or (fmin - act_obj) / (
                        abs(fmin) + epsilon
                    ) > min_improve_rel:
                        fmin = act_obj
                        best_var = temp[i]
                        best_dir = i
                        best_dist = dist
                        improve = True
                        best_path = generate_initialization(
                            m_solved, starting_initialization=False, model_name='best'
                        )
                else:
                    # We want slightly worse solutions if the distance is larger
                    if (
                        ((act_obj - fmin) < abs_tol)
                        or ((act_obj - fmin) / (abs(fmin) + epsilon) < rel_tol)
                    ) and dist >= best_dist:
                        fmin = act_obj
                        best_var = temp[i]
                        best_dir = i
                        best_dist = dist
                        improve = True
                        best_path = generate_initialization(
                            m_solved, starting_initialization=False, model_name='best'
                        )

            if time.perf_counter() - current_time > timelimit:  # current
                break

    if global_tee:
        print()
        print('New best neighbor:', best_var)
    return fmin, best_var, best_dir, improve, evaluation_time, ns_evaluated, best_path


def do_line_search(
    start: list,
    fmin: float,
    direction: list,
    model_function,
    model_args: dict,
    ext_dict: dict,
    ext_logic,
    mip_transformation: bool = False,
    transformation: str = 'bigm',
    subproblem_solver: str = 'knitro',
    subproblem_solver_options: dict = {},
    min_allowed: dict = {},
    max_allowed: dict = {},
    iter_timelimit: float = 10,
    timelimit: float = 3600,
    current_time: float = 0,
    gams_output: bool = False,
    tee: bool = False,
    global_tee: bool = False,
    rel_tol: float = 1e-3,
    global_evaluated: list = [],
    init_path=None,
):
    """
    Function that moves in a given "best direction" and evaluates the new moved point

    Args:
        start (list): Point of that is to be moved
        fmin (float): Objective at actual point
        direction: moving direction
        model_function (function): function that returns GDP model to be solved
        model_args (dict): Contains the argument values needed for model_function
        ext_dict (dict): Dictionary with Boolean variables to be reformulated (keys) and their corresponding ordered sets (values)
        ext_logic (function): Function that returns a list of lists of the form [a,b], where a is an expressions of the reformulated Boolean variables and b is an equivalent Boolean or indicator variable (b<->a)
        mip_transformation (bool): Whether to solve the enumeration using the external variables applied to the MIP problem insed of the GDP
        transformation (str): Which transformation to apply to the GDP
        subproblem_solver (str): MINLP or NLP solver algorithm
        subproblem_solver_options (dict): MINLP or NLP solver algorithm options
        min_allowed (dict): In keys contains external variables and in items their respective lower bounds
        max_allowed (dict): In keys contains external variables and in items their respective upper bounds
        iter_timelimit (float): time limit in seconds for the solve statement for each iteration
        current_time (float): Current time in global algorithm
        gams_output (bool): Determine keeping or not GAMS files
        tee (bool): Display iteration output
        global_tee (gool): display D-SDA iteration output
        rel_tol (float): Relative optimality tolerance
        global_evaluated (list): list with points already evaluated
        init_path (os.path): path to initialization file
    Returns:
        fmin (int): Type int and gives the moved point objective
        best_var (list): Type list and gives the moved point
        moved (bool): Type bool and shows if an improvement was made while line searching
        ls_time (float): Total solver-statement time only
        ls_evaluated (list): evaluations in line search
        new_path (os.path): path of best json file
    """
    # Global Tolerance parameters
    epsilon = 1e-10
    min_improve = 1e-5
    min_improve_rel = 1e-3

    # Initialize
    ls_evaluated = []
    ls_time = 0
    best_var = start
    moved = False
    new_path = init_path

    # Line search in given direction
    moved_point = list(map(sum, zip(list(start), list(direction))))
    checked = 0
    for j in range(len(moved_point)):  # Check if within bounds
        if (
            moved_point[j] >= min_allowed[j + 1]
            and moved_point[j] <= max_allowed[j + 1]
        ):
            checked += 1

    if checked == len(moved_point):  # Solve model
        if moved_point not in global_evaluated:
            m = model_function(**model_args)
            m_init = initialize_model(m, json_path=init_path)
            if (
                mip_transformation
            ):  # If you want a MIP reformulation, go ahead and use it'
                m_init, ext_dict = extvars_gdp_to_mip(
                    m=m, gdp_dict_extvar=ext_dict, transformation=transformation
                )
            m_fixed = external_ref(
                m=m_init,
                x=moved_point,
                extra_logic_function=ext_logic,
                dict_extvar=ext_dict,
                mip_ref=mip_transformation,
                tee=False,
            )

            t_remaining = min(
                iter_timelimit, timelimit - (time.perf_counter() - current_time)
            )
            if t_remaining < 0:
                return fmin, best_var, moved, ls_time, ls_evaluated, new_path
            m_solved = solve_subproblem(
                m=m_fixed,
                subproblem_solver=subproblem_solver,
                subproblem_solver_options=subproblem_solver_options,
                timelimit=t_remaining,
                gams_output=gams_output,
                tee=tee,
            )
            ls_time += m_solved.dsda_usertime
            ls_evaluated.append(moved_point)

            if m_solved.dsda_status == 'Optimal':  # Check status
                if global_tee:
                    print(
                        'Evaluated:',
                        moved_point,
                        '   |   Objective:',
                        round(pe.value(m_solved.obj), 5),
                        '   |   Global Time:',
                        round(time.perf_counter() - current_time, 2),
                    )
                act_obj = pe.value(m_solved.obj)
                # Return moved point
                if (fmin - act_obj) > min_improve or (fmin - act_obj) / (
                    abs(fmin) + epsilon
                ) > min_improve_rel:
                    fmin = act_obj
                    best_var = moved_point
                    moved = True
                    new_path = generate_initialization(
                        m_solved, starting_initialization=False, model_name='best'
                    )

    return fmin, best_var, moved, ls_time, ls_evaluated, new_path


def solve_with_dsda(
    model_function,
    model_args: dict,
    starting_point: list,
    ext_dict,
    ext_logic,
    mip_transformation: bool = False,
    transformation: str = 'bigm',
    k: str = 'Infinity',
    provide_starting_initialization: bool = True,
    feasible_model: str = '',
    subproblem_solver: str = 'knitro',
    subproblem_solver_options: dict = {},
    iter_timelimit: float = 10,
    timelimit: float = 3600,
    gams_output: bool = False,
    tee: bool = False,
    global_tee: bool = True,
    rel_tol: float = 1e-3,
):
    """
    Function that computes Discrete-Steepest Descend Algorithm
    Args:
        model_function (function): function that returns GDP model to be solved
        model_args (dict): Contains the argument values needed for model_function
        starting_point (list): Feasible external variable initial point
        ext_dict (dict): Dictionary with Boolean variables to be reformulated (keys) and their corresponding ordered sets (values). Both keys and values are pyomo objects.
        ext_logic (function): Function that returns a list of lists of the form [a,b], where a is an expressions of the reformulated Boolean variables and b is an equivalent Boolean or indicator variable (b<->a).
        mip_transformation (bool): Whether to solve the enumeration using the external variables applied to the MIP problem instead of the GDP
        transformation (str): Which transformation to apply to the GDP
        k (string): Type of neighborhood ('2' or 'Infinity)
        provide_intialization (bool): If an existing json file is provided with a feasible initialization of starting_point
        subproblem_solver (str): MINLP or NLP solver algorithm
        subproblem_solver_options (dict): MINLP or NLP solver algorithm options
        iter_timelimi (float): time limit in seconds for the solve statement for each iteration
        timelimit (float): time limit in seconds for the algorithm
        gams_output (bool): Determine keeping or not GAMS files
        tee (bool): Display iteration output
        global_tee (bool): Display D-SDA output
        rel_tol (float): Relative optimality tolerance
    Returns:
        m2_solved (pyomo.ConcreteModel): Solved Pyomo Model
        route (list): List containing points evaluated in throughout iteration
        obj_route (list): List containing objectives evaluated in throughout iteration

    """

    if global_tee:
        print('\nStarting D-SDA with k =', k)
        print(
            '--------------------------------------------------------------------------'
        )

    # Initialize
    route = []
    obj_route = []
    global_evaluated = []
    ext_var = starting_point

    # Check if  feasible initialization is provided
    m = model_function(**model_args)
    dict_extvar, num_ext_var, min_allowed, max_allowed = get_external_information(
        m, ext_dict
    )
    if len(starting_point) != num_ext_var:
        print(
            "The size of the initialization vector must be equal to " + str(num_ext_var)
        )

    t_start = time.perf_counter()
    dsda_usertime = 0
    if provide_starting_initialization:
        m_init = initialize_model(
            m, from_feasible=True, feasible_model=feasible_model, json_path=None
        )
    else:
        m_init = m

    if mip_transformation:  # If you want a MIP reformulation, go ahead and use it'
        m_init, dict_extvar = extvars_gdp_to_mip(
            m=m, gdp_dict_extvar=dict_extvar, transformation=transformation
        )

    m_fixed = external_ref(
        m=m_init,
        x=ext_var,
        extra_logic_function=ext_logic,
        dict_extvar=dict_extvar,
        mip_ref=mip_transformation,
        tee=False,
    )

    # Solve for initialization
    m_solved = solve_subproblem(
        m=m_fixed,
        subproblem_solver=subproblem_solver,
        subproblem_solver_options=subproblem_solver_options,
        timelimit=iter_timelimit,
        gams_output=gams_output,
        tee=tee,
    )
    dsda_usertime += m_solved.dsda_usertime
    fmin = pe.value(m_solved.obj)
    if global_tee:
        print('Initializing...')
        print(
            'Evaluated:',
            ext_var,
            '   |   Objective:',
            round(fmin, 5),
            '   |   Global Time:',
            round(time.perf_counter() - t_start, 2),
        )

    # m_solved.pprint()
    best_path = generate_initialization(m_solved)

    route.append(ext_var)
    obj_route.append(fmin)
    global_evaluated.append(ext_var)

    # Define neighborhood
    if k == '2':
        neighborhood = neighborhood_k_eq_2(len(ext_var))
    elif k == 'Infinity':
        neighborhood = neighborhood_k_eq_inf(len(ext_var))
    else:
        return "Enter a valid neighborhood ('Infinity' or '2')"

    looking_in_neighbors = True

    # Look in neighbors (outer cycle)
    while looking_in_neighbors:
        if time.perf_counter() - t_start > timelimit:
            break

        # Find neighbors of the actual point
        neighbors = find_actual_neighbors(
            ext_var, neighborhood, min_allowed=min_allowed, max_allowed=max_allowed
        )

        if time.perf_counter() - t_start > timelimit:
            break

        (
            fmin,
            best_var,
            best_dir,
            improve,
            eval_time,
            ns_evaluated,
            best_path,
        ) = evaluate_neighbors(
            ext_vars=neighbors,
            fmin=fmin,
            model_function=model_function,
            model_args=model_args,
            ext_dict=dict_extvar,
            ext_logic=ext_logic,
            mip_transformation=mip_transformation,
            transformation=transformation,
            subproblem_solver=subproblem_solver,
            subproblem_solver_options=subproblem_solver_options,
            iter_timelimit=iter_timelimit,
            timelimit=timelimit,
            current_time=t_start,
            gams_output=gams_output,
            tee=tee,
            global_tee=global_tee,
            rel_tol=rel_tol,
            global_evaluated=global_evaluated,
            init_path=best_path,
        )

        dsda_usertime += eval_time
        global_evaluated = global_evaluated + ns_evaluated

        # Stopping condition in case there is no improvement amongst neighbors
        if improve:
            line_searching = True
            route.append(best_var)
            obj_route.append(fmin)
            if global_tee and time.perf_counter() - t_start < timelimit:
                print()
                print('Line search in direction:', neighborhood[best_dir])

            # If improvement was made start line search (inner cycle)
            while line_searching:
                if time.perf_counter() - t_start > timelimit:
                    break

                (
                    fmin,
                    best_var,
                    moved,
                    ls_time,
                    ls_evaluated,
                    best_path,
                ) = do_line_search(
                    start=best_var,
                    fmin=fmin,
                    direction=neighborhood[best_dir],
                    model_function=model_function,
                    model_args=model_args,
                    ext_dict=dict_extvar,
                    ext_logic=ext_logic,
                    mip_transformation=mip_transformation,
                    transformation=transformation,
                    subproblem_solver=subproblem_solver,
                    min_allowed=min_allowed,
                    max_allowed=max_allowed,
                    iter_timelimit=iter_timelimit,
                    timelimit=timelimit,
                    current_time=t_start,
                    gams_output=gams_output,
                    tee=tee,
                    global_tee=global_tee,
                    rel_tol=rel_tol,
                    global_evaluated=global_evaluated,
                    init_path=best_path,
                )
                global_evaluated = global_evaluated + ls_evaluated
                dsda_usertime += ls_time

                if time.perf_counter() - t_start > timelimit:
                    break

                # Stopping condition in case no movement was done
                if moved:
                    route.append(best_var)
                    obj_route.append(fmin)
                else:
                    ext_var = best_var
                    line_searching = False
                    if global_tee:
                        print()
                        print('New best point:', best_var)

        else:
            looking_in_neighbors = False

    t_end = round(time.perf_counter() - t_start, 2)

    # Generate final solved model
    m2 = model_function(**model_args)
    m2_solved = initialize_model(m2, json_path=best_path)
    m2_solved.dsda_time = t_end
    m2_solved.dsda_usertime = dsda_usertime
    if t_end > timelimit:
        m2_solved.dsda_status = 'maxTimeLimit'
    else:
        m2_solved.dsda_status = 'optimal'

    # Print results
    if global_tee:
        print(
            '--------------------------------------------------------------------------'
        )
        print('Objective:', round(fmin, 5))
        print('External variables:', route[-1])
        print('Execution time [s]:', t_end)
        print('User time [s]:', round(dsda_usertime, 5))

    return m2_solved, route, obj_route


def visualize_dsda(
    route: list = [],
    feas_x: list = [],
    feas_y: list = [],
    objs: list = [],
    k: str = '?',
    ext1_name: str = 'External variable 1',
    ext2_name: str = 'External variable 2',
):
    """
    Function that plots Discrete-Steepest Descend Algorithm for two external variables
    Args:
        route (list): List containing points evaluated in throughout iteration
        feas_x (list): List containing x-axis position of feasible points
        feas_y (list): List containing y-axis position of feasible points
        objs (list): List containing objective function of feasible points
        k (str): Type of neighborhood
        ext1_name (str): External variable 1 name
        ext2_name (str): External variable 2 name

    Returns:
        None: Directly visualizes the DSDA progression
    """

    # Plot feasible points
    # Rename feas_x and feas_y tp X1 and X2 (to possibly filter)
    X1, X2 = feas_x, feas_y
    # Create colormap
    cm = plt.cm.get_cmap('viridis_r')

    # Function to draw arrows between points in route
    def drawArrow(A, B):
        """
        Draws an arrow from point A to point B
        
        Args:
            A (list): Point A
            B (list): Point B
        Returns:
            None: Directly draws the arrow
        """
        plt.arrow(
            A[0],
            A[1],
            B[0] - A[0],
            B[1] - A[1],
            width=0.00005,
            head_width=0.15,
            head_length=0.08,
            color='black',
            shape='full',
        )

    # Draws scatter plot of feasible points
    sc = plt.scatter(X1, X2, s=80, c=objs, cmap=cm)
    cbar = plt.colorbar(sc)
    cbar.set_label('Objective function', rotation=90)
    title_string = 'D-SDA with k = ' + k
    plt.title(title_string)
    plt.xlabel(ext1_name)
    plt.ylabel(ext2_name)
    plt.show()

    # Plot complete route
    for i in range(len(route) - 1):
        drawArrow(route[i], route[i + 1])

    


def solve_complete_external_enumeration(
    model_function,
    model_args: dict,
    ext_dict: dict,
    ext_logic,
    mip_transformation: bool = False,
    transformation: str = 'bigm',
    feasible_model: str = '',
    points: list = [],
    subproblem_solver: str = 'knitro',
    subproblem_solver_options: dict = {},
    iter_timelimit: float = 10,
    timelimit: float = None,
    gams_output: bool = False,
    tee: bool = False,
    global_tee: bool = True,
    export_csv: bool = False,
):
    """
    Function that computes complete enumeration using the external variable reformulation
    Args:
        model_function (function): function that returns GDP model to be solved
        model_args (dict): Contains the argument values needed for model_function
        ext_dict (dict): Dictionary with Boolean variables to be reformulated (keys) and their corresponding ordered sets (values). Both keys and values are pyomo objects.
        ext_logic: Function that returns a list of lists of the form [a,b], where a is an expressions of the reformulated Boolean variables and b is an equivalent Boolean or indicator variable (b<->a).
        mip_transformation (bool): Whether to solve the enumeration using the external variables applied to the MIP problem insed of the GDP
        transformation (str): Which transformation to apply to the GDP
        feasible_model (str): feasible model name to initialize
        points (list): list of points to carry on enumeration
        subproblem_solver (str): MINLP or NLP solver algorithm
        subproblem_solver_options (dict): MINLP or NLP solver algorithm options
        iter_timelimit (float): time limit in seconds for the solve statement for each iteration
        timelimit (float): time limit in seconds for the algorithm
        gams_output (bool): Determine keeping or not GAMS files
        tee (bool): Display iteration output
        global_tee (bool): Display D-SDA output
        export_csv (bool): Export answer to a csv file
    Returns:
        m2_solved (pyomo.ConcreteModel): Solved Pyomo Model

    """
    results = {}
    feasibles = {}
    t_start = time.perf_counter()
    csv_columns = ['Point', 'x', 'y', 'Objective', 'Status', 'Time', 'Global_Time']
    dict_data = []
    csv_file = (
        'compl_enum_' + str(feasible_model) + '_' + str(subproblem_solver) + '.csv'
    )

    m = model_function(**model_args)
    dict_extvar, num_ext_var, min_allowed, max_allowed = get_external_information(
        m, ext_dict, tee=global_tee
    )

    bounds = []
    for i in range(1, num_ext_var + 1):
        bounds.append(list(range(min_allowed[i], max_allowed[i] + 1)))

    if len(points) == 0:
        points = list(it.product(*bounds))

    if timelimit is None:
        timelimit = 1.5 * iter_timelimit * len(points)

    if global_tee:
        print('\nStarting Complete Enumeration of External Variables')
        print(
            '----------------------------------------------------------------------------------------------'
        )

    for i in points:
        new_result = {}
        m = model_function(**model_args)
        m_init = initialize_model(
            m=m, from_feasible=True, feasible_model=feasible_model, json_path=None
        )
        if mip_transformation:  # If you want a MIP reformulation, go ahead and use it
            csv_file = (
                'compl_enum_'
                + str(feasible_model)
                + '_'
                + str(subproblem_solver)
                + '_'
                + transformation
                + '.csv'
            )
            m_init, dict_extvar = extvars_gdp_to_mip(
                m=m, gdp_dict_extvar=dict_extvar, transformation=transformation
            )
        m_fixed = external_ref(
            m=m_init,
            x=list(i),
            extra_logic_function=ext_logic,
            dict_extvar=dict_extvar,
            mip_ref=mip_transformation,
            tee=False,
        )
        t_remaining = min(iter_timelimit, timelimit - (time.perf_counter() - t_start))
        if t_remaining < 0:  # No time remaining for optimization
            break
        m_solved = solve_subproblem(
            m=m_fixed,
            subproblem_solver=subproblem_solver,
            subproblem_solver_options=subproblem_solver_options,
            timelimit=t_remaining,
            gams_output=gams_output,
            tee=tee,
        )

        results[i] = (m_solved.dsda_status, pe.value(m_solved.obj))
    
        if global_tee:
            print(
                'Evaluated:',
                list(i),
                ' |   Objective:',
                round(pe.value(m_solved.obj), 5),
                ' |   Global Time:',
                round(time.perf_counter() - t_start, 2),
                ' |   Status:',
                m_solved.dsda_status,
            )

        if m_solved.dsda_status == 'Optimal':
            feasibles[i] = float(pe.value(m_solved.obj))

        if export_csv:
            dir_path = os.path.dirname(os.path.abspath(__file__))
            csv_file = os.path.join(dir_path, "../../results", csv_file)
            if m_solved.dsda_status != 'FBBT_Infeasible':
                new_result = {
                    'Point': list(i),
                    'x': i[0],
                    'y': i[1],
                    'Objective': pe.value(m_solved.obj),
                    'Status': m_solved.dsda_status,
                    'Time': m_solved.results.solver.user_time,
                    'Global_Time': time.perf_counter() - t_start,
                }
                dict_data.append(new_result)
            try:
                with open(csv_file, 'w') as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
                    writer.writeheader()
                    for data in dict_data:
                        writer.writerow(data)
            except IOError:
                print("I/O error")

        if time.perf_counter() - t_start > timelimit:
            break

    int_feasibles = {}
    for i in feasibles:
        if not isnan(feasibles[i]):
            int_feasibles[i] = feasibles[i]

    # If there are feasible integer combinations resolve with the best found
    if int_feasibles:
        minimum = min(int_feasibles, key=int_feasibles.get)
        m2 = model_function(**model_args)
        m2_init = initialize_model(
            m=m2, from_feasible=True, feasible_model=feasible_model, json_path=None
        )
        if mip_transformation:  # If you want a MIP reformulation, go ahead and use it
            m2_init, dict_extvar = extvars_gdp_to_mip(
                m=m2_init, gdp_dict_extvar=dict_extvar, transformation=transformation
            )
        m2_fixed = external_ref(
            m=m2_init,
            x=list(minimum),
            extra_logic_function=ext_logic,
            dict_extvar=dict_extvar,
            mip_ref=mip_transformation,
            tee=False,
        )
        m2_solved = solve_subproblem(
            m=m2_fixed,
            subproblem_solver=subproblem_solver,
            subproblem_solver_options=subproblem_solver_options,
            timelimit=iter_timelimit,
            gams_output=gams_output,
            tee=tee,
        )
        if (
            not mip_transformation
        ):  # Error generating json file with MINLP fixed problems
            _ = generate_initialization(m_solved)

        t_end = time.perf_counter() - t_start
        m2_solved.total_time = t_end

        print(m2_solved.results)
        if export_csv:
            final = {
                'Point': list(minimum),
                'x': minimum[0],
                'y': minimum[1],
                'Objective': pe.value(m2_solved.obj),
                'Status': 'Final',
                'Time': m2_solved.results.solver.user_time,
                'Global_Time': t_end,
            }
            dict_data.append(final)
            try:
                with open(csv_file, 'w') as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
                    writer.writeheader()
                    for data in dict_data:
                        writer.writerow(data)
            except IOError:
                print("I/O error")

        if global_tee:
            print(
                '----------------------------------------------------------------------------------------------'
            )
            print('Objective:', round(pe.value(m2_solved.obj), 5))
            print('External variables:', list(minimum))
            print('Execution time [s]:', round(t_end, 2))

        return m2_solved

    return None
