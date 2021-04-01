import copy
import pyomo.environ as pe


m = pe.ConcreteModel()

m.N = pe.RangeSet(1, 3)
m.YF = pe.BooleanVar(m.N)

def _logical_constraint(m):
    return pe.exactly(1,m.YF)

m.logical_constraint = pe.LogicalConstraint(rule=_logical_constraint)

m2 = copy.deepcopy(m)

pe.TransformationFactory('core.logical_to_linear').apply_to(m)

pe.TransformationFactory('core.logical_to_linear').apply_to(m2)
