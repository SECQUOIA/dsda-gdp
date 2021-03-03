import pyomo.environ as pe
from pyomo.core.base.misc import display
from pyomo.opt.base.solvers import SolverFactory

# The objective of this minimal working example is to show how dependent variables are currently being declared and initialized

# Model declaration
m = pe.ConcreteModel('MWE_autoinit')

# Parameters
m.x = pe.Param(initialize=2)
m.y = pe.Param(initialize=5)

# Independent variable declaration
m.Beta1 = pe.Var(within=pe.Reals, initialize=1)
m.Beta2 = pe.Var(within=pe.Reals, initialize=1)

# Dependent variable initialization in equation 1
y_hat_init = pe.value(m.Beta1) + pe.value(m.Beta2)*m.x

# Dependent variable declaration
m.y_hat = pe.Var(within=pe.Reals, initialize=y_hat_init)

# Constriant (Equation 1)


@m.Constraint()
def Equation1(m):
    return m.y_hat == m.Beta1 + m.Beta2*m.x

# Objective


@m.Objective()
def obj(m):
    return (m.y - m.y_hat)**2


# Solver
opt = SolverFactory('ipopt')
opt_success = opt.solve(m)

# Is it posible to avoid the explicit initialization declaration done in line 20?
# Could y_hat be initializated directly using Equation1?
