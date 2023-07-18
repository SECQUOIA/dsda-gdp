# os module provides a portable way of using operating system dependent functionality. 
# It allows you to interface with the underlying operating system that Python is running on.
import os

# Pyomo is a Python-based, open-source optimization modeling language. 
# You use it to define problems from a variety of areas in a format that is easy to read and understand.
# Here, you're importing the whole module as "pe" for convenience.
import pyomo.environ as pe

# The 'display' function from 'pyomo.core.base.misc' is used to print the solution of the model to the console.
from pyomo.core.base.misc import display

# This is a function that can be used to convert boolean variables to binary variables. 
# In other words, it's used to perform a logical to linear transformation.
from pyomo.core.plugins.transform.logical_to_linear import \
    update_boolean_vars_from_binary

# Disjunct and Disjunction are parts of the Generalized Disjunctive Programming (GDP) module of Pyomo. 
# A disjunction is a constraint that ensures that at least one of its disjuncts is satisfied. 
# A disjunct is a logical condition that represents a particular mode of operation in the system.
from pyomo.gdp import Disjunct, Disjunction

# SolverFactory is a factory class that provides a general interface for creating instances of solvers.
# A solver is a piece of software that takes a mathematical model and solves it.
from pyomo.opt.base.solvers import SolverFactory

# This function is used to create a concrete model for the optimization problem.
def build_small_batch():
    # Number of potential parallel units (range from 1 to 3)
    NK = 3 

    # Create a concrete model
    m = pe.ConcreteModel()

    # Initialize sets: products, stages of production, and potential number of parallel units
    m.i = pe.Set(initialize=['a', 'b'])  
    m.j = pe.Set(initialize=['mixer', 'reactor','centrifuge']) 
    m.k = pe.RangeSet(NK)    

    # Initialize parameters: horizon time, lower and upper bounds for batch unit sizes
    m.h = pe.Param(initialize=6000)
    m.vlow = pe.Param(initialize=250)  
    m.vupp = pe.Param(initialize=2500)

    # Initialize the demand for each product
    m.q = pe.Param(m.i, initialize={'a': 200000, 'b': 150000})

    # Initialize the cost coefficients for each production stage
    m.alpha = pe.Param(
        m.j, initialize={'mixer': 250, 'reactor': 500, 'centrifuge': 340})

    # Initialize the cost exponents for each production stage
    m.beta = pe.Param(
        m.j, initialize={'mixer': 0.6, 'reactor': 0.6, 'centrifuge': 0.6})

    # Initialize the number of parallel units using a function to calculate the coefficients
    def coeff_init(m, k):
        return pe.log(k)
    m.coeff = pe.Param(m.k, initialize=coeff_init)

    # Initialize the size factor for each product at each stage
    s_init = {('a', 'mixer'): 2, ('a', 'reactor'): 3, ('a', 'centrifuge'): 4,
              ('b', 'mixer'): 4, ('b', 'reactor'): 6, ('b', 'centrifuge'): 3}
    m.s = pe.Param(m.i, m.j, initialize=s_init)

    # Initialize the processing time for each product at each stage
    t_init = {('a', 'mixer'): 8, ('a', 'reactor'): 20, ('a', 'centrifuge'): 4,
              ('b', 'mixer'): 10, ('b', 'reactor'): 12, ('b', 'centrifuge'): 3}
    m.t = pe.Param(m.i, m.j, initialize=t_init)

    # Initialize variables: stage existence, activation of coefficient, volume of stage,
    # batch size of product, cycle time of product, and number of units in parallel stage
    m.Y = pe.BooleanVar(m.k, m.j)
    m.coeffval = pe.Var(m.k, m.j,  within=pe.NonNegativeReals, bounds=(0, pe.log(NK)))
    m.v = pe.Var(m.j, within=pe.NonNegativeReals, bounds=(pe.log(m.vlow), pe.log(m.vupp)))
    m.b = pe.Var(m.i, within=pe.NonNegativeReals)
    m.tl = pe.Var(m.i, within=pe.NonNegativeReals)
    m.n = pe.Var(m.j, within=pe.NonNegativeReals)

    # Define constraints: volume requirement in each stage, cycle time for each product,
    # production time, number of units, limit on parallel units.

    # Constraints

    # Volume requirement in each stage for each product
    @m.Constraint(m.i, m.j)
    def vol(m, i, j):
        return m.v[j] >= pe.log(m.s[i, j]) + m.b[i]

    # Cycle time for each product in each stage
    @m.Constraint(m.i, m.j)
    def cycle(m, i, j):
        return m.n[j] + m.tl[i] >= pe.log(m.t[i, j])

    # Constraint for total production time (should not exceed horizon time)
    @m.Constraint()
    def time(m):
        return sum(m.q[i]*pe.exp(m.tl[i]-m.b[i]) for i in m.i) <= m.h

    # The number of units in each stage is equal to the sum of coefficients
    @m.Constraint(m.j)
    def units(m, j):
        return m.n[j] == sum(m.coeffval[k, j] for k in m.k)

    # Logical constraint to ensure that only one choice for parallel units is feasible for each stage
    @m.LogicalConstraint(m.j)
    def lim(m, j):
        return pe.exactly(1, m.Y[1, j], m.Y[2, j], m.Y[3, j])

    # Disjunction rules for unit existence and non-existence
    def build_existence_equations(disjunct, k, j):
        m = disjunct.model()
        @disjunct.Constraint()
        def coeffval_act(disjunct):  # Activation of coefficient value if unit exists
            return m.coeffval[k, j] == m.coeff[k]

    def build_not_existence_equations(disjunct, k, j):
        m = disjunct.model()
        @disjunct.Constraint()
        def coeffval_deact(disjunct):  # Deactivation of coefficient value if unit does not exist
            return m.coeffval[k, j] == 0

    # Create disjuncts for unit existence and non-existence
    m.Y_exists = Disjunct(m.k, m.j, rule=build_existence_equations)
    m.Y_not_exists = Disjunct(m.k, m.j, rule=build_not_existence_equations)

    # Create disjunction to represent unit existence or non-existence
    @m.Disjunction(m.k, m.j)
    def Y_exists_or_not(m, k, j):
        return [m.Y_exists[k, j], m.Y_not_exists[k, j]]

    # Associate Boolean variables with with disjunction
    for k in m.k:
        for j in m.j:
            m.Y[k, j].associate_binary_var(m.Y_exists[k, j].indicator_var)

    # Objective function to minimize cost
    def obj_rule(m):
        return sum(m.alpha[j]*(pe.exp(m.n[j] + m.beta[j]*m.v[j])) for j in m.j)

    # Define the objective
    m.obj = pe.Objective(rule=obj_rule, sense=pe.minimize)

    return m


def external_ref(m, x, logic_expr=None):
    # Create a dictionary to hold external reference for each stage
    ext_var = {}
    p = 0
    for j in m.j:
        # Assign each stage j to its corresponding external reference from x
        ext_var[j] = x[p]
        p = p+1

    # Go through each unit k and stage j
    for k in m.k:
        for j in m.j:
            # If unit k is the selected unit for stage j according to ext_var
            if k == ext_var[j]:
                # Fix the boolean variable Y[k,j] to True, meaning unit k is chosen for stage j
                m.Y[k, j].fix(True)
                # Fix the indicator variable of Y_exists[k, j] to True, because unit k exists
                m.Y_exists[k, j].indicator_var.fix(True)  
                # Fix the indicator variable of Y_not_exists[k, j] to False, because unit k does not exist
                m.Y_not_exists[k, j].indicator_var.fix(False)  
            else:
                # If unit k is not the selected unit for stage j according to ext_var
                # Fix the boolean variable Y[k,j] to False, meaning unit k is not chosen for stage j
                m.Y[k, j].fix(False)
                # Fix the indicator variable of Y_exists[k, j] to False, because unit k does not exist
                m.Y_exists[k, j].indicator_var.fix(False)  
                # Fix the indicator variable of Y_not_exists[k, j] to True, because unit k exists
                m.Y_not_exists[k, j].indicator_var.fix(True)  

    # Applying transformations
    # Convert logical constraints to linear constraints
    pe.TransformationFactory('core.logical_to_linear').apply_to(m)
    # Fix the state of all disjuncts to match their current indicator_var values
    pe.TransformationFactory('gdp.fix_disjuncts').apply_to(m)
    # Deactivate constraints that do not affect the feasible region nor the optimal solution
    pe.TransformationFactory('contrib.deactivate_trivial_constraints').apply_to(
        m, tmp=False, ignore_infeasible=True)

    return m



def solve_with_minlp(m, transformation='bigm', minlp='baron', timelimit=10):
    # Convert logical constraints to linear constraints
    pe.TransformationFactory('core.logical_to_linear').apply_to(m)
    # Apply the selected transformation (by default, the Big-M method) to the model
    transformation_string = 'gdp.' + transformation
    pe.TransformationFactory(transformation_string).apply_to(m)

    # Set up the directory for automatically generated GAMS files
    dir_path = os.path.dirname(os.path.abspath(__file__))
    gams_path = os.path.join(dir_path, "gamsfiles/")
    if not(os.path.exists(gams_path)):
        print('Directory for automatically generated files ' +
              gams_path + ' does not exist. We will create it')
        os.makedirs(gams_path)

    # Configure the solver
    solvername = 'gams'
    # Create a solver instance (by default, the BARON solver is used)
    opt = SolverFactory(solvername, solver=minlp)
    # Solve the model
    m.results = opt.solve(m, tee=True,
                          # Uncomment the following lines if you want to save GAMS models
                          # keepfiles=True,
                          # tmpdir=gams_path,
                          # symbolic_solver_labels=True,
                          add_options=[
                              'option reslim = ' + str(timelimit) + ';'
                              'option optcr = 0.0;'
                              # Uncomment the following lines to setup IIS computation of BARON through option file
                              # 'GAMS_MODEL.optfile = 1;'
                              # '\n'
                              # '$onecho > baron.opt \n'
                              # 'CompIIS 1 \n'
                              # '$offecho'
                              # 'display(execError);'
                          ])
    # Update the Boolean variables in the model based on the binary variable values in the solution
    update_boolean_vars_from_binary(m)
    return m



if __name__ == "__main__":
    # Build the initial mathematical model
    m = build_small_batch()
    # Solve the model using the defined 'bigm' transformation method, the 'baron' solver and a time limit of 120 seconds
    # m_solved = solve_with_minlp(m, transformation='bigm', minlp='baron', timelimit=120)

    # External reference test (this test can be deleted)
    # Generate a new model by specifying the number of parallel units for each stage directly (i.e., [1,2,3]), 
    # and applying logical to linear transformation, GDP fix disjuncts and deactivation of trivial constraints
    newmodel = external_ref(m, [1, 2, 3], logic_expr=None)

