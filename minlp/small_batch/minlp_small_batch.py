import pyomo.environ as pe
from pyomo.core.base.misc import display
from pyomo.opt.base.solvers import SolverFactory
import os

def build_small_batch(NK=3):

    # Model
    m = pe.ConcreteModel()

    # Sets
    m.i = pe.Set(initialize=['a','b'])  # Set of products
    m.j = pe.Set(initialize=['mixer', 'reactor', 'centrifuge']) # Set of stages
    m.k = pe.RangeSet(NK)    # Set of potential number of parallel units

    # Parameters and Scalars

    m.h = pe.Param(initialize=6000) # Horizon time  (available time hrs)
    m.vlow = pe.Param(initialize=250) # Lower bound for size of batch unit
    m.vupp = pe.Param(initialize=2500)  # Upper bound for size of batch unit

    m.q = pe.Param(m.i, initialize={'a':200000, 'b':150000})    # Demand of product i
    m.alpha = pe.Param(m.j, initialize={'mixer':250, 'reactor':500, 'centrifuge':340})  # Cost coefficient for batch units
    m.beta = pe.Param(m.j, initialize={'mixer':0.6, 'reactor':0.6, 'centrifuge':0.6})   # Cost exponent for batch units

    def coeff_init(m,k):
        return pe.log(k)

    m.coeff = pe.Param(m.k, initialize=coeff_init)  # Represent number of parallel units

    s_init = {('a','mixer'):2, ('a','reactor'):3, ('a','centrifuge'):4, 
              ('b','mixer'):4, ('b','reactor'):6, ('b','centrifuge'):3}

    m.s = pe.Param(m.i, m.j, initialize=s_init)   # Size factor for product i in stage j (kg per l)

    t_init = {('a','mixer'):8, ('a','reactor'):20, ('a','centrifuge'):4, 
              ('b','mixer'):10, ('b','reactor'):12, ('b','centrifuge'):3}

    m.t = pe.Param(m.i, m.j, initialize=t_init)   # Processing time of product i in batch j   (hrs)

    # Variables
    m.y = pe.Var(m.k, m.j, within=pe.Binary)    # Stage existence
    m.v = pe.Var(m.j, within=pe.NonNegativeReals, bounds=(pe.log(m.vlow),pe.log(m.vupp))) # Volume of stage j 
    m.b = pe.Var(m.i, within=pe.NonNegativeReals) # Batch size of product i
    m.tl = pe.Var(m.i, within=pe.NonNegativeReals)  # Cycle time of product i
    m.n = pe.Var(m.j, within=pe.NonNegativeReals)   # Number of units in parallel stage j

    # Constraints

    # Volume requirement in stage j
    @m.Constraint(m.i, m.j)
    def vol(m,i,j): 
        return m.v[j] >= pe.log(m.s[i,j]) + m.b[i]

    # Cycle time for each product i
    @m.Constraint(m.i, m.j)
    def cycle(m,i,j): 
        return m.n[j] + m.tl[i] >= pe.log(m.t[i,j])

    # Constraint for production time
    @m.Constraint()
    def time(m): 
        return sum(m.q[i]*pe.exp(m.tl[i]-m.b[i]) for i in m.i) <= m.h

    # Relating number of units to 0-1 variables
    @m.Constraint(m.j)
    def units(m,j): 
        return m.n[j] == sum(m.coeff[k]*m.y[k,j] for k in m.k)

    # Only one choice for parallel units is feasible
    @m.Constraint(m.j)
    def lim(m,j): 
        return 1 == sum(m.y[k,j] for k in m.k)

    # Objective
    def obj_rule(m):
        return sum(m.alpha[j]*(pe.exp(m.n[j] + m.beta[j]*m.v[j])) for j in m.j)

    m.obj = pe.Objective(rule=obj_rule, sense=pe.minimize)

    return m


def solve_minlp(m, minlp='baron', timelimit=10):

    # Solution step
    dir_path = os.path.dirname(os.path.abspath(__file__))
    gams_path = os.path.join(dir_path, "gamsfiles/")
    if not(os.path.exists(gams_path)):
        print('Directory for automatically generated files ' +
              gams_path + ' does not exist. We will create it')
        os.makedirs(gams_path)

    solvername = 'gams'
    opt = SolverFactory(solvername, solver=minlp)
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

    return m

if __name__ == "__main__":
    NK = 3
    m = build_small_batch(NK=NK)
    m_solved = solve_minlp(m, minlp='baron', timelimit=120)
    
