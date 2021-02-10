from pyomo.core.base.misc import display
import pyomo.environ as py
from pyomo.opt.base.solvers import SolverFactory

#INPUTS
NT=5 #Size of the superstructure
InitialNumerOfReactors=1 #initialization for yf
InitialLocationOfRecycle=1 #initialization for yr

#PYOMO MODEL
m=py.AbstractModel(name="MINLP_REACTOR")

#SETS
m.I=py.Set(initialize=['A','B']) #Set of components (index i)
m.N=py.RangeSet(1,NT,1) #Set of units in the superstructure (index n,j)
#m.J=py.SetOf(m.N) #Set of units in the superstructure (index n,j)

#PARAMETERS 
m.k=py.Param(initialize=2) #kientic constant [L/(mol*s)]
m.order1=py.Param(initialize=1)#partial order of reacton 1
m.order2=py.Param(initialize=1)#partial order of reaction 2
m.QF0=py.Param(initialize=1) #Inlet volumetric flow [L/s]
C0_Def={'A':0.99,'B':0.01}
m.C0=py.Param(m.I,initialize=C0_Def) #Initial concentration of reagents [mol/L]
def F0_Def(m,I):
    return m.C0[I]*m.QF0
m.F0=py.Param(m.I,initialize=F0_Def) #Inlet molar flow [mol/s]



#BINARY VARIABLES
def yf_Init(m,N): #initialization
    if N==InitialNumerOfReactors:
        return 1
    else:
        return 0
m.yf=py.Var(m.N,within=py.Binary,initialize=yf_Init) #Existence of an unreacted feed in unit n
def yr_Init(m,N):#initialization
    if N==InitialLocationOfRecycle:
        return 1
    else:
        return 0
m.yr=py.Var(m.N,within=py.Binary,initialize=yr_Init) #Existence of recycle flow in unit n
def yp_Init(m,N):#initialization
    return 1-sum(m.yf[J] for J in m.N if J>=1 and J<=N-1)
m.yp=py.Var(m.N,within=py.NonNegativeReals,initialize=yp_Init) #Unit operation in n


#REAL VARIABLES
m.Q=py.Var(m.N,within=py.NonNegativeReals,bounds=(0,10)) #Outlet flow rate of the superstricture unit [L/s]







instance=m.create_instance()
instance.display()




