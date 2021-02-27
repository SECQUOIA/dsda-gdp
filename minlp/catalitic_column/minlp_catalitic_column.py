import pyomo.environ as pe
import networkx as nx
import matplotlib.pyplot as plt
from math import (sqrt, pi)
from pyomo.core.base.misc import display
from pyomo.opt.base.solvers import SolverFactory
import os

def minlp_catalitic_column(NT=22,  visualize=False):

    FT = 2
    KT = 3

    # PYOMO MODEL
    m = pe.ConcreteModel(name='minlp_catalitic_column')

    # ______________________________ Sections 1-4 ______________________________
    # Main variables and sets for definig the system
    # Hydraulic equations calculation

    # Sets
    m.I = pe.Set(initialize=['iButene', 'Ethanol', 'nButene','ETBE'])  # Set of components
    m.F = pe.RangeSet(1, FT)  # Set of feeds
    m.N = pe.RangeSet(1, NT)  # Set of all stages in the column

    # Variables
    m.L = pe.Var(m.N,  within=pe.NonNegativeReals)  # Flow of liquid [mol/min]
    m.V = pe.Var(m.N,  within=pe.NonNegativeReals)  # Flow of vapor [mol/min]
    m.x = pe.Var(m.I, m.N,  within=pe.NonNegativeReals)  # Molar composition of liquid [*]
    m.y = pe.Var(m.I, m.N,  within=pe.NonNegativeReals)  # Molar composition of vapor [*]
    m.Temp = pe.Var(m.N,  within=pe.NonNegativeReals)   # Operation temperature [K]
    m.P = pe.Var(m.N,  within=pe.NonNegativeReals)  # Stage pressure [bar]
    m.Z = pe.Var(m.N,  within=pe.NonNegativeReals)  # Compressibility coefficient [*]
    m.RR = pe.Var(within=pe.NonNegativeReals)   # Reflux ratio [*]
    m.Qc = pe.Var(within=pe.NonNegativeReals)   # Condensator duty [kJ/min]
    m.Qr = pe.Var(within=pe.NonNegativeReals)   # Reboiler duty [kJ/min]
    m.BR = pe.Var(within=pe.NonNegativeReals)   # Boil up [*]

    # Hydraulic parameters
    m.da = pe.Param(initialize=0.002)   # Plate hole diameter [m]
    m.ep = pe.Param(initialize=0.002)   # Plate thickness [m]
    m.pitch = pe.Param(initialize=0.009)   # Distance between plate holes [m]
    m.Sfactor = pe.Param(initialize=0.15)   # Safety height factor
    m.poro = pe.Param(initialize=0.907*sqrt(m.da/m.pitch)) # Plate pororsity [*]
    m.K0 = pe.Param(initialize=(880.6-(67.7*m.da/m.ep)+(7.32*((m.da/m.ep)**2))-(0.338*((m.da/m.ep)**3)))*0.001) # Hole coefficient [*]

    # Hydraulic variables
    m.D = pe.Var(within=pe.NonNegativeReals)   # Column diameter [m]
    m.hw = pe.Var(within=pe.NonNegativeReals)   # Weir height [m]
    m.HS = pe.Var(within=pe.NonNegativeReals)   # Plate height [m]
    m.Htotal = pe.Var(within=pe.NonNegativeReals)   # Total column height [m]
    m.At = pe.Var(within=pe.NonNegativeReals)   # Active area [m**2]
    m.Ad = pe.Var(within=pe.NonNegativeReals)   # Weir area [m**2]
    m.Lw = pe.Var(within=pe.NonNegativeReals)   # Weir lenght [m]
    m.A0 = pe.Var(within=pe.NonNegativeReals)   # Holed area [m**2]

    # Hydraulic constraints
    @m.Constraint()
    def EqHwmin(m):
        return m.hw >= 0.05*m.HS

    @m.Constraint()
    def EqHwmax(m):
        return m.hw <= m.HS/3

    @m.Constraint()
    def EqAt(m):
        return m.At == pe.sqrt(m.D/2)*(pi-(1.854590-0.96))
    
    @m.Constraint()
    def EqAd(m):
        return m.Ad == pe.sqrt(m.D/2)*(0.5*(1.854590-0.96))

    @m.Constraint()
    def EqLw(m):
        return m.Lw == 0.8*m.D
    
    @m.Constraint()
    def EqA0(m):
        return m.A0 == m.At*m.poro
    
    # Butene feed
    m.FB = pe.Param(initialize=5.774)   # Butene flow in feed [mol/min]
    zb_ibutane = 30
    zb_init = {'iButene':zb_ibutane, 'Ethanol':0, 'nButene':100-zb_ibutane, 'ETBE':0}
    m.zb = pe.Param(m.I, initialize=zb_init)

    # Ethanol feed
    m.FE = pe.Param(initialize=1.7118)   # Ethanol flow in feed [mol-h]
    ze_init = {'iButene':0, 'Ethanol':100, 'nButene':0, 'ETBE':0}
    m.ze = pe.Param(m.I, initialize=ze_init)

    # Operation parameters
    m.Pop = pe.Param(initialize=9.5)   # Condenser pressure [bar]
    m.TaliB = pe.Param(initialize=323)   # Butene feed temperature [K]
    m.TaliE = pe.Param(initialize=342.38)   # Ethanol feed temperature [K]
    m.xBetbe = pe.Param(initialize=83)   # Desired composition of ETBE in bottom [*]
    m.MCR = pe.Param(initialize=1)   # Constant flow keep in reboiler and condenser [mol]
    m.cR = pe.Param(initialize=0.00008314)   # Ideal gas constant [m**3*bar/K*mol]
    m.hour = pe.Param(initialize=60)   # Seconds in an hour [60]

    # Composition restriction in bottoms
    @m.Constraint(m.N)
    def pureza0(m,n):
        if n == NT:
            return m.x['ETBE',n] >= m.xBetbe
        else:
            return pe.Constraint.Skip

    # ______________________________ Section 5 ______________________________
    # Saturation pressures usinn Antoine equation

    # Constants for expanded Antoine equation
    C1a_init = {'iButene':66.4970745, 'Ethanol':61.7910745, 'nButene':40.3230745, 'ETBE':52.67507454}
    m.C1a = pe.Param(m.I, initialize=C1a_init)
    C2a_init = {'iButene':-4634.1, 'Ethanol':-7122.3, 'nButene':-4019.2, 'ETBE':-5820.2}
    m.C2a = pe.Param(m.I, initialize=C2a_init)
    C3a_init = {'iButene':0, 'Ethanol':0, 'nButene':0, 'ETBE':0}
    m.C3a = pe.Param(m.I, initialize=C3a_init)
    C4a_init = {'iButene':0, 'Ethanol':0, 'nButene':0, 'ETBE':0}
    m.C4a = pe.Param(m.I, initialize=C4a_init)
    C5a_init = {'iButene':-8.9575, 'Ethanol':-7.1424, 'nButene':-4.5229, 'ETBE':-6.1343}
    m.C5a = pe.Param(m.I, initialize=C5a_init)
    C6a_init = {'iButene':1.3413*10**-5, 'Ethanol':2.8853*10**-6, 'nButene':4.8833*10**-17, 'ETBE':2.1405*10**-17}
    m.C6a = pe.Param(m.I, initialize=C6a_init)
    C7a_init = {'iButene':2, 'Ethanol':2, 'nButene':6, 'ETBE':6}
    m.C7a = pe.Param(m.I, initialize=C7a_init)

    # Antoine equation
    m.Psat = pe.Var(m.I, m.N, within=pe.NonNegativeReals)   # Saturation pressure [bar]
    @m.Constraint(m.I, m.N)
    def EqPsat(m,i,n):
        return m.Psat[i,n] == pe.exp(m.C1a[i] + (m.C2a[i]/(m.Temp[n] + m.C3a[i])) + (m.C4a[i]*m.Temp[n]) + ((m.C5a[i]*pe.log(m.Temp[n])) + (m.C6a[i]*(m.Temp[n]**m.C7a[i]))))

    # ______________________________ Section 6 ______________________________
    # Calculation of liquid density using IK-CAPI equation
    # Calculation of liquid density using critic DIPPR equation
    # Calculation of gas density using corrected ideal gas equation
    
    # Constants for DIPPR equation
    MW_init = {'iButene':56.10752, 'Ethanol':46.06904, 'nButene':56.10752, 'ETBE':102.17656}
    m.MW = pe.Param(m.I, initialize=MW_init)    # Molecular weight [kg/kmol]
    Tcrit_init = {'iButene':417.9, 'Ethanol':516.2, 'nButene':419.6, 'ETBE':509.4}
    m.Tcrit = pe.Param(m.I, initialize=Tcrit_init) # Critic temperature [K]
    Pcrit_init = {'iButene':38.98675, 'Ethanol':60.35675, 'nButene':39.18675, 'ETBE':28.32675}
    m.Pcrit = pe.Param(m.I, initialize=Pcrit_init) # Critic pressure [bar]
    C1rh_init = {'iButene':8.9711123119, 'Ethanol':-2.932961888*10**-2, 'nButene':5.956235579, 'ETBE':-1.323678817*10**-1}
    m.C1rh = pe.Param(m.I, initialize=C1rh_init)
    C2rh_init = {'iButene':0, 'Ethanol':6.9361857406*10**-4, 'nButene':0, 'ETBE':2.1486345729*10**-3}
    m.C2rh = pe.Param(m.I, initialize=C2rh_init)
    C3rh_init = {'iButene':0, 'Ethanol':-1.962897037*10**-6, 'nButene':0, 'ETBE':-6.092181735*10**-6}
    m.C3rh = pe.Param(m.I, initialize=C3rh_init)
    C4rh_init = {'iButene':0, 'Ethanol':2.089632106*10**-9, 'nButene':0, 'ETBE':6.4627035532*10**-9}
    m.C4rh = pe.Param(m.I, initialize=C4rh_init)
    C5rh_init = {'iButene':0, 'Ethanol':0, 'nButene':0, 'ETBE':0}
    m.C5rh = pe.Param(m.I, initialize=C5rh_init)
    C6rh_init = {'iButene':-1.4666609*10**-10, 'Ethanol':0, 'nButene':-9.3717935*10**-11, 'ETBE':0}
    m.C6rh = pe.Param(m.I, initialize=C6rh_init)
    C7rh_init = {'iButene':1.286186216*10**-12, 'Ethanol':0, 'nButene':8.150339357*10**-13, 'ETBE':0}
    m.C7rh = pe.Param(m.I, initialize=C7rh_init)
    C8rh_init = {'iButene':-4.33826109*10**-15, 'Ethanol':0, 'nButene':-2.72421122*10**-15, 'ETBE':0}
    m.C8rh = pe.Param(m.I, initialize=C8rh_init)
    C9rh_init = {'iButene':6.619652613*10**-18, 'Ethanol':0, 'nButene':4.115761136*10**-18, 'ETBE':0}
    m.C9rh = pe.Param(m.I, initialize=C9rh_init)
    C10rh_init = {'iButene':-3.8362103001*10**-21, 'Ethanol':0, 'nButene':-2.3593237507*10**-21, 'ETBE':0}
    m.C10rh = pe.Param(m.I, initialize=C10rh_init)
    C1r_init = {'iButene':1.1446, 'Ethanol':1.6288, 'nButene':1.0877, 'ETBE':0.66333}
    m.C1r = pe.Param(m.I, initialize=C1r_init)
    C2r_init = {'iButene':0.2724, 'Ethanol':0.27469, 'nButene':2.6454*10**-1, 'ETBE':2.6135*10**-1}
    m.C2r = pe.Param(m.I, initialize=C2r_init)
    C3r_init = {'iButene':0.28172, 'Ethanol':0.23178, 'nButene':0.2843, 'ETBE':0.28571}
    m.C3r = pe.Param(m.I, initialize=C3r_init)
    C4r_init = {'iButene':0, 'Ethanol':0, 'nButene':0, 'ETBE':0}
    m.C4r = pe.Param(m.I, initialize=C4r_init)

    m.Tcritm = pe.Var(m.N, within=pe.NonNegativeReals)
    @m.Constraint(m.N)
    def EqTcritm(m,n):
        return m.Tcritm[n] == (pe.sqrt(sum((m.x[i,n]/100)*m.Tcrit[i]/(m.Pcrit[i]**0.5) for i in m.I))) / (sum((m.x[i,n]/100)*m.Tcrit[i]/m.Pcrit[i] for i in m.I))
    
    m.rho = pe.Var(m.I, m.N, within=pe.NonNegativeReals) # Liquid molar density [mol/m**3]
    @m.Constraint(m.I, m.N)
    def Eqrho(m,i,n):
        return m.rho[i,n] == (m.C1r[i]/(m.C2r[i]**(1+((1-(m.Temp[n]/m.Tcritm[n]))**m.C4r[i]))))*1000

    m.rhoV = pe.Var(m.N, within=pe.NonNegativeReals) # Vapor molar density [mol/m**3]
    @m.Constraint(m.N)
    def EqrhoV(m,n):
        return m.rhoV[n] == m.P[n]/(0.00008314*m.Temp[n]*(m.Z[n]))

    # ______________________________ Section 7 ______________________________
    # Calculation of superficial tension using critic DIPPR equation

    # Constants for DIPPR equation
    C1sig_init = {'iButene':0.05544, 'Ethanol':0.03764, 'nButene':0.055945, 'ETBE':0.071885}
    m.C1sig = pe.Param(m.I, initialize=C1sig_init)
    C2sig_init = {'iButene':1.2453, 'Ethanol':-2.157*10**-5, 'nButene':1.2402, 'ETBE':2.1204}
    m.C2sig = pe.Param(m.I, initialize=C2sig_init)
    C3sig_init = {'iButene':0, 'Ethanol':1.025*10**-7, 'nButene':0, 'ETBE':-1.5583}
    m.C3sig = pe.Param(m.I, initialize=C3sig_init)
    C4sig_init = {'iButene':0, 'Ethanol':0, 'nButene':0, 'ETBE':0.76657}
    m.C4sig = pe.Param(m.I, initialize=C4sig_init)

    m.sigma = pe.Var(m.N, within=pe.NonNegativeReals) # Liquid-vapor superficial tension [N/m]
    @m.Constraint(m.N)
    def Eqsigma(m,n):
        return m.sigma[n] == sum((m.x[i,n]/100)*m.C1sig[i]*(1-(m.Temp[n]/m.Tcritm[n]))**(m.C2sig[i]+m.C3sig[i]*(m.Temp[n]/m.Tcritm[n])+m.C4sig[i]*((m.Temp[n]/m.Tcritm[n]))**2) for i in m.I)

    # ______________________________ Section 8 ______________________________
    # Calculation of activity coefficient using NRTL model

    a_nrtl_init = {(i,i2):0 for i in m.I for i2 in m.I}
    m.a_nrtl = pe.Param(m.I, m.I, initialize=a_nrtl_init)
    b_nrtl_init = {('iButene','iButene'):0, ('iButene','Ethanol'):623.5810010, ('iButene','nButene'):107.526499, ('iButene','ETBE'):219.73407,
                   ('Ethanol','iButene'):141.9632130, ('Ethanol','Ethanol'):0, ('Ethanol','nButene'):164.57256, ('Ethanol','ETBE'):187.104064,
                   ('nButene','iButene'):-93.24546420, ('nButene','Ethanol'):595.5299820, ('nButene','nButene'):0, ('nButene','ETBE'):226.373398,
                   ('ETBE','iButene'):-172.59152, ('ETBE','Ethanol'):344.481315, ('ETBE','nButene'):-177.88565, ('ETBE','ETBE'):0}
    m.b_nrtl = pe.Param(m.I, m.I, initialize=b_nrtl_init)
    c_nrtl_init = {(i,i2):0.3 for i in m.I for i2 in m.I}
    for i in m.I:
        c_nrtl_init[i,i] = 0
    m.c_nrtl = pe.Param(m.I, m.I, initialize=c_nrtl_init)
    
    def alfa_nrtl_init(m, i, i2):
        if i != i2:
            return m.c_nrtl[i,i2]
        else:
            return pe.Param.Skip

    m.alfa_nrtl = pe.Param(m.I, m.I, initialize=alfa_nrtl_init, within=pe.Any)

    m.tao_nrtl = pe.Var(m.I, m.I, m.N, within=pe.Reals)
    @m.Constraint(m.I, m.I, m.N)
    def Eq_tao_nrtl(m,i,i2,n):
        return m.tao_nrtl[i,i2,n] == m.a_nrtl[i,i2] + (m.b_nrtl[i,i2]/m.Temp[n])

    m.g_nrtl = pe.Var(m.I, m.I, m.N, within=pe.Reals)
    @m.Constraint(m.I, m.I, m.N)
    def Eq_g_nrtl(m,i,i2,n):
        if i != i2:
            return m.g_nrtl[i,i2,n] == pe.exp(-m.alfa_nrtl[i,i2]*m.tao_nrtl[i,i2,n])
        else:
            return pe.Constraint.Skip

    m.gamma = pe.Var(m.I, m.N, within=pe.Reals)
    @m.Constraint(m.I, m.N)
    def Eqgamma(m,comp,n):
        return m.gamma[comp,n] == pe.exp(sum(m.x[comp1,n]*m.tao_nrtl[comp1,comp,n]*
        m.g_nrtl[comp1,comp,n] for comp1 in m.I)/sum(m.x[comp1,n]*
        m.g_nrtl[comp1,comp,n] for comp1 in m.I)+sum(m.x[comp1,n]*
        m.g_nrtl[comp,comp1,n]/sum(m.x[comp2,n]*
        m.g_nrtl[comp2,comp1,n] for comp2 in m.I)*(m.tao_nrtl[comp,comp1,n]-
        sum(m.x[comp2,n]*m.tao_nrtl[comp2,comp1,n]*
        m.g_nrtl[comp2,comp1,n] for comp2 in m.I)/sum(m.x[comp3,n]*
        m.g_nrtl[comp3,comp1,n] for comp3 in m.I)) for comp1 in m.I))

    # ______________________________ Section 9 ______________________________
    # Chemical reaction

    Nu_init = {'iButene':-1, 'Ethanol':-1, 'nButene':0, 'ETBE':1}
    m.Nu = pe.Param(m.I, initialize=Nu_init)    # Stoichiometry coeffients [*]
    m.mcat = pe.Param(initialize=0.4)   # Catalizer mass [kg]
    m.Ketbe = pe.Var(m.N, within=pe.Reals) # Equilibrium constant [*]
    @m.Constraint(m.N)
    def EqKetbe(m,n):
        if n != NT and n != 1:
            return m.Ketbe[n] == pe.exp(10.387+4060.59/(m.Temp[n])
            -2.89055*pe.log(m.Temp[n])-0.01915144*m.Temp[n]
            +0.0000528586*(m.Temp[n]**2)-0.0000000532977*(m.Temp[n]**3))
        else:
            return pe.Constraint.Skip
    
    m.Krate = pe.Var(m.N, within=pe.NonNegativeReals)   # Reaction advance rate [mol/kg_cat*min]
    @m.Constraint(m.N)
    def EqKrate(m,n):
        if n != NT and n != 1:
            return m.Krate[n] == 7.41816*10**15*pe.exp(-60400/(8.314*m.Temp[n]))*m.hour/3600
        else:
            return pe.Constraint.Skip

    m.Ka = pe.Var(m.N, within=pe.NonNegativeReals)  # Adsorption rate
    @m.Constraint(m.N)
    def EqKa(m,n):
        if n != NT and n != 1:
            return m.Ka[n]==pe.exp(-1.0707+1323.1/m.Temp[n])
        else:
            return pe.Constraint.Skip

    m.Rx = pe.Var(m.N, within=pe.Reals)  # Reaction rate [mol/kg_cat*min]
    @m.Constraint(m.N)
    def EqRx(m,n):
        if n != NT and n != 1:
            return m.Rx[n]*((1+m.Ka[n]*m.gamma['Ethanol',n]*m.x['Ethanol',n]/100)**3)*m.Ketbe[n] == (m.Krate[n]*(m.gamma['Ethanol',n]*m.x['Ethanol',n]/100))*((m.Ketbe[n]*m.gamma['iButene',n]*m.x['iButene',n]/100*m.gamma['Ethanol',n]*m.x['Ethanol',n]/100)-(m.gamma['ETBE',n]*m.x['ETBE',n]/100))
        else:
            return pe.Constraint.Skip

    # ______________________________ Section 10 ______________________________
    # Phi calculation

    Omega_init = {'iButene':0.19484, 'Ethanol':0.643558, 'nButene':0.184495, 'ETBE':0.316231}
    m.Omega = pe.Param(m.I, initialize=Omega_init)    # Acentric factor [*]
    TcritSRK_init = {'iButene':417.9, 'Ethanol':514, 'nButene':419.5, 'ETBE':509.4}
    m.TcritSRK = pe.Param(m.I, initialize=TcritSRK_init)    # Critic temperature for Soave-Redlich-Kwong ecuation [K]

    def mEOS_init(m,i):
        return 0.48508+1.55171*m.Omega[i]-0.15613*pe.sqrt(m.Omega[i])
    m.mEOS = pe.Param(m.I, initialize=mEOS_init)

    def biEOS_init(m,i):
        return 0.08664*0.00008314*m.TcritSRK[i]/m.Pcrit[i]
    m.biEOS = pe.Param(m.I, initialize=biEOS_init)

    m.alphaEOS = pe.Var(m.I, m.N, within=pe.NonNegativeReals)
    @m.Constraint(m.I, m.N)
    def EqAlphaEOS(m,i,n):
        return m.alphaEOS[i,n] == pe.sqrt(1+m.mEOS[i]*(1-(m.Temp[n]/m.Tcritm[n])**(1/2)))

    m.aiEOS = pe.Var(m.I, m.N, within=pe.NonNegativeReals)
    @m.Constraint(m.I, m.N)
    def EqaiEOS(m,i,n):
        return m.aiEOS[i,n] == m.alphaEOS[i,n]*0.42747*(pe.sqrt(0.00008314*m.TcritSRK[i]))/m.Pcrit[i]

    m.bEOS = pe.Var(m.N, within=pe.NonNegativeReals)
    @m.Constraint(m.N)
    def EqbEOS(m,n):
        return m.bEOS[n] == sum((m.y[i,n]/100)*m.biEOS[i] for i in m.I)

    m.aEOS = pe.Var(m.N, within=pe.NonNegativeReals)
    @m.Constraint(m.N)
    def EqaEOS(m,n):
        return m.aEOS[n] == sum(sum((m.y[i,n]/100)*(m.y[i2,n]/100)*(m.aiEOS[i,n]*m.aiEOS[i2,n])**0.5 for i2 in m.I) for i in m.I)

    @m.Constraint(m.N)
    def EqVaporZ(m,n):
        return (m.Z[n])**3-(m.Z[n])**2+(m.Z[n])*((m.aEOS[n]*m.P[n]/((0.00008314*m.Temp[n])**2))-(m.bEOS[n]*m.P[n]/(0.00008314*m.Temp[n]))-(m.bEOS[n]*m.P[n]/(0.00008314*m.Temp[n]))**2)-((m.aEOS[n]*m.P[n]/((0.00008314*m.Temp[n])**2)))*(m.bEOS[n]*m.P[n]/(0.00008314*m.Temp[n])) == 0
    
    m.phi = pe.Var(m.I, m.N, within=pe.NonNegativeReals)
    @m.Constraint(m.I, m.N)
    def EqPhi(m,i,n):
        return m.phi[i,n] == pe.exp(((m.Z[n])-1)*m.biEOS[i]/m.bEOS[n]-pe.log((m.Z[n])-m.bEOS[n])-(m.aEOS[n]/m.bEOS[n])*(2*((m.aiEOS[i,n]/m.aEOS[n])**(1/2))-m.biEOS[i]/m.bEOS[n])*pe.log(((m.Z[n])-m.bEOS[n])/(m.Z[n])))

    # ______________________________ Section 11 ______________________________
    # Entalphy calculation

    # Cp constants [kJ/mol*K]
    C1c_init = {'iButene':0.016052191, 'Ethanol':0.00901418, 'nButene':-0.00299356, 'ETBE':-0.014651654}
    m.C1c = pe.Param(m.I, initialize=C1c_init)
    C2c_init = {'iButene':0.000280432, 'Ethanol':0.000214071, 'nButene':0.000353198, 'ETBE':0.000698631}
    m.C2c = pe.Param(m.I, initialize=C2c_init)
    C3c_init = {'iButene':-0.00000010914988, 'Ethanol':-0.000000083903472, 'nButene':-0.00000019904047, 'ETBE':-0.00000044791741}
    m.C3c = pe.Param(m.I, initialize=C3c_init)
    C4c_init = {'iButene':0.0000000000090979164, 'Ethanol':0.0000000000013732704, 'nButene':0.000000000044631288, 'ETBE':0.00000000011636811}
    m.C4c = pe.Param(m.I, initialize=C4c_init)
    C5c_init = {'iButene':0, 'Ethanol':0, 'nButene':0, 'ETBE':0}
    m.C5c = pe.Param(m.I, initialize=C5c_init)
    C6c_init = {'iButene':0, 'Ethanol':0, 'nButene':0, 'ETBE':0}
    m.C6c = pe.Param(m.I, initialize=C6c_init)
    m.Tref = pe.Param(initialize=298.15)    # Reference temperature [K]
    Hform_init = {'iButene':-16.9147, 'Ethanol':-234.963, 'nButene':-0.125604, 'ETBE':-313.9}
    m.Hform = pe.Param(m.I, initialize=Hform_init)  # Formation enthalphy [kJ/mol]
    Tb_init = {'iButene':341.7, 'Ethanol':421.9, 'nButene':342.6, 'ETBE':438.8}
    m.Tb = pe.Param(m.I, initialize=Tb_init)  # Component boiling temperature @ P=9.5bar [K]
    # Enthalphy calculation Integal(CpdT)
    m.HVi = pe.Var(m.I, m.N, within=pe.Reals)
    m.HV = pe.Var(m.N, within=pe.Reals)
    @m.Constraint(m.I, m.N)
    def EqHVi(m,i,n):
        return m.HVi[i,n] == ((m.C1c[i]*(m.Temp[n]-m.Tref)) + ((m.C2c[i]/2)*((m.Temp[n]**2)-(m.Tref**2)))
                            + ((m.C3c[i]/3)*((m.Temp[n]**3)-(m.Tref**3))) + ((m.C4c[i]/4)*((m.Temp[n]**4)-(m.Tref**4)))
                            + ((m.C5c[i]/5)*((m.Temp[n]**5)-(m.Tref**5))) + ((m.C6c[i]/6)*((m.Temp[n]**6)-(m.Tref**6))) + m.Hform[i]
                            + (8.314/1000)*m.Temp[n]*(m.Z[n]-1)+(1+m.mEOS[i])*((m.aEOS[n]**0.5)/m.bEOS[n])*pe.log(m.Z[n]/(m.Z[n]+(m.bEOS[n]*m.P[n]/(0.00008314*m.Temp[n])))))
    
    @m.Constraint(m.N)
    def EqHV(m,n):
        return m.HV[n]==sum(m.HVi[i,n]*m.y[i,n]/100 for i in m.I)

    # Vaporization entalphy constants [kJ/mol]
    C1v_init = {'iButene':32.614, 'Ethanol':55.789, 'nButene':33.774, 'ETBE':45.29}
    m.C1v = pe.Param(m.I, initialize=C1v_init)
    C2v_init = {'iButene':0.38073, 'Ethanol':0.31245, 'nButene':0.5107, 'ETBE':0.27343}
    m.C2v = pe.Param(m.I, initialize=C2v_init)
    C3v_init = {'iButene':0, 'Ethanol':0, 'nButene':-0.17304, 'ETBE':0.21645}
    m.C3v = pe.Param(m.I, initialize=C3v_init)
    C4v_init = {'iButene':0, 'Ethanol':0, 'nButene':0.05181, 'ETBE':-0.11756}
    m.C4v = pe.Param(m.I, initialize=C4v_init)
    C5v_init = {'iButene':0, 'Ethanol':0, 'nButene':0, 'ETBE':0}
    m.C5v = pe.Param(m.I, initialize=C5v_init)

    def Tred_init(m,i):
        return m.Tb[i]/m.Tcrit[i]
    m.Tred = pe.Param(m.I, initialize=Tred_init)    # Reduced temperature [*]

    def alphaEOSb_init(m,i):
        return (1+m.mEOS[i]*(1-(m.Tb[i]/m.Tcrit[i])**(1/2)))**2
    m.alphaEOSb = pe.Param(m.I, initialize=alphaEOSb_init)

    def aiEOSb_init(m,i):
        return m.alphaEOSb[i]*0.42747*((0.00008314*m.TcritSRK[i])**2)/m.Pcrit[i]

    m.aiEOSb = pe.Param(m.I, initialize=aiEOSb_init)

    m.Zboil = pe.Var(m.I, m.N, within=pe.NonNegativeReals)
    @m.Constraint(m.I, m.N)
    def VaporZb(m,i,n):
        return (m.Zboil[i,n])**3-(m.Zboil[i,n])**2+(m.Zboil[i,n])*((m.aiEOSb[i]*m.P[n]/((0.00008314*m.Tb[i])**2))-(m.biEOS[i]*m.P[n]/(0.00008314*m.Tb[i]))-(m.biEOS[i]*m.P[n]/(0.00008314*m.Tb[i]))**2)-((m.aiEOSb[i]*m.P[n]/((0.00008314*m.Tb[i])**2)))*(m.biEOS[i]*m.P[n]/(0.00008314*m.Tb[i])) == 0

    # Vaporization enthalpy
    def DHVap_init(m,i):
        return ( m.C1v[i]*( (1-m.Tred[i])**( m.C2v[i] + (m.C3v[i]*m.Tred[i]) + (m.C4v[i]*(m.Tred[i]**2)) + (m.C5v[i]*(m.Tred[i]**3)) ) ) )

    m.DHVap = pe.Param(m.I, initialize=DHVap_init)    

    def HVib_init(m,i):
        return ( (m.C1c[i]*(m.Tb[i]-m.Tref)) + ((m.C2c[i]/2)*((m.Tb[i]**2)-(m.Tref**2))) + ((m.C3c[i]/3)*((m.Tb[i]**3)-(m.Tref**3))) + ((m.C4c[i]/4)*((m.Tb[i]**4)-(m.Tref**4))) + ((m.C5c[i]/5)*((m.Tb[i]**5)-(m.Tref**5))) + ((m.C6c[i]/6)*((m.Tb[i]**6)-(m.Tref**6))) + m.Hform[i])

    m.HVib = pe.Param(m.I, initialize=HVib_init)
    m.depHvib = pe.Var(m.I, m.N, within=pe.Reals)
    @m.Constraint(m.I, m.N)
    def EqdepHvib(m,i,n):
        return m.depHvib[i,n] == (8.314/1000)*m.Tb[i]*(m.Zboil[i,n]-1)+(1+m.mEOS[i])*((m.aiEOSb[i]**0.5)/m.biEOS[i])*pe.log(m.Zboil[i,n]/(m.Zboil[i,n]+(m.biEOS[i]*m.P[n]/(0.00008314*m.Tb[i]))))

    # Cp contants of liquid [kJ/mol**K]
    C1l_init = {'iButene':0.08768, 'Ethanol':0.10264, 'nButene':0.18205, 'ETBE':0.11096}
    m.C1l = pe.Param(m.I, initialize=C1l_init)
    C2l_init = {'iButene':0.0002171, 'Ethanol':-0.00013963, 'nButene':-0.001611, 'ETBE':0.00031422}
    m.C2l = pe.Param(m.I, initialize=C2l_init)
    C3l_init = {'iButene':-9.15300*10**-7, 'Ethanol':-3.03410*10**-8, 'nButene':1.19630*10**-5, 'ETBE':1.74800*10**-7}
    m.C3l = pe.Param(m.I, initialize=C3l_init)
    C4l_init = {'iButene':2.2660*10**-9, 'Ethanol':2.0386*10**-9, 'nButene':-3.7454*10**-8, 'ETBE':0}
    m.C4l = pe.Param(m.I, initialize=C4l_init)
    C5l_init = {'iButene':0, 'Ethanol':0, 'nButene':4.5027*10**-11, 'ETBE':0}
    m.C5l = pe.Param(m.I, initialize=C5l_init)
    m.HLi = pe.Var(m.I, m.N, within=pe.Reals)
    m.HL = pe.Var(m.N, within=pe.Reals)
    @m.Constraint(m.I, m.N)
    def EqHLi(m,i,n):
        return m.HLi[i,n] == m.HVib[i]-m.DHVap[i]+((m.C1l[i]*(m.Temp[n]-m.Tb[i])) + ((m.C2l[i]/2)*((m.Temp[n]**2)-(m.Tb[i]**2)))+((m.C3l[i]/3)*((m.Temp[n]**3)-(m.Tb[i]**3))) + ((m.C4l[i]/4)*((m.Temp[n]**4)-(m.Tb[i]**4)))+((m.C5l[i]/5)*((m.Temp[n]**5)-(m.Tb[i]**5))))+m.depHvib[i,n]

    @m.Constraint(m.N)
    def EqHL(m,n):
        return m.HL[n] == sum(m.HLi[i,n]*m.x[i,n]/100 for i in m.I)

    # ______________________________ Section 12 ______________________________
    # Entalphy in feed calculation

    def HV_b_init(m,i):
        return ( (m.C1c[i]*(m.TaliB-m.Tref)) + ((m.C2c[i]/2)*((m.TaliB**2)-(m.Tref**2))) + ((m.C3c[i]/3)*((m.TaliB**3)-(m.Tref**3))) + ((m.C4c[i]/4)*((m.TaliB**4)-(m.Tref**4))) + ((m.C5c[i]/5)*((m.TaliB**5)-(m.Tref**5))) + ((m.C6c[i]/6)*((m.TaliB**6)-(m.Tref**6))) + m.Hform[i])

    m.HV_b = pe.Param(m.I, initialize=HV_b_init)    # Vapor entalphy of feed [kJ/mol]

    def Tred_b_init(m,i):
        return m.TaliB/m.Tcrit[i]

    m.Tred_b = pe.Param(m.I, initialize=Tred_b_init)    # Reduced temperature of feed [*]

    def DHVap_b_init(m,i):
        return ( m.C1v[i]*( (1-m.Tred_b[i])**( m.C2v[i] + (m.C3v[i]*m.Tred_b[i]) + (m.C4v[i]*(m.Tred_b[i]**2)) + (m.C5v[i]*(m.Tred_b[i]**3)) ) ) )

    m.DHVap_b = pe.Param(m.I, initialize=DHVap_b_init)    # Vaporization entalphy of feed [kJ/mol]

    def HL_b_init(m,i):
        return m.HV_b[i]-m.DHVap_b[i]

    m.HL_b = pe.Param(m.I, initialize=HL_b_init)    # Liquid entalphy of feed [kJ/mol]

    def alphaEOSbut_init(m,i):
        return (1+m.mEOS[i]*(1-(m.TaliB/m.Tcrit[i])**(1/2)))**2

    m.alphaEOSbut = pe.Param(m.I, initialize=alphaEOSbut_init)

    def aiEOSbut_init(m,i):
        return m.alphaEOSbut[i]*0.42747*((0.00008314*m.TcritSRK[i])**2)/m.Pcrit[i]

    m.aiEOSbut = pe.Param(m.I, initialize=aiEOSbut_init)

    m.aEOSbut = pe.Param(initialize=(sum((m.zb[i]/100)*m.biEOS[i] for i in m.I)))
    m.bEOSbut = pe.Param(initialize=(sum(sum((m.zb[i]/100)*(m.zb[i2]/100)*(m.aiEOSbut[i]*m.aiEOSbut[i2])**0.5 for i2 in m.I) for i in m.I)))
    m.Zbut = pe.Var(m.N, within=pe.NonNegativeReals)

    @m.Constraint(m.N)
    def VaporZbut(m,n):
        if n != NT and n != 1:
            return (m.Zbut[n])**3-(m.Zbut[n])**2+(m.Zbut[n])*((m.aEOSbut*m.P[n]/((0.00008314*m.TaliB)**2))-(m.bEOSbut*m.P[n]/(0.00008314*m.TaliB))-(m.bEOSbut*m.P[n]/(0.00008314*m.TaliB))**2)-((m.aEOSbut*m.P[n]/((0.00008314*m.TaliB)**2)))*(m.bEOSbut*m.P[n]/(0.00008314*m.TaliB)) == 0
        else:
            return pe.Constraint.Skip

    m.HFB = pe.Var(m.N, within=pe.Reals)    # Butene entalpthy in feed
    @m.Constraint(m.N)
    def EqHFB(m,n):
        if n != NT and n != 1:
            return m.HFB[n] == sum((m.zb[i]/100)*(m.HL_b[i]+(8.314/1000)*m.TaliB*(m.Zbut[n]-1)+(1+m.mEOS[i])*((m.aEOSbut**0.5)/m.bEOSbut)*pe.log(m.Zbut[n]/(m.Zbut[n]+(m.bEOSbut*m.P[n]/(0.00008314*m.TaliB))))) for i in m.I)
        else:
            return pe.Constraint.Skip

    def HV_e_init(m,i):
        return ( (m.C1c[i]*(m.TaliE-m.Tref)) + ((m.C2c[i]/2)*((m.TaliE**2)-(m.Tref**2))) + ((m.C3c[i]/3)*((m.TaliE**3)-(m.Tref**3))) + ((m.C4c[i]/4)*((m.TaliE**4)-(m.Tref**4))) + ((m.C5c[i]/5)*((m.TaliE**5)-(m.Tref**5))) + ((m.C6c[i]/6)*((m.TaliE**6)-(m.Tref**6))) + m.Hform[i])

    m.HV_e = pe.Param(m.I, initialize=HV_e_init)    # Vapor entalphy of feed [kJ/mol]

    def Tred_e_init(m,i):
        return m.TaliE/m.Tcrit[i]

    m.Tred_e = pe.Param(m.I, initialize=Tred_e_init)    # Reduced temperature of feed [K]

    def DHVap_e_init(m,i):
        return ( m.C1v[i]*( (1-m.Tred_e[i])**( m.C2v[i] + (m.C3v[i]*m.Tred_e[i]) + (m.C4v[i]*(m.Tred_e[i]**2)) + (m.C5v[i]*(m.Tred_e[i]**3)) ) ) )

    m.DHVap_e = pe.Param(m.I, initialize=DHVap_e_init)    # Vaporization entalphy of feed [kJ/mol]

    def HL_e_init(m,i):
        return m.HV_e[i]-m.DHVap_e[i]

    m.HL_e = pe.Param(m.I, initialize=HL_e_init)    # Liquid entalphy of feed [kJ/mol]

    def alphaEOSeth_init(m,i):
        return (1+m.mEOS[i]*(1-(m.TaliE/m.Tcrit[i])**(1/2)))**2

    m.alphaEOSeth = pe.Param(m.I, initialize=alphaEOSeth_init)

    def aiEOSeth_init(m,i):
        return m.alphaEOSeth[i]*0.42747*((0.00008314*m.TcritSRK[i])**2)/m.Pcrit[i]

    m.aiEOSeth = pe.Param(m.I, initialize=aiEOSeth_init)
    
    m.aEOSeth = pe.Param(initialize=(sum(sum((m.ze[i]/100)*(m.ze[i2]/100)*(m.aiEOSeth[i]*m.aiEOSeth[i2])**0.5 for i2 in m.I) for i in m.I)))
    m.bEOSeth = pe.Param(initialize=(sum((m.ze[i]/100)*m.biEOS[i] for i in m.I)))
    m.Zeth = pe.Var(m.N, within=pe.NonNegativeReals)
    @m.Constraint(m.N)
    def VaporZeth(m,n):
        if n != NT and n != 1:
            return (m.Zeth[n])**3-(m.Zeth[n])**2+(m.Zeth[n])*((m.aEOSeth*m.P[n]/((0.00008314*m.TaliE)**2))-(m.bEOSeth*m.P[n]/(0.00008314*m.TaliE))-(m.bEOSeth*m.P[n]/(0.00008314*m.TaliE))**2)-((m.aEOSeth*m.P[n]/((0.00008314*m.TaliE)**2)))*(m.bEOSeth*m.P[n]/(0.00008314*m.TaliE)) == 0
        else:
            return pe.Constraint.Skip

    m.HFE = pe.Var(m.N, within=pe.Reals)    # Ethanol entalpthy in feed
    @m.Constraint(m.N)
    def EqHFE(m,n):
        if n != NT and n != 1:
            return m.HFE[n] == sum((m.ze[i]/100)*(m.HL_e[i]+(8.314/1000)*m.TaliE*(m.Zeth[n]-1)+(1+m.mEOS[i])*((m.aEOSeth**0.5)/m.bEOSeth)*pe.log(m.Zeth[n]/(m.Zeth[n]+(m.bEOSeth*m.P[n]/(0.00008314*m.TaliE))))) for i in m.I)
        else:
            return pe.Constraint.Skip

    # ______________________________ Section 13 ______________________________
    # Parameter, constraint and binary variable definition

    m.CASE = pe.Param(initialize=0) # 1 if there is equilibrium in reactive stages, 0 otherwise
    m.yc = pe.Var(m.N, within=pe.Binary)    # 1 if in stage n there is catalizer, 0 otherwise
    m.yr = pe.Var(m.N, within=pe.Binary)    # 1 if in stage n there is reflux, 0 otherwise  (parameter?)
    m.yb = pe.Var(m.N, within=pe.Binary)    # 1 if in stage n there is boilup, 0 otherwise
    m.par = pe.Var(m.N, within=pe.NonNegativeReals)    # 1 if stage n exists, 0 otherwise   (real?)

    @m.Constraint(m.N)
    def eqpar(m,n):
        if n == 1:
            return m.par[n] == 1
        elif n == NT:
            return m.par[n] == 1
        else:
            return m.par[n] == 1 - (1-sum(m.yr[j] for j in range(2,n+1)) - (sum(m.yb[j] for j in range(2,n))))

    m.ye = pe.Var(m.N, within=pe.NonNegativeReals)    # 1 if in stage n there is equilibrium, 0 otherwise   (real?)

    @m.Constraint(m.N)
    def eqyeq(m,n):
        if n != NT and n != 1:
            return m.ye[n] == m.par[n]*(1-m.yc[n])*m.CASE+m.par[n]*(1-m.CASE)
        else:
            return pe.Constraint.Skip

    return m

    m.yf = pe.Var(m.N, m.F, within=pe.Binary)    # 1 if in stage n there is feed f, 0 otherwise

    # ______________________________ Section 14 ______________________________
    # Logic constraints

    m.cmej = pe.Param(initialize=1)
    m.NCmax = pe.Param(initialize=KT)   # Maximum number of reactive stages

    @m.Constraint(m.N)
    def logic1(m,n):    # The boilup is below the reflux stage
        if n != NT and n != 1:
            return m.cmej*sum(m.yr[j] for j in range(2,n+1)) >= m.yb[n]*m.cmej
        else:
            return pe.Constraint.Skip

    @m.Constraint()
    def logic2(m):    # There is only one reflux stage
        return m.cmej*sum(m.yr[j] for j in range(2,NT)) == 1*m.cmej

    @m.Constraint()
    def logic3(m):    # There is only one boil-up stage
        return m.cmej*sum(m.yb[j] for j in range(2,NT)) == 1*m.cmej

    @m.Constraint(m.F)
    def logic4(m,f):    # There is only one stage per feed
        return m.cmej*sum(m.yf[j,f] for j in range(2,NT)) == 1*m.cmej

    @m.Constraint()
    def logic6(m):    # There is maximum number of catalitic stages
        return m.cmej*sum(m.yc[j] for j in range(2,NT)) == m.NCmax*m.cmej

    @m.Constraint(m.N, m.F)
    def logic7(m,n,f):    # The feed stages are below reflux stage
        if n != NT and n != 1:
            return m.cmej*sum(m.yr[j] for j in range(2,n+1)) >= m.yf[n,f]*m.cmej
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.N, m.F)
    def logic8(m,n,f):    # The boil-up stage is below feed stages
        if n != NT and n != 1:
            return m.cmej*sum(m.yf[j,f] for j in range(2,n+1)) >= m.yb[n]*m.cmej
        else:
            return pe.Constraint.Skip
    
    @m.Constraint(m.N)
    def logic9(m,n):    # The ethanol feed is above the butenes feed
        if n != NT and n != 1:
            return m.cmej*sum(m.yf[j,1] for j in range(2,n+1)) >= m.yf[n,2]*m.cmej
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.N)
    def logic10(m,n):    # The catalitic stages are below the ethanol feed stage
        if n != NT and n != 1:
            return m.cmej*sum(m.yf[j,1] for j in range(2,n+1)) >= m.yc[n]*m.cmej
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.N)
    def logic11(m,n):    # The catalitic stages are above the butenes feed stage
        if n != NT and n != 1:
            return m.cmej*(sum(m.yf[j,2] for j in range(2,n+1)) - m.yf[n,2] )<= m.cmej*(1-m.yc[n])
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.N)
    def logic12(m,n):    # The catalitic stages are below reflux stage
        if n != NT and n != 1:
            return m.cmej*sum(m.yr[j] for j in range(2,n+1)) >= m.yc[n]*m.cmej
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.N)
    def logic13(m,n):    # The catalitic stages are above boil-up stage
        if n != NT and n != 1:
            return m.cmej*(sum(m.yb[j] for j in range(2,n+1)) - m.yb[n] )<= m.cmej*(1-m.yc[n])
        else:
            return pe.Constraint.Skip

    # ______________________________ Section 15 ______________________________
    # Condenser equations

    # Steady-state equations
    @m.Constraint()
    def BalMassC0(m):
        return 0 == m.V[2]-m.V[1]*(1+m.RR)

    @m.Constraint(m.I)
    def BalMassPartialC0(m,i):
        return 0 == m.V[2]*m.y[i,2]-m.V[1]*m.x[i,1]*(1+m.RR)

    @m.Constraint()
    def SumC0(m):
        return sum(m.y[i,1]-m.x[i,1] for i in m.I) == 0
        
    @m.Constraint(m.I)
    def EquilibriumC0(m,i):
        return m.y[i,1]*m.P[1]*m.phi[i] == m.Psat[i,1]*m.gamma[i,1]*m.x[i,1]

    @m.Constraint()
    def BalEnergyC0(m):
        return 0 == m.V[2]*m.HV[2]-m.V[1]*(1+m.RR)*m.HL[1]-m.QC

    # Fixed liquid flow 
    @m.Constraint()
    def fixedL(m):
        return m.L[1] == 0

    
    

    
    

    

    return m


if __name__ == "__main__":
    NT = 22
    results = minlp_catalitic_column(NT)