import pyomo.environ as pe
import networkx as nx
import matplotlib.pyplot as plt
from math import (sqrt, pi)
from pyomo.core.base.misc import display
from pyomo.opt.base.solvers import SolverFactory
import os

def minlp_extractive_column(NT=30,  visualize=False):

    # PYOMO MODEL
    m = pe.ConcreteModel(name='minlp_extractive_column')

    # ______________________________ Sections 1-4 ______________________________
    # Main variables and sets for defining the system
    # Hydraulic equations calculation

    # Sets
    m.I = pe.Set(initialize=['Water', 'Ethanol', 'Glycerol'])  # Set of components
    m.N = pe.RangeSet(1, NT)  # Set of all stages in the column

    # Variables
    m.L = pe.Var(m.N,  within=pe.NonNegativeReals, bounds=(0, 200))  # Flow of liquid [mol/hr]
    m.V = pe.Var(m.N,  within=pe.NonNegativeReals, bounds=(0, 200))  # Flow of vapor [mol/hr]
    m.x = pe.Var(m.I, m.N,  within=pe.NonNegativeReals, bounds=(0, 100))  # Molar composition of liquid [*]
    m.y = pe.Var(m.I, m.N,  within=pe.NonNegativeReals, bounds=(0, 100))  # Molar composition of vapor [*]
    m.Temp = pe.Var(m.N,  within=pe.NonNegativeReals, bounds=(200, 417.89))   # Operation temperature [K]
    m.P = pe.Var(m.N,  within=pe.NonNegativeReals)  # Stage pressure [atm]
    m.Z = pe.Var(m.N,  within=pe.NonNegativeReals, bounds=(0.5, 1.3))  # Compressibility coefficient [*]
    m.RR = pe.Var(within=pe.NonNegativeReals, bounds=(1, 10))   # Reflux ratio [*]
    m.Qc = pe.Var(within=pe.NonNegativeReals, bounds=(0, 900))   # Condensator duty [kJ/hr]
    m.Qr = pe.Var(within=pe.NonNegativeReals, bounds=(100, 400))   # Reboiler duty [kJ/hr]
    m.BR = pe.Var(within=pe.NonNegativeReals, bounds=(0, 10))   # Boil up [*]

    # Hydraulic parameters
    m.d_hole = pe.Param(initialize=0.0127)   # Hole diameter [m]
    m.tray_t = pe.Param(initialize=0.002)   # Plate thickness [m]
    m.hw = pe.Param(initialize=0.0254)   # Weir height [m]
    m.Lw = pe.Param(initialize=0.578)   # Weir lenght [m]
    m.HS = pe.Param(initialize=0.61)   # Plate height [m]
    m.Sfactor = pe.Param(initialize=0.15)   # Safety height factor
    m.K0 = pe.Param(initialize=(880.6-(67.7*m.d_hole/m.tray_t)+(7.32*((m.d_hole/m.tray_t)**2))-(0.338*((m.d_hole/m.tray_t)**3)))*10**-3) # Hole coefficient [*]

    # Hydraulic variables
    m.D = pe.Var(within=pe.NonNegativeReals, bounds=(1.23, 2))   # Column diameter [m]
    m.Htotal = pe.Var(within=pe.NonNegativeReals, bounds=(0, 10))   # Total column height [m]
    m.At = pe.Var(within=pe.NonNegativeReals, bounds=(0.299, 0.5))   # Active area [m**2]
    m.Ad = pe.Var(within=pe.NonNegativeReals, bounds=(0.037, 0.05))   # Weir area [m**2]
    m.A0 = pe.Var(within=pe.NonNegativeReals, bounds=(0.012, 0.02))   # Holed area [m**2]
    m.poro = pe.Var(within=pe.NonNegativeReals, bounds=(0.907*pe.sqrt(m.d_hole/(0.12*0.05)), 1))   # Plate porosity [*]
    m.pitch = pe.Var(within=pe.NonNegativeReals, bounds=(0.072, 0.1))   # Distance between plate holes [m]
    m.A_col = pe.Var(within=pe.NonNegativeReals)   # Cross section area [m**2]

    # Hydraulic constraints
    @m.Constraint()
    def EqAt(m):
        return m.At == pe.sqrt(m.D/2)*(pi-(1.854590-0.96))
    
    @m.Constraint()
    def EqAd(m):
        return m.Ad == pe.sqrt(m.D/2)*(0.5*(1.854590-0.96))

    @m.Constraint()
    def EqPitch(m):
        return m.pitch == pe.sqrt(m.D/2)*pi*0.12

    @m.Constraint()
    def EqA0(m):
        return m.A0 == m.At*m.poro

    @m.Constraint()
    def EqPoro(m):
        return m.poro == 0.907*pe.sqrt(m.d_hole/m.pitch)

    @m.Constraint()
    def EqDiam(m):
        return m.A_col == pe.sqrt(m.D/2)*pi

    # Azeotropic feed
    m.Faz = pe.Param(initialize=100)   # Azeotrope feed [kmol/hr]
    zaz_Eth = 85
    zAz_init = {'Water':100-zaz_Eth, 'Ethanol':zaz_Eth, 'Glycerol':0}
    m.zAz = pe.Param(m.I, initialize=zAz_init)  # Molar composition

    # Glycerol feed
    m.FG = pe.Param(initialize=52)   # Glycerol feed [kmol/hr]
    zg_init = {'Water':0, 'Ethanol':0, 'Glycerol':100}
    m.zg = pe.Param(m.I, initialize=zg_init)

    # Operation parameters
    m.Pop = pe.Param(initialize=1)   # Pressure at condenser [atm]
    m.Tali1 = pe.Param(initialize=305)   # Glycerol feed temperature [K]
    m.Tali2 = pe.Param(initialize=293.15)   # Azeotrope feed temperature [K]
    m.xReth = pe.Param(initialize=99.5)   # Desired composition of Ethanol in bottom [*]
    m.R = pe.Param(initialize=8.31447215)   # Ideal gas constant [kJ/K*kmol]

    # Composition restriction in bottoms
    @m.Constraint()
    def pureza0(m):
        return m.x['Ethanol',1] >= m.xReth

    # ______________________________ Section 5 ______________________________
    # Saturation pressures usinn Antoine equation

    # Constants for expanded Antoine equation
    C1a_init = {'Water':62.12291155, 'Ethanol':61.77791155, 'Glycerol':88.45991155}
    m.C1a = pe.Param(m.I, initialize=C1a_init)
    C2a_init = {'Water':-7258.200000, 'Ethanol':-7122.300000, 'Glycerol':-13808.00000}
    m.C2a = pe.Param(m.I, initialize=C2a_init)
    C3a_init = {'Water':0, 'Ethanol':0, 'Glycerol':0}
    m.C3a = pe.Param(m.I, initialize=C3a_init)
    C4a_init = {'Water':0, 'Ethanol':0, 'Glycerol':0}
    m.C4a = pe.Param(m.I, initialize=C4a_init)
    C5a_init = {'Water':-7.303700000, 'Ethanol':-7.142400000, 'Glycerol':-10.08800000}
    m.C5a = pe.Param(m.I, initialize=C5a_init)
    C6a_init = {'Water':4.16530000*10**-6, 'Ethanol':2.88530000*10**-6, 'Glycerol':3.5712000*10**-19}
    m.C6a = pe.Param(m.I, initialize=C6a_init)
    C7a_init = {'Water':2, 'Ethanol':2, 'Glycerol':6}
    m.C7a = pe.Param(m.I, initialize=C7a_init)

    # Antoine equation
    m.Psat = pe.Var(m.I, m.N, within=pe.NonNegativeReals, bounds=(0, 100))   # Saturation pressure [atm]
    @m.Constraint(m.I, m.N)
    def EqPsat(m,i,n):
        return m.Psat[i,n] == 1/1.01325*pe.exp(m.C1a[i] + (m.C2a[i]/(m.Temp[n]+m.C3a[i])) + (m.C4a[i]*m.Temp[n]) + (m.C5a[i]*pe.log(m.Temp[n]) + (m.C6a[i]*(m.Temp[n]**m.C7a[i]))) )

    # ______________________________ Section 6 ______________________________
    # Calculation of liquid density using IK-CAPI equation
    # Calculation of liquid density using critic DIPPR equation
    # Calculation of gas density using corrected ideal gas equation

    # Constants for DIPPR equation
    MW_init = {'Water':18.01528, 'Ethanol':46.06904, 'Glycerol':92.09382}
    m.MW = pe.Param(m.I, initialize=MW_init)    # Molecular weight [kg/kmol]
    Tcrit_init = {'Water':647.096, 'Ethanol':514.00, 'Glycerol':850.00}
    m.Tcrit = pe.Param(m.I, initialize=Tcrit_init) # Critic temperature [K]
    Pcrit_init = {'Water':220.6351, 'Ethanol':63.00, 'Glycerol':75.00}
    m.Pcrit = pe.Param(m.I, initialize=Pcrit_init) # Critic pressure [bar]
    C1r_init = {'Water':17.863, 'Ethanol':1.6288, 'Glycerol':0.92382}
    m.C1r = pe.Param(m.I, initialize=C1r_init)
    C2r_init = {'Water':58.606, 'Ethanol':0.27469, 'Glycerol':0.24386}
    m.C2r = pe.Param(m.I, initialize=C2r_init)
    C3r_init = {'Water':-95.396, 'Ethanol':515, 'Glycerol':850}
    m.C3r = pe.Param(m.I, initialize=C3r_init)
    C4r_init = {'Water':213.89, 'Ethanol':0.23178, 'Glycerol':0.22114}
    m.C4r = pe.Param(m.I, initialize=C4r_init)
    m.C5r = pe.Param(initialize=-141.26)

    m.Tcritm = pe.Var(m.N, within=pe.NonNegativeReals, bounds=(417.9, 600))
    @m.Constraint(m.N)
    def EqTcritm(m,n):
        return m.Tcritm[n] == (pe.sqrt(sum((m.x[i,n]/100)*m.Tcrit[i]/(m.Pcrit[i]**0.5) for i in m.I)))/(sum((m.x[i,n]/100)*m.Tcrit[i]/m.Pcrit[i] for i in m.I))
    
    m.rho = pe.Var(m.I, m.N, within=pe.NonNegativeReals, bounds=(15000, 10000)) # Liquid molar density [mol/m**3]
    @m.Constraint(m.I, m.N)
    def Eqrho12(m,i,n):
        if i != 'Water':
            return m.rho[i,n] == ( m.C1r[i]/(m.C2r[i]**(1+((1-(m.Temp[n]/m.Tcritm[n]))**m.C4r[i]))) )*1000
        else:
            return m.rho[i,n] == ( m.C1r[i]+(m.C2r[i]*(1-(m.Temp[n]/m.Tcritm[n]))**(0.35))+(m.C3r[i]*(1-(m.Temp[n]/m.Tcritm[n]))**(2/3))+(m.C4r[i]*(1-(m.Temp[n]/m.Tcritm[n]))) + (m.C5r*(1-(m.Temp[n]/m.Tcritm[n]))**(4/3)) )*1000

    m.rhoV = pe.Var(m.N, within=pe.NonNegativeReals, bounds=(60, 500)) # Vapor molar density [mol/m**3]
    @m.Constraint(m.N)
    def EqrhoV(m,n):
        return m.rhoV[n] == m.P[n]/(m.R/101325*m.Temp[n])

    m.Qliq = pe.Var(m.N, within=pe.NonNegativeReals, bounds=(6, 100)) # Liquid volumetric flow [m**3/hr]
    @m.Constraint(m.N)
    def EqQliq(m,n):
        return m.Qliq[n] == m.L[n]/sum(m.rho[i,n]*m.x[i,n]/100 for i in m.I)

    m.Qvap = pe.Var(m.N, within=pe.NonNegativeReals, bounds=(4100, 50000)) # Vapor volumetric flow [m**3/hr]
    @m.Constraint(m.N)
    def EqQvap(m,n):
        return m.Qvap[n] == m.V[n]/m.rhoV[n]

    
    


    return m




if __name__ == "__main__":
    NT = 30
    m = minlp_extractive_column(NT)
