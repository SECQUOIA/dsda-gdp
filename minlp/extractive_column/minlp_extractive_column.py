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

    

    return m




if __name__ == "__main__":
    NT = 30
    m = minlp_extractive_column(NT)
