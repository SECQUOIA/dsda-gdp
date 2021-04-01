import pyomo.environ as pe
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from math import (sqrt, pi, exp, log)
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
    # Set of components
    m.I = pe.Set(initialize=['Water', 'Ethanol', 'Glycerol'])
    m.N = pe.RangeSet(1, NT)  # Set of all stages in the column

    # Operation pressure
    m.Pop = pe.Param(initialize=1)   # Pressure at condenser [atm]

    # Dependences
    # D / Temp / x / y / L / V

    x_init = {}
    y_init = {}
    Temp_init = {}

    Temp_init[1]=351.39
    Temp_init[2]=368.92
    Temp_init[3]=369.42
    Temp_init[4]=370.12
    Temp_init[5]=371.02
    Temp_init[6]=372.08
    Temp_init[7]=373
    Temp_init[8]=374
    Temp_init[9]=375
    Temp_init[10]=376
    Temp_init[11]=376
    Temp_init[12]=374
    Temp_init[13]=357
    Temp_init[14]=358
    Temp_init[15]=362
    Temp_init[16]=376
    Temp_init[30]=446.39

    x_init[('Water', 1)]=0.4862
    x_init[('Water', 2)]=1.03
    x_init[('Water', 3)]=2.26
    x_init[('Water', 4)]=3.86
    x_init[('Water', 5)]=5.78
    x_init[('Water', 6)]=7.87
    x_init[('Water', 7)]=9.89
    x_init[('Water', 8)]=11.68
    x_init[('Water', 9)]=13.11
    x_init[('Water', 10)]=14.19
    x_init[('Water', 11)]=15.05
    x_init[('Water', 12)]=16.40
    x_init[('Water', 13)]=19.69
    x_init[('Water', 14)]=30.62
    x_init[('Water', 15)]=53.69
    x_init[('Water', 16)]=68.51
    x_init[('Water', 30)]=21.88

    x_init[('Ethanol', 1)]=99.5
    x_init[('Ethanol', 2)]=30.12
    x_init[('Ethanol', 3)]=28.91
    x_init[('Ethanol', 4)]=27.33
    x_init[('Ethanol', 5)]=25.47
    x_init[('Ethanol', 6)]=23.48
    x_init[('Ethanol', 7)]=21.58
    x_init[('Ethanol', 8)]=19.95
    x_init[('Ethanol', 9)]=18.67
    x_init[('Ethanol', 10)]=17.76
    x_init[('Ethanol', 11)]=17.31
    x_init[('Ethanol', 12)]=18.48
    x_init[('Ethanol', 13)]=53.61
    x_init[('Ethanol', 14)]=42.56
    x_init[('Ethanol', 15)]=19.42
    x_init[('Ethanol', 16)]=3.44
    x_init[('Ethanol', 30)]=0.11
    
    x_init[('Glycerol', 1)]=0.01
    x_init[('Glycerol', 2)]=68.83
    x_init[('Glycerol', 3)]=68.82
    x_init[('Glycerol', 4)]=68.79
    x_init[('Glycerol', 5)]=68.74
    x_init[('Glycerol', 6)]=68.64
    x_init[('Glycerol', 7)]=68.51
    x_init[('Glycerol', 8)]=68.36
    x_init[('Glycerol', 9)]=68.2
    x_init[('Glycerol', 10)]=68.04
    x_init[('Glycerol', 11)]=67.62
    x_init[('Glycerol', 12)]=65.11
    x_init[('Glycerol', 13)]=26.68
    x_init[('Glycerol', 14)]=26.81
    x_init[('Glycerol', 15)]=26.87
    x_init[('Glycerol', 16)]=28.04
    x_init[('Glycerol', 30)]=77.99

    for n in range(17,30):
        x_init[('Water', n)]=68.51
        x_init[('Ethanol', n)]=3.44
        x_init[('Glycerol', n)]=28.04
        Temp_init[n]=376.2923


    for n in m.N:
        #x_init[('Water', n)] = 1.03159004
        y_init[('Water', n)] = 0.48
        #x_init[('Ethanol', n)] = 30.1285211
        y_init[('Ethanol', n)] = 99.5
        #x_init[('Glycerol', n)] = 68.8398887
        y_init[('Glycerol', n)] = 0.01

    # Variables
    m.L = pe.Var(m.N,  within=pe.NonNegativeReals, bounds=(
        0, 10**6), initialize=75500)  # Flow of liquid [mol/hr]
    m.V = pe.Var(m.N,  within=pe.NonNegativeReals, bounds=(
        0, 10**6), initialize=90300)  # Flow of vapor [mol/hr]
    m.x = pe.Var(m.I, m.N,  within=pe.NonNegativeReals, bounds=(
        0, 100), initialize=x_init)  # Molar composition of liquid [*]
    m.y = pe.Var(m.I, m.N,  within=pe.NonNegativeReals, bounds=(
        0, 100), initialize=y_init)  # Molar composition of vapor [*]
    m.Temp = pe.Var(m.N,  within=pe.NonNegativeReals, bounds=(
        300, 500), initialize=368.92)   # Operation temperature [K]
    m.P = pe.Var(m.N,  within=pe.NonNegativeReals, bounds=(
        m.Pop, 2), initialize=1)  # Stage pressure [atm]
    # Compressibility coefficient [*]
    m.Z = pe.Var(m.N,  within=pe.NonNegativeReals)
    m.RR = pe.Var(within=pe.NonNegativeReals, bounds=(
        10**-4, 10), initialize=0.058)   # Reflux ratio [*]
    m.Qc = pe.Var(within=pe.NonNegativeReals, bounds=(
        6*10**4, 2*10**8), initialize=3635108)   # Condensator duty [kJ/hr]
    m.Qr = pe.Var(within=pe.NonNegativeReals, bounds=(
        10**5, 10**8), initialize=6144494)   # Reboiler duty [kJ/hr]
    m.BR = pe.Var(within=pe.NonNegativeReals, bounds=(
        0, 20), initialize=2)   # Boil up [*]

    # Hydraulic parameters
    m.d_hole = pe.Param(initialize=0.0127)   # Hole diameter [m]
    m.tray_t = pe.Param(initialize=0.002)   # Plate thickness [m]
    m.hw = pe.Param(initialize=0.0254)   # Weir height [m]
    m.Lw = pe.Param(initialize=0.578)   # Weir lenght [m]
    m.HS = pe.Param(initialize=0.61)   # Plate height [m]
    m.Sfactor = pe.Param(initialize=0.15)   # Safety height factor
    m.K0 = pe.Param(initialize=(880.6-(67.7*m.d_hole/m.tray_t)+(7.32*((m.d_hole/m.tray_t)
                                                                      ** 2))-(0.338*((m.d_hole/m.tray_t)**3)))*10**-3)  # Hole coefficient [*]

    # Hydraulic variables
    m.D = pe.Var(within=pe.NonNegativeReals, bounds=(0.01, 10), initialize=0.91)   # Column diameter [m]
    m.pitch = pe.Var(within=pe.NonNegativeReals, initialize=0.07943)   # Distance between plate holes [m]
    m.At = pe.Var(within=pe.NonNegativeReals, initialize=0.4735)   # Active area [m**2]
    m.Ad = pe.Var(within=pe.NonNegativeReals, initialize=0.0942)   # Weir area [m**2]
    m.poro = pe.Var(within=pe.NonNegativeReals, initialize=0.02317)   # Plate porosity [*]
    m.A0 = pe.Var(within=pe.NonNegativeReals, initialize=pe.value(m.At)*pe.value(m.poro))   # Holed area [m**2]
    m.A_col = pe.Var(within=pe.NonNegativeReals, initialize=2) # Cross section area [m**2]

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
    m.FAz = pe.Param(initialize=100)   # Azeotrope feed [kmol/hr]
    zaz_Eth = 85
    zAz_init = {'Water': 100-zaz_Eth, 'Ethanol': zaz_Eth, 'Glycerol': 0}
    m.zAz = pe.Param(m.I, initialize=zAz_init)  # Molar composition

    # Glycerol feed
    m.FG = pe.Param(initialize=52)   # Glycerol feed [kmol/hr]
    zg_init = {'Water': 0, 'Ethanol': 0, 'Glycerol': 100}
    m.zg = pe.Param(m.I, initialize=zg_init)

    # Operation parameters
    m.Tali1 = pe.Param(initialize=305)   # Glycerol feed temperature [K]
    m.Tali2 = pe.Param(initialize=293.15)   # Azeotrope feed temperature [K]
    # Desired composition of Ethanol in bottom [*]
    m.xReth = pe.Param(initialize=99.5)
    m.R = pe.Param(initialize=8.31447215)   # Ideal gas constant [kJ/K*kmol]

    # Composition restriction in bottoms
    @m.Constraint()
    def pureza0(m):
        return m.x['Ethanol', 1] >= m.xReth

    # ______________________________ Section 5 ______________________________
    # Saturation pressures usinn Antoine equation

    # Constants for expanded Antoine equation
    C1a_init = {'Water': 62.12291155,
                'Ethanol': 61.77791155, 'Glycerol': 88.45991155}
    m.C1a = pe.Param(m.I, initialize=C1a_init)
    C2a_init = {'Water': -7258.200000, 'Ethanol': -
                7122.300000, 'Glycerol': -13808.00000}
    m.C2a = pe.Param(m.I, initialize=C2a_init)
    C3a_init = {'Water': 0, 'Ethanol': 0, 'Glycerol': 0}
    m.C3a = pe.Param(m.I, initialize=C3a_init)
    C4a_init = {'Water': 0, 'Ethanol': 0, 'Glycerol': 0}
    m.C4a = pe.Param(m.I, initialize=C4a_init)
    C5a_init = {'Water': -7.303700000, 'Ethanol': -
                7.142400000, 'Glycerol': -10.08800000}
    m.C5a = pe.Param(m.I, initialize=C5a_init)
    C6a_init = {'Water': 4.16530000*10**-6,
                'Ethanol': 2.88530000*10**-6, 'Glycerol': 3.5712000*10**-19}
    m.C6a = pe.Param(m.I, initialize=C6a_init)
    C7a_init = {'Water': 2, 'Ethanol': 2, 'Glycerol': 6}
    m.C7a = pe.Param(m.I, initialize=C7a_init)

    # Antoine equation
    Psat_init = {}
    for i in m.I:
        for n in m.N:
            Psat_init[(i, n)] = 1/1.01325*exp(m.C1a[i] + (m.C2a[i]/(pe.value(m.Temp[n])+m.C3a[i])) + (m.C4a[i] *
                                                                                                      pe.value(m.Temp[n])) + (m.C5a[i]*log(pe.value(m.Temp[n])) + (m.C6a[i]*(pe.value(m.Temp[n])**m.C7a[i]))))

    m.Psat = pe.Var(m.I, m.N, within=pe.NonNegativeReals,
                    initialize=Psat_init)   # Saturation pressure [atm]

    @m.Constraint(m.I, m.N)
    def EqPsat(m, i, n):
        return m.Psat[i, n] == 1/1.01325*pe.exp(m.C1a[i] + (m.C2a[i]/(m.Temp[n]+m.C3a[i])) + (m.C4a[i]*m.Temp[n]) + (m.C5a[i]*pe.log(m.Temp[n]) + (m.C6a[i]*(m.Temp[n]**m.C7a[i]))))

    # ______________________________ Section 6 ______________________________
    # Calculation of liquid density using IK-CAPI equation
    # Calculation of liquid density using critic DIPPR equation
    # Calculation of gas density using corrected ideal gas equation

    # Constants for DIPPR equation
    MW_init = {'Water': 18.01528, 'Ethanol': 46.06904, 'Glycerol': 92.09382}
    m.MW = pe.Param(m.I, initialize=MW_init)    # Molecular weight [kg/kmol]
    Tcrit_init = {'Water': 647.096, 'Ethanol': 514.00, 'Glycerol': 850.00}
    m.Tcrit = pe.Param(m.I, initialize=Tcrit_init)  # Critic temperature [K]
    Pcrit_init = {'Water': 220.6351, 'Ethanol': 63.00, 'Glycerol': 75.00}
    m.Pcrit = pe.Param(m.I, initialize=Pcrit_init)  # Critic pressure [bar]
    C1r_init = {'Water': 17.863, 'Ethanol': 1.6288, 'Glycerol': 0.92382}
    m.C1r = pe.Param(m.I, initialize=C1r_init)
    C2r_init = {'Water': 58.606, 'Ethanol': 0.27469, 'Glycerol': 0.24386}
    m.C2r = pe.Param(m.I, initialize=C2r_init)
    C3r_init = {'Water': -95.396, 'Ethanol': 515, 'Glycerol': 850}
    m.C3r = pe.Param(m.I, initialize=C3r_init)
    C4r_init = {'Water': 213.89, 'Ethanol': 0.23178, 'Glycerol': 0.22114}
    m.C4r = pe.Param(m.I, initialize=C4r_init)
    m.C5r = pe.Param(initialize=-141.26)

    Tcritm_init = {}
    for n in m.N:
        Tcritm_init[n] = sum((pe.value(m.x[i, n])/100)*m.Tcrit[i] for i in m.I)

    m.Tcritm = pe.Var(m.N, within=pe.NonNegativeReals, bounds=(
        510, 860), initialize=Tcritm_init)  # 745

    @m.Constraint(m.N)
    def EqTcritm(m, n):
        return m.Tcritm[n] == (sum((m.x[i, n]/100)*m.Tcrit[i] for i in m.I))

    rho_init = {}
    for n in m.N:
        for i in m.I:
            if i != 'Water':
                rho_init[(i, n)] = (
                    (m.C1r[i]/(m.C2r[i]**(1+((1-(pe.value(m.Temp[n])/pe.value(m.Tcritm[n])))**m.C4r[i]))))*1000)
            else:
                rho_init[(i, n)] = ((m.C1r[i]+(m.C2r[i]*(1-(pe.value(m.Temp[n])/pe.value(m.Tcritm[n])))**(0.35))+(m.C3r[i]*(1-(pe.value(m.Temp[n])/pe.value(m.Tcritm[n])))
                                                                                                                  ** (2/3))+(m.C4r[i]*(1-(pe.value(m.Temp[n])/pe.value(m.Tcritm[n])))) + (pe.value(m.C5r)*(1-(pe.value(m.Temp[n])/pe.value(m.Tcritm[n])))**(4/3)))*1000)

    m.rho = pe.Var(m.I, m.N, within=pe.NonNegativeReals, bounds=(
        1000, 70000), initialize=rho_init)  # Liquid molar density [mol/m**3]

    @m.Constraint(m.I, m.N)
    def Eqrho12(m, i, n):
        if i != 'Water':
            return m.rho[i, n] == (m.C1r[i]/(m.C2r[i]**(1+((1-(m.Temp[n]/m.Tcritm[n]))**m.C4r[i]))))*1000
        else:
            return m.rho[i, n] == (m.C1r[i]+(m.C2r[i]*(1-(m.Temp[n]/m.Tcritm[n]))**(0.35))+(m.C3r[i]*(1-(m.Temp[n]/m.Tcritm[n]))**(2/3))+(m.C4r[i]*(1-(m.Temp[n]/m.Tcritm[n]))) + (m.C5r*(1-(m.Temp[n]/m.Tcritm[n]))**(4/3)))*1000

    rhoV_init = {}
    for n in m.N:
        rhoV_init[n] = pe.value(m.P[n])/(m.R/101325*pe.value(m.Temp[n]))

    # Vapor molar density [mol/m**3]
    m.rhoV = pe.Var(m.N, within=pe.NonNegativeReals, initialize=rhoV_init)

    @m.Constraint(m.N)
    def EqrhoV(m, n):
        return m.rhoV[n] == m.P[n]/(m.R/101325*m.Temp[n])

    Qliq_init = {}
    for n in m.N:
        Qliq_init[n] = pe.value(
            m.L[n])/sum(pe.value(m.rho[i, n])*pe.value(m.x[i, n])/100 for i in m.I)

    # Liquid volumetric flow [m**3/hr]
    m.Qliq = pe.Var(m.N, within=pe.NonNegativeReals, initialize=Qliq_init)

    @m.Constraint(m.N)
    def EqQliq(m, n):
        return m.Qliq[n] == m.L[n]/sum(m.rho[i, n]*m.x[i, n]/100 for i in m.I)

    Qvap_init = {}
    for n in m.N:
        Qvap_init[n] = pe.value(m.V[n])/pe.value(m.rhoV[n])

    # Vapor volumetric flow [m**3/hr]
    m.Qvap = pe.Var(m.N, within=pe.NonNegativeReals, initialize=Qvap_init)

    @m.Constraint(m.N)
    def EqQvap(m, n):
        return m.Qvap[n] == m.V[n]/m.rhoV[n]

    # ______________________________ Section 7 ______________________________
    # Calculation of superficial tension using critic DIPPR equation

    # Constants for DIPPR equation
    C1sig_init = {'Water': 0.17766,
                  'Ethanol': 0.0241005, 'Glycerol': 0.0645335}
    m.C1sig = pe.Param(m.I, initialize=C1sig_init)
    C2sig_init = {'Water': 2.567, 'Ethanol': -
                  7.75658*10**-5, 'Glycerol': -5.38024*10**-5}
    m.C2sig = pe.Param(m.I, initialize=C2sig_init)
    C3sig_init = {'Water': -3.3377, 'Ethanol': -
                  1.025*10**-7, 'Glycerol': -2.1558*10**-7}
    m.C3sig = pe.Param(m.I, initialize=C3sig_init)
    C4sig_init = {'Water': 1.9699, 'Ethanol': 0, 'Glycerol': 0}
    m.C4sig = pe.Param(m.I, initialize=C4sig_init)

    sigma_init = {}
    for n in m.N:
        sigma_init[n] = sum((pe.value(m.x[i, n])/100)*m.C1sig[i]*(1-(pe.value(m.Temp[n])/pe.value(m.Tcritm[n])))**(m.C2sig[i]+m.C3sig[i]*(
            pe.value(m.Temp[n])/pe.value(m.Tcritm[n]))+m.C4sig[i]*((pe.value(m.Temp[n])/pe.value(m.Tcritm[n])))**2) for i in m.I)

    # Liquid-vapor superficial tension [N/m]
    m.sigma = pe.Var(m.N, within=pe.NonNegativeReals, initialize=sigma_init)

    @m.Constraint(m.N)
    def Eqsigma(m, n):
        return m.sigma[n] == sum((m.x[i, n]/100)*m.C1sig[i]*(1-(m.Temp[n]/m.Tcritm[n]))**(m.C2sig[i]+m.C3sig[i]*(m.Temp[n]/m.Tcritm[n])+m.C4sig[i]*((m.Temp[n]/m.Tcritm[n]))**2) for i in m.I)

    # ______________________________ Sections 8-9 ______________________________
    # Calculation of activity coefficient using NRTL model

    a_nrtl_init = {('Water', 'Water'): 0, ('Water', 'Ethanol'): 3.4578, ('Water', 'Glycerol'): -1.2515,
                   ('Ethanol', 'Water'): -0.8009, ('Ethanol', 'Ethanol'): 0, ('Ethanol', 'Glycerol'): 0,
                   ('Glycerol', 'Water'): -0.7318, ('Glycerol', 'Ethanol'): 0, ('Glycerol', 'Glycerol'): 0}
    m.a_nrtl = pe.Param(m.I, m.I, initialize=a_nrtl_init)
    b_nrtl_init = {('Water', 'Water'): 0, ('Water', 'Ethanol'): -586.0809, ('Water', 'Glycerol'): 272.6075,
                   ('Ethanol', 'Water'): 246.18, ('Ethanol', 'Ethanol'): 0, ('Ethanol', 'Glycerol'): 442.713,
                   ('Glycerol', 'Water'): 170.9167, ('Glycerol', 'Ethanol'): 36.139, ('Glycerol', 'Glycerol'): 0}
    m.b_nrtl = pe.Param(m.I, m.I, initialize=b_nrtl_init)
    c_nrtl_init = {(i, i2): 0.3 for i in m.I for i2 in m.I}
    for i in m.I:
        c_nrtl_init[i, i] = 0
    m.c_nrtl = pe.Param(m.I, m.I, initialize=c_nrtl_init)

    def alfa_nrtl_init(m, i, i2):
        if i != i2:
            return m.c_nrtl[i, i2]
        else:
            return pe.Param.Skip

    m.alfa_nrtl = pe.Param(m.I, m.I, initialize=alfa_nrtl_init, within=pe.Any)

    tao_nrtl_init = {}
    for n in m.N:
        for i in m.I:
            for i2 in m.I:
                if i == i2:
                    tao_nrtl_init[(i, i2, n)] = 0
                else:
                    tao_nrtl_init[(i, i2, n)] = m.a_nrtl[i, i2] + \
                        (m.b_nrtl[i, i2]/pe.value(m.Temp[n]))

    m.tao_nrtl = pe.Var(m.I, m.I, m.N, within=pe.Reals,
                        initialize=tao_nrtl_init)

    @m.Constraint(m.I, m.I, m.N)
    def Eq_tao_nrtl(m, i, i2, n):
        if i == i2:
            return m.tao_nrtl[i, i2, n] == 0
        else:
            return m.tao_nrtl[i, i2, n] == m.a_nrtl[i, i2] + (m.b_nrtl[i, i2]/m.Temp[n])

    g_nrtl_init = {}
    for n in m.N:
        for i in m.I:
            for i2 in m.I:
                if i == i2:
                    g_nrtl_init[(i, i2, n)] = 1
                else:
                    g_nrtl_init[(i, i2, n)] = pe.exp(-m.alfa_nrtl[i,
                                                                  i2]*pe.value(m.tao_nrtl[i, i2, n]))

    m.g_nrtl = pe.Var(m.I, m.I, m.N, within=pe.Reals, initialize=g_nrtl_init)

    @m.Constraint(m.I, m.I, m.N)
    def Eq_g_nrtl(m, i, i2, n):
        if i == i2:
            return m.g_nrtl[i, i2, n] == 1
        else:
            return m.g_nrtl[i, i2, n] == pe.exp(-m.alfa_nrtl[i, i2]*m.tao_nrtl[i, i2, n])

    gamma_init = {}
    for n in m.N:
        for comp in m.I:
            gamma_init[comp, n] = exp(sum(pe.value(m.x[comp1, n])*pe.value(m.tao_nrtl[comp1, comp, n])*pe.value(m.g_nrtl[comp1, comp, n]) for comp1 in m.I)/sum(pe.value(m.x[comp1, n])*pe.value(m.g_nrtl[comp1, comp, n]) for comp1 in m.I)+sum(pe.value(m.x[comp1, n])*pe.value(m.g_nrtl[comp, comp1, n])/sum(pe.value(m.x[comp2, n])*pe.value(
                m.g_nrtl[comp2, comp1, n]) for comp2 in m.I)*(pe.value(m.tao_nrtl[comp, comp1, n])-sum(pe.value(m.x[comp2, n])*pe.value(m.tao_nrtl[comp2, comp1, n])*pe.value(m.g_nrtl[comp2, comp1, n]) for comp2 in m.I)/sum(pe.value(m.x[comp3, n])*pe.value(m.g_nrtl[comp3, comp1, n]) for comp3 in m.I)) for comp1 in m.I))

    m.gamma = pe.Var(m.I, m.N, within=pe.Reals, initialize=gamma_init)

    @m.Constraint(m.I, m.N)
    def Eqgamma(m, comp, n):
        return m.gamma[comp, n] == pe.exp(sum(m.x[comp1, n]*m.tao_nrtl[comp1, comp, n]*m.g_nrtl[comp1, comp, n] for comp1 in m.I)/sum(m.x[comp1, n]*m.g_nrtl[comp1, comp, n] for comp1 in m.I)+sum(m.x[comp1, n]*m.g_nrtl[comp, comp1, n]/sum(m.x[comp2, n]*m.g_nrtl[comp2, comp1, n] for comp2 in m.I)*(m.tao_nrtl[comp, comp1, n]-sum(m.x[comp2, n]*m.tao_nrtl[comp2, comp1, n]*m.g_nrtl[comp2, comp1, n] for comp2 in m.I)/sum(m.x[comp3, n]*m.g_nrtl[comp3, comp1, n] for comp3 in m.I)) for comp1 in m.I))

    # ______________________________ Section 10 ______________________________
    # Phi calculation

    Omega_init = {'Water': 0.344861, 'Ethanol': 0.643558, 'Glycerol': 0.51269}
    m.Omega = pe.Param(m.I, initialize=Omega_init)    # Acentric factor [*]
    m.phi = pe.Param(m.I, m.N,  initialize=0.985)

    # ______________________________ Section 11 ______________________________
    # Entalphy calculation

    # Cp constants [kJ/mol*K]
    C1c_init = {'Water': 35.86, 'Ethanol': 6.638866920,
                'Glycerol': 12.02602320}
    m.C1c = pe.Param(m.I, initialize=C1c_init)
    C2c_init = {'Water': -0.0229,
                'Ethanol': 0.2249913420, 'Glycerol': 0.4250839170}
    m.C2c = pe.Param(m.I, initialize=C2c_init)
    C3c_init = {'Water': 6*10**-5, 'Ethanol': -
                1.0887479*10**-4, 'Glycerol': -2.8316157*10**-4}
    m.C3c = pe.Param(m.I, initialize=C3c_init)
    C4c_init = {'Water': -4*10**-8, 'Ethanol': 1.83461616 *
                10**-8, 'Glycerol': 7.58362604*10**-8}
    m.C4c = pe.Param(m.I, initialize=C4c_init)
    m.Tref = pe.Param(initialize=298.15)    # Reference temperature [K]
    m.Hscale = pe.Param(initialize=1000)    # Scaling factor for entalphy

    HVi_init = {}
    for n in m.N:
        for i in m.I:
            HVi_init[i, n] = ((m.C1c[i]*(pe.value(m.Temp[n])-m.Tref)) + ((m.C2c[i]/2)*((pe.value(m.Temp[n])**2)-(m.Tref**2))) + (
                (m.C3c[i]/3)*((pe.value(m.Temp[n])**3)-(m.Tref**3))) + ((m.C4c[i]/4)*((pe.value(m.Temp[n])**4)-(m.Tref**4))))/m.Hscale

    m.HVi = pe.Var(m.I, m.N, within=pe.Reals, initialize=HVi_init)

    @m.Constraint(m.I, m.N)
    def EqHVi(m, i, n):
        return m.HVi[i, n] == ((m.C1c[i]*(m.Temp[n]-m.Tref)) + ((m.C2c[i]/2)*((m.Temp[n]**2)-(m.Tref**2))) + ((m.C3c[i]/3)*((m.Temp[n]**3)-(m.Tref**3))) + ((m.C4c[i]/4)*((m.Temp[n]**4)-(m.Tref**4))))/m.Hscale

    HV_init = {}
    for n in m.N:
        HV_init[n] = sum(pe.value(m.HVi[i, n]) *
                         pe.value(m.y[i, n])/100 for i in m.I)

    m.HV = pe.Var(m.N, within=pe.Reals, initialize=HV_init)

    @m.Constraint(m.N)
    def EqHV(m, n):
        return m.HV[n] == sum(m.HVi[i, n]*m.y[i, n]/100 for i in m.I)

    # Vaporization entalphy constants [kJ/mol]
    C1v_init = {'Water': 51546.000,
                'Ethanol': 55789.00, 'Glycerol': 1.106700*10**5}
    m.C1v = pe.Param(m.I, initialize=C1v_init)
    C2v_init = {'Water': 0.2840200, 'Ethanol': 0.312450, 'Glycerol': 0.4831900}
    m.C2v = pe.Param(m.I, initialize=C2v_init)
    C3v_init = {'Water': -0.1584300, 'Ethanol': 0, 'Glycerol': 0}
    m.C3v = pe.Param(m.I, initialize=C3v_init)
    C4v_init = {'Water': 0.237500, 'Ethanol': 0, 'Glycerol': 0}
    m.C4v = pe.Param(m.I, initialize=C4v_init)
    Tb_init = {'Water': 100, 'Ethanol': 78.29, 'Glycerol': 287.85}
    # Component boiling temperature @1atm [C]
    m.Tb = pe.Param(m.I, initialize=Tb_init)

    def Tred_init(m, i):
        return m.Tb[i]/m.Tcrit[i]
    m.Tred = pe.Param(m.I, initialize=Tred_init)  # Reduced temperature [*]

    DHvap_init = {}
    for n in m.N:
        for i in m.I:
            DHvap_init[i, n] = (m.C1v[i]*((1-(pe.value(m.Temp[n])/pe.value(m.Tcrit[i]))))**(m.C2v[i]+m.C3v[i]*(
                pe.value(m.Temp[n])/pe.value(m.Tcrit[i]))+(m.C4v[i]*(pe.value(m.Temp[n])/pe.value(m.Tcrit[i])))))/m.Hscale

    # Vaporization entalphy [kJ/mol]
    m.DHvap = pe.Var(m.I, m.N, within=pe.Reals, initialize=DHvap_init)

    @m.Constraint(m.I, m.N)
    def EqdHvap(m, i, n):
        return m.DHvap[i, n] == (m.C1v[i]*((1-(m.Temp[n]/m.Tcrit[i])))**(m.C2v[i]+m.C3v[i]*(m.Temp[n]/m.Tcrit[i])+(m.C4v[i]*(m.Temp[n]/m.Tcrit[i]))))/m.Hscale

    HLi_init = {}
    for n in m.N:
        for i in m.I:
            HLi_init[i, n] = pe.value(m.HVi[i, n])-pe.value(m.DHvap[i, n])

    # Liquid phase entalphy [kJ/mol]
    m.HLi = pe.Var(m.I, m.N, within=pe.Reals, initialize=HLi_init)

    @m.Constraint(m.I, m.N)
    def EqHLi(m, i, n):
        return m.HLi[i, n] == m.HVi[i, n]-m.DHvap[i, n]

    HL_init = {}
    for n in m.N:
        HL_init[n] = sum(pe.value(m.HLi[i, n]) *
                         pe.value(m.x[i, n])/100 for i in m.I)

    m.HL = pe.Var(m.N, within=pe.Reals, initialize=HL_init)

    @m.Constraint(m.N)
    def EqHL(m, n):
        return m.HL[n] == sum(m.HLi[i, n]*m.x[i, n]/100 for i in m.I)

    # ______________________________ Section 12 ______________________________
    # Entalphy in feed calculation

    # Azeotropic feed
    def HV_Az_init(m, i):
        return ((m.C1c[i]*(m.Tali2-m.Tref)) + ((m.C2c[i]/2)*((m.Tali2**2)-(m.Tref**2))) + ((m.C3c[i]/3)*((m.Tali2**3)-(m.Tref**3))) + ((m.C4c[i]/4)*((m.Tali2**4)-(m.Tref**4))))/m.Hscale

    # Vapor enthalphy in feed [kJ/mol]
    m.HV_Az = pe.Param(m.I, initialize=HV_Az_init)

    def Tred_Az_init(m, i):
        return m.Tali2/m.Tcrit[i]

    # Reduced temperature in feed [*]
    m.Tred_Az = pe.Param(m.I, initialize=Tred_Az_init)

    def DHVap_Az_init(m, i):
        return (m.C1v[i]*((1-m.Tred_Az[i])**(m.C2v[i] + (m.C3v[i]*m.Tred_Az[i]) + (m.C4v[i]*(m.Tred_Az[i]**2)))))/m.Hscale

    # Vaporization enthalphy in feed [kJ/mol]
    m.DHVap_Az = pe.Param(m.I, initialize=DHVap_Az_init)

    def HL_Az_init(m, i):
        return m.HV_Az[i]-m.DHVap_Az[i]

    # Liquid phase enthalphy in feed [kJ/mol]
    m.HL_Az = pe.Param(m.I, initialize=HL_Az_init)

    HFAz_init = {}
    for n in m.N:
        if n != NT and n != 1:
            HFAz_init[n] = sum((m.zAz[i]/100)*(m.HL_Az[i]) for i in m.I)
        else:
            HFAz_init[n] = 0

    # Azeotropic feed enthalpy
    m.HFAz = pe.Param(m.N, initialize=HFAz_init)

    # Glycerol feed
    def HV_g_init(m, i):
        return ((m.C1c[i]*(m.Tali1-m.Tref)) + ((m.C2c[i]/2)*((m.Tali1**2)-(m.Tref**2))) + ((m.C3c[i]/3)*((m.Tali1**3)-(m.Tref**3))) + ((m.C4c[i]/4)*((m.Tali1**4)-(m.Tref**4))))/m.Hscale

    # Vapor enthalphy in feed [kJ/mol]
    m.HV_g = pe.Param(m.I, initialize=HV_g_init)

    def Tred_g_init(m, i):
        return m.Tali1/m.Tcrit[i]

    # Reduced temperature in feed [*]
    m.Tred_g = pe.Param(m.I, initialize=Tred_g_init)

    def DHVap_g_init(m, i):
        return (m.C1v[i]*((1-m.Tred_g[i])**(m.C2v[i] + (m.C3v[i]*m.Tred_g[i]) + (m.C4v[i]*(m.Tred_g[i]**2)))))/m.Hscale

    # Vaporization enthalphy in feed [kJ/mol]
    m.DHVap_g = pe.Param(m.I, initialize=DHVap_g_init)

    def HL_g_init(m, i):
        return m.HV_g[i]-m.DHVap_g[i]

    # Liquid phase enthalphy in feed [kJ/mol]
    m.HL_g = pe.Param(m.I, initialize=HL_g_init)

    HFG_init = {}
    for n in m.N:
        if n != NT and n != 1:
            HFG_init[n] = sum((m.zg[i]/100)*(m.HL_g[i]) for i in m.I)
        else:
            HFG_init[n] = 0

    m.HFG = pe.Param(m.N, initialize=HFG_init)  # Azeotropic feed enthalpy

    # ______________________________ Section 13 ______________________________
    # Parameter, constraint and binary variable definition

    yfG_init = {n: 0 for n in m.N}
    yfG_init[2] = 1

    yfAz_init = {n: 0 for n in m.N}
    yfAz_init[13] = 1

    yb_init = {n: 0 for n in m.N}
    yb_init[16] = 1

    yr_init = {n: 0 for n in m.N}
    yr_init[2] = 1

    par_init = {n: 0 for n in m.N}
    for j in range(2, 17):
        par_init[j] = 1

    # 1 if in stage n there is reflux, 0 otherwise (Binary?)
    m.yr = pe.Var(m.N, within=pe.Binary, initialize=yr_init)
    # 1 if in stage n there is boilup, 0 otherwise
    m.yb = pe.Var(m.N, within=pe.Binary, initialize=yb_init)
    # 1 if stage n exists, 0 otherwise
    m.par = pe.Var(m.N, within=pe.Binary, initialize=par_init)
    # 1 if in stage n there is Glycerol feed, 0 otherwise
    m.yfG = pe.Var(m.N, within=pe.Binary, initialize=yfG_init)
    # 1 if in stage n there is Azeotrope feed, 0 otherwise
    m.yfAz = pe.Var(m.N, within=pe.Binary, initialize=yfAz_init)

    #@m.Constraint(m.N)
    #def eqpar(m, n):
    #    if n == 1:
    #        return m.par[n] == 1
    #    elif n == NT:
    #        return m.par[n] == 1
    #    else:
    #        return m.par[n] == 1 - (1-sum(m.yr[j] for j in range(2, n+1)) - (sum(m.yb[j] for j in range(2, n))))

    # ______________________________ Section 14 ______________________________
    # Logic constraints

    m.cmej = pe.Param(initialize=1)

    @m.Constraint(m.N)
    def logic1(m, n):    # The boilup is below the reflux stage
        if n != NT and n != 1:
            return m.cmej*sum(m.yr[j] for j in range(2, n+1)) >= m.yb[n]*m.cmej
        else:
            return pe.Constraint.Skip

    @m.Constraint()
    def logic2(m):    # There is only one reflux stage
        return m.cmej*sum(m.yr[j] for j in range(2, NT)) == 1*m.cmej

    @m.Constraint()
    def logic3(m):    # There is only one boil-up stage
        return m.cmej*sum(m.yb[j] for j in range(2, NT)) == 1*m.cmej

    @m.Constraint()
    def logic4(m):    # There is only one feed of Glycerol
        return m.cmej*sum(m.yfG[j] for j in range(2, NT)) == 1*m.cmej

    @m.Constraint()
    def logic5(m):    # There is only one feed of Azeotrope
        return m.cmej*sum(m.yfAz[j] for j in range(2, NT)) == 1*m.cmej

    @m.Constraint(m.N)
    def logic6(m, n):    # Glycerol feed is below the reflux
        if n != NT and n != 1:
            return m.cmej*sum(m.yr[j] for j in range(2, n+1)) >= m.yfG[n]*m.cmej
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.N)
    def logic7(m, n):    # Azetropic feed is below the reflux
        if n != NT and n != 1:
            return m.cmej*sum(m.yr[j] for j in range(2, n+1)) >= m.yfAz[n]*m.cmej
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.N)
    def logic8(m, n):    # Boil-up stage is below Glycerol feed stage
        if n != NT and n != 1:
            return m.cmej*sum(m.yfG[j] for j in range(2, n+1)) >= m.yb[n]*m.cmej
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.N)
    def logic9(m, n):    # Boil-up stage is below Azeotropic feed stage
        if n != NT and n != 1:
            return m.cmej*sum(m.yfAz[j] for j in range(2, n+1)) >= m.yb[n]*m.cmej
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.N)
    def logic10(m, n):    # Glycerol feed is above the Azeotropic feed
        if n != NT and n != 1:
            return m.cmej*sum(m.yfG[j] for j in range(2, n+1)) >= m.yfAz[n]*m.cmej
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.N)
    def logic11(m, n):    # Glycerol feed stage is below reflux
        if n != NT and n != 1:
            return m.cmej*sum(m.yr[j] for j in range(2, n+1)) >= m.yfG[n]*m.cmej
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.N)
    def logic12(m, n):    # Azeotropic feed stage is below reflux
        if n != NT and n != 1:
            return m.cmej*sum(m.yr[j] for j in range(2, n+1)) >= m.yfAz[n]*m.cmej
        else:
            return pe.Constraint.Skip

    # ______________________________ Section 15 ______________________________
    # Condenser equations

    # Initial conditions for steady-state
    @m.Constraint()
    def BalMassC0(m):
        return 0 == m.V[2]-m.V[1]*(1+m.RR)

    @m.Constraint(m.I)
    def BalMassPartialC0(m, i):
        return 0 == m.V[2]*m.y[i, 2]-m.V[1]*m.x[i, 1]*(1+m.RR)

    @m.Constraint()
    def SumC0(m):
        return sum(m.y[i, 1]-m.x[i, 1] for i in m.I) == 0

    @m.Constraint(m.I)
    def EquilibriumC0(m, i):
        return m.y[i, 1]*m.P[1]*m.phi[i, 1] == m.Psat[i, 1]*m.gamma[i, 1]*m.x[i, 1]

    @m.Constraint()
    def BalEnergyC0(m):
        return 0 == m.V[2]*m.HV[2]-m.V[1]*(1+m.RR)*m.HL[1]-m.Qc

    # Fixed liquid flow
    @m.Constraint()
    def fixedL(m):
        return m.L[1] == 0

    # ______________________________ Section 16 ______________________________
    # Initial point column equations

    # Initial conditions for steady-state
    @m.Constraint(m.N)
    def BalMass0(m, n):
        if n != NT and n != 1:
            return 0 == m.yfG[n]*m.FG*1000 + m.yfAz[n]*m.FAz*1000 + m.RR*m.V[1]*m.yr[n] + m.BR*m.L[NT]*m.yb[n] + m.L[n-1] + m.V[n+1] - m.L[n] - m.V[n]
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.I, m.N)
    def BalMassPartial0(m, i, n):
        if n != NT and n != 1:
            return 0 == m.yfG[n]*m.FG*m.zg[i]*1000 + m.yfAz[n]*m.FAz*m.zAz[i]*1000 + m.RR*m.V[1]*m.yr[n]*m.x[i, 1] + m.BR*m.L[NT]*m.yb[n]*m.y[i, NT] + m.L[n-1]*m.x[i, n-1] + m.V[n+1]*m.y[i, n+1] - m.L[n]*m.x[i, n] - m.V[n]*m.y[i, n]
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.N)
    def Sum0(m, n):
        if n != NT and n != 1:
            return sum(m.x[i, n]-m.y[i, n] for i in m.I) == 0
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.N)
    def BalEnergy0(m, n):
        if n != NT and n != 1:
            return 0 == m.yfG[n]*m.FG*m.HFG[n]*1000 + m.yfAz[n]*m.FAz*m.HFAz[n]*1000 + m.RR*m.V[1]*m.yr[n]*m.HL[1] + m.BR*m.L[NT]*m.yb[n]*m.HV[NT] + m.L[n-1]*m.HL[n-1] + m.V[n+1]*m.HV[n+1] - m.L[n]*m.HL[n] - m.V[n]*m.HV[n]
        else:
            return pe.Constraint.Skip

    # Equilibrium relations
    @m.Constraint(m.I, m.N)
    def Equilibrium10(m, i, n):
        if n != NT and n != 1:
            return 0 == m.par[n]*((m.y[i, n]*m.P[n]*m.phi[i, 1])-(m.Psat[i, n]*m.gamma[i, n]*m.x[i, n]))
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.I, m.N)
    def Equilibrium20(m, i, n):
        if n != NT and n != 1 and i != 'Water':
            return sum(m.yr[j] for j in range(2, n+1))*(1-m.par[n])*(m.y[i, n]-m.y[i, n+1]) == 0
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.I, m.N)
    def Equilibrium30(m, i, n):
        if n != NT and n != 1 and i != 'Water':
            return (1-sum(m.yr[j] for j in range(2, n+1)))*(1-m.par[n])*(m.x[i, n]-m.x[i, n-1]) == 0
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.N)
    def Equilibrium40(m, n):
        if n != NT and n != 1:
            return 0 == (1-m.par[n])*(m.V[n]-m.V[n+1]-m.BR*m.L[NT]*m.yb[n])
        else:
            return pe.Constraint.Skip

    # ______________________________ Section 17 ______________________________
    # Reboiler equations

    # Initial conditions for steady-state
    @m.Constraint()
    def BalMassR0(m):
        return 0 == m.L[NT-1]-m.L[NT]*(1+m.BR)

    @m.Constraint(m.I)
    def BalMassPartialR0(m, i):
        return 0 == m.L[NT-1]*m.x[i, NT-1]-m.L[NT]*(m.x[i, NT]+m.BR*m.y[i, NT])

    @m.Constraint()
    def SumR0(m):
        return sum(m.y[i, NT]-m.x[i, NT] for i in m.I) == 0

    @m.Constraint(m.I)
    def EquilibriumR0(m, i):
        return m.y[i, NT]*m.P[NT]*m.phi[i, NT] == m.Psat[i, NT]*m.gamma[i, NT]*m.x[i, NT]

    @m.Constraint()
    def BalEnergyR0(m):
        return 0 == m.Qr+m.L[NT-1]*m.HL[NT-1]-m.L[NT]*m.HL[NT]-m.BR*m.L[NT]*m.HV[NT]

    # Fixed vapor flow in last stage
    @m.Constraint()
    def fixedV(m):
        return m.V[NT] == 0

    # ______________________________ Section 18 ______________________________
    # Hydraulic relations por all internal stages

    # Define vapor velocity
    far_init = {}
    for n in m.N:
        if n != NT and n != 1:
            far_init[n] = (0.981*exp(-0.411*(pe.value(m.Qvap[n])/3600)/pe.value(m.At)*(pe.value(m.rhoV[n])*sum(pe.value(m.MW[i])*pe.value(m.y[i, n])/100 for i in m.I)/1000)**(0.5)))
        else:
            far_init[n] = 0

    m.far = pe.Var(m.N, within=pe.NonNegativeReals, initialize=far_init)     # Aeration factor [*]

    @m.Constraint(m.N)
    def Eqfa(m, n):
        if n != NT and n != 1:
            return m.par[n]*(m.far[n]) == m.par[n]*(0.981*pe.exp(-0.411*(m.Qvap[n]/3600)/m.At*(m.rhoV[n]*sum(m.MW[i]*m.y[i, n]/100 for i in m.I)/1000)**(0.5)))
        else:
            return pe.Constraint.Skip

    hD_init = {}
    for n in m.N:
        if n != NT and n != 1:
            hD_init[n] = 0.6*((pe.value(m.Qliq[n])/3600/m.Lw)**(2/3))
        else:
            hD_init[n] = 0

    # Liquid height over divisor [m]
    m.hD = pe.Var(m.N, within=pe.NonNegativeReals, initialize=hD_init)
    @m.Constraint(m.N)
    def EqhD(m, n):
        if n != NT and n != 1:
            return (m.hD[n]) == 0.6*((m.Qliq[n]/3600/m.Lw)**(2/3))
        else:
            return pe.Constraint.Skip

    uhv_init = {}
    for n in m.N:
        if n != NT and n != 1:
            uhv_init[n] = ((pe.value(m.V[n])/(pe.value(m.rhoV[n]))))/(3600*pe.value(m.A0))
        else:
            uhv_init[n] = 0

    # Vapor velocity in hole [m/s]
    m.uhv = pe.Var(m.N, within=pe.NonNegativeReals, initialize=uhv_init)

    @m.Constraint(m.N)
    def Equhv(m, n):
        if n != NT and n != 1:
            return m.par[n]*(m.uhv[n]*m.A0) == m.par[n]*((m.V[n]/(m.rhoV[n])))/3600
        else:
            return pe.Constraint.Skip

    unv_init = {}
    for n in m.N:
        if n != NT and n != 1:
            unv_init[n] = ((pe.value(m.V[n])/(pe.value(m.rhoV[n]))))/(3600*pe.value(m.At))
        else:
            unv_init[n] = 0

    # Vapor velocity in plate  [m/s]
    m.unv = pe.Var(m.N, within=pe.NonNegativeReals, initialize=unv_init)

    @m.Constraint(m.N)
    def Equnv(m, n):
        if n != NT and n != 1:
            return m.par[n]*m.unv[n]*m.At == m.par[n]*((m.V[n]/(m.rhoV[n])))/3600
        else:
            return pe.Constraint.Skip

    ul_init = {}
    for n in m.N:
        if n != NT and n != 1:
            ul_init[n] = ((pe.value(m.L[n])/sum(pe.value(m.rho[i, n])*pe.value(m.x[i, n])/100 for i in m.I)))/(3600*pe.value(m.Ad))
        else:
            ul_init[n] = 0

    # Liquid velocity in weir  [m/s]
    m.ul = pe.Var(m.N, within=pe.NonNegativeReals, initialize=ul_init)
    @m.Constraint(m.N)
    def Equl(m, n):
        if n != NT and n != 1:
            return m.par[n]*m.ul[n]*m.Ad == m.par[n]*((m.L[n]/sum(m.rho[i, n]*m.x[i, n]/100 for i in m.I)))/3600
        else:
            return pe.Constraint.Skip

    # Liquid load
    m.consmach = pe.Param(initialize=1*10**-20)

    hcl_init = {}
    for n in m.N:
        if n != NT and n != 1:
            hcl_init[n] =  ((0.157*((pe.value(m.poro)+m.consmach)**(-0.791))/(1+1.04E-4*(((((pe.value(m.L[n])+m.consmach)/sum(pe.value(m.rho[i, n])*pe.value(m.x[i, n])/100 for i in m.I)))/3600/m.Lw)**(-0.59))*((pe.value(m.poro)+m.consmach)**(-1.791))))*(m.d_hole**0.833)*(996/(sum(pe.value(m.rho[i, n])*pe.value(m.x[i, n])/100 for i in m.I)*sum(m.MW[i]*pe.value(m.x[i, n])/100 for i in m.I)/1000))**(0.5*(1-0.91*m.d_hole/(pe.value(m.poro)+m.consmach))))
        else:
            hcl_init[n] = 0 

    # Liquid height in spray regime  [m]
    m.hcl = pe.Var(m.N, within=pe.NonNegativeReals, initialize=hcl_init)
 
    @m.Constraint(m.N)
    def Eqhcl(m, n):
        if n != NT and n != 1:
            return m.par[n]*m.hcl[n] == m.par[n]*((0.157*((m.poro+m.consmach)**(-0.791))/(1+1.04E-4*(((((m.L[n]+m.consmach)/sum(m.rho[i, n]*m.x[i, n]/100 for i in m.I)))/3600/m.Lw)**(-0.59))*((m.poro+m.consmach)**(-1.791))))*(m.d_hole**0.833)*(996/(sum(m.rho[i, n]*m.x[i, n]/100 for i in m.I)*sum(m.MW[i]*m.x[i, n]/100 for i in m.I)/1000))**(0.5*(1-0.91*m.d_hole/(m.poro+m.consmach))))
        else:
            return pe.Constraint.Skip

    m.Csbf = pe.Var(m.N, within=pe.NonNegativeReals, initialize=0.15)

    @m.Constraint(m.N)
    def EqCsbf(m, n):
        if n != NT and n != 1:
            return m.par[n]*(m.Csbf[n]) == m.par[n]*0.15
        else:
            return pe.Constraint.Skip

    # Pressure drop

    DPL_init = {}
    for n in m.N:
        if n != NT and n != 1:
            DPL_init[n] =  ((pe.value(m.far[n])*9.81*(sum(pe.value(m.rho[i, n])*pe.value(m.x[i, n])/100 for i in m.I)*sum(pe.value(m.MW[i])*pe.value(m.x[i, n])/100 for i in m.I)/1000)*(pe.value(m.hD[n])+pe.value(m.hw)))/101325)
        else:
            DPL_init[n] = 0

    # Pressure drop due to liquid presence [atm]
    m.DPL = pe.Var(m.N, within=pe.NonNegativeReals, initialize=DPL_init)
    
    @m.Constraint(m.N)
    def EqDPL(m, n):
        if n != NT and n != 1:
            return (m.DPL[n]) == ((m.far[n]*9.81*(sum(m.rho[i, n]*m.x[i, n]/100 for i in m.I)*sum(m.MW[i]*m.x[i, n]/100 for i in m.I)/1000)*(m.hD[n]+m.hw))/101325)
        else:
            return pe.Constraint.Skip

    # Pressure drop due to hole presence (dry) [atm]
    m.DPS = pe.Var(m.N, within=pe.NonNegativeReals, initialize=0.0001)

    @m.Constraint(m.N)
    def EqDPS(m, n):
        if n != NT and n != 1:
            return m.DPS[n] == 0.0001
        else:
            return pe.Constraint.Skip

    DPq_init = {}
    for n in m.N:
        if n != NT and n != 1:
            DPq_init[n] = (1.62*((sum(pe.value(m.rho[i,n])*pe.value(m.x[i,n])/100 for i in m.I)*sum(pe.value(m.MW[i])*pe.value(m.x[i,n])/100 for i in m.I)/1000))/(pe.sqrt(pe.value(m.Lw)*pe.value(m.hw)))*(pe.sqrt((pe.value(m.L[n])/sum(pe.value(m.rho[i,n])*pe.value(m.x[i,n])/100 for i in m.I)/3600))+pe.sqrt((pe.value(m.V[n])/(pe.value(m.rhoV[n]))/3600))))/101325
        
        else:
            DPq_init[n] = 0
    
    # Pressure drop at weir [atm]
    m.DPq = pe.Var(m.N, within=pe.NonNegativeReals, initialize=48.59)

    @m.Constraint(m.N)
    def EqDPq(m, n):
        if n != NT and n != 1:
            return m.DPq[n]*101325 == 1.62*((sum(m.rho[i,n]*m.x[i,n]/100 for i in m.I)*sum(m.MW[i]*m.x[i,n]/100 for i in m.I)/1000))/(pe.sqrt(m.Lw*m.hw))*(pe.sqrt((m.L[n]/sum(m.rho[i,n]*m.x[i,n]/100 for i in m.I)/3600))+pe.sqrt((m.V[n]/(m.rhoV[n])/3600)))
        else:
            return pe.Constraint.Skip

    DP_init = {}
    for n in m.N:
        if n != NT and n != 1:
            DP_init[n] = (pe.value(m.DPS[n])+pe.value(m.DPL[n]))
        else:
            DP_init[n] = 0

    # Total pressure drop [atm]
    m.DP = pe.Var(m.N, within=pe.NonNegativeReals, initialize=DP_init)

    # Define stage pressure
    @m.Constraint()
    def EqDPR(m):
        return m.DP[NT] == m.DP[NT-1]
    
    @m.Constraint(m.N)
    def EqDP(m, n):
        if n != NT and n != 1:
            return (m.DP[n]) == (m.DPS[n]+m.DPL[n])
        else:
            return pe.Constraint.Skip

    @m.Constraint(m.N)
    def EqP(m, n):
        if n != NT and n != 1:
            return m.P[n] == m.P[n-1]+m.par[n]*m.DP[n]
        else:
            return pe.Constraint.Skip

    @m.Constraint()
    def EqPC(m):
        return m.P[1] == m.Pop

    @m.Constraint()
    def EqPR(m):
        return m.P[NT] == m.P[NT-1]

    # Undesired effects in the column
    downF_init = {}
    for n in m.N:
        if n != NT and n != 1:
            downF_init[n] = ((pe.value(m.hD[n])+((pe.value(m.DP[n])*101325+pe.value(m.DPq[n])))/(9.81*(((sum(pe.value(m.rho[i, n])*pe.value(m.x[i, n])/100 for i in m.I)*sum(pe.value(m.MW[i])*pe.value(m.x[i, n])/100 for i in m.I)/1000))-(pe.value(m.rhoV[n])*sum(pe.value(m.MW[i])*pe.value(m.y[i, n])/100 for i in m.I)/1000))))-(pe.value(m.HS)))
        else:
            downF_init[n] = 0

    m.downF = pe.Var(m.N, within=pe.Reals, initialize=downF_init)

    @m.Constraint(m.N)
    def DownFlood(m, n):
        if n != NT and n != 1:
            return m.downF[n]*m.par[n] == ((m.hD[n]+((m.DP[n]*101325+m.DPq[n]))/(9.81*(((sum(m.rho[i, n]*m.x[i, n]/100 for i in m.I)*sum(m.MW[i]*m.x[i, n]/100 for i in m.I)/1000))-(m.rhoV[n]*sum(m.MW[i]*m.y[i, n]/100 for i in m.I)/1000))))-(m.HS))*m.par[n]
        else:
            return pe.Constraint.Skip

    # Entrainment flooding
    EntrainFloodVVar_init = {}
    for n in m.N:
        if n != NT and n != 1:
            EntrainFloodVVar_init[n] = (pe.value(m.unv[n])-(pe.value(m.Csbf[n])*(((((sum(pe.value(m.rho[i, n])*pe.value(m.x[i, n])/100 for i in m.I)*sum(pe.value(m.MW[i])*pe.value(m.x[i, n])/100 for i in m.I)/1000))-(pe.value(m.rhoV[n])*sum(pe.value(m.MW[i])*pe.value(m.y[i, n])/100 for i in m.I)/1000)))/(pe.value(m.rhoV[n])*sum(pe.value(m.MW[i])*pe.value(m.y[i, n])/100 for i in m.I)/1000))**0.5))
        else:
            EntrainFloodVVar_init[n] = 0

    m.EntrainFloodVVar = pe.Var(m.N, within=pe.Reals, initialize=EntrainFloodVVar_init)

    @m.Constraint(m.N)
    def EntrainFloodV(m, n):
        if n != NT and n != 1:
            return m.EntrainFloodVVar[n]*m.par[n] == m.par[n]*((m.unv[n])-(m.Csbf[n]*(((((sum(m.rho[i, n]*m.x[i, n]/100 for i in m.I)*sum(m.MW[i]*m.x[i, n]/100 for i in m.I)/1000))-(m.rhoV[n]*sum(m.MW[i]*m.y[i, n]/100 for i in m.I)/1000)))/(m.rhoV[n]*sum(m.MW[i]*m.y[i, n]/100 for i in m.I)/1000))**0.5))
        else:
            return pe.Constraint.Skip

    EntrainFloodLVar_init = {}
    for n in m.N:
        if n != NT and n != 1:
            EntrainFloodLVar_init[n] = (pe.value(m.ul[n])-((pe.value(m.sigma[n])*9.81*(((sum(pe.value(m.rho[i, n])*pe.value(m.x[i, n])/100 for i in m.I)*sum(pe.value(m.MW[i])*pe.value(m.x[i, n])/100 for i in m.I)/1000))-(pe.value(m.rhoV[n])*sum(pe.value(m.MW[i])*pe.value(m.y[i, n])/100 for i in m.I)/1000))/((sum(pe.value(m.rho[i, n])*pe.value(m.x[i, n])/100 for i in m.I)*sum(pe.value(m.MW[i])*pe.value(m.x[i, n])/100 for i in m.I)/1000)**2))**(1/4)))
        else:
            EntrainFloodLVar_init[n] = 0

    m.EntrainFloodLVar = pe.Var(m.N, within=pe.Reals, initialize=EntrainFloodLVar_init)
    @m.Constraint(m.N)
    def EntrainFloodL(m, n):
        if n != NT and n != 1:
            return m.EntrainFloodLVar[n]*m.par[n] == m.par[n]*((m.ul[n])-((m.sigma[n]*9.81*(((sum(m.rho[i, n]*m.x[i, n]/100 for i in m.I)*sum(m.MW[i]*m.x[i, n]/100 for i in m.I)/1000))-(m.rhoV[n]*sum(m.MW[i]*m.y[i, n]/100 for i in m.I)/1000))/((sum(m.rho[i, n]*m.x[i, n]/100 for i in m.I)*sum(m.MW[i]*m.x[i, n]/100 for i in m.I)/1000)**2))**(1/4)))
        else:
            return pe.Constraint.Skip

    # Weeping
    @m.Constraint(m.N)
    def Weep(m, n):
        if n != NT and n != 1:
            return 0 >= (((0.68+0.12)/(((m.rhoV[n]*sum(m.MW[i]*m.y[i, n]/100 for i in m.I)/1000)/((sum(m.rho[i, n]*m.x[i, n]/100 for i in m.I)*sum(m.MW[i]*m.x[i, n]/100 for i in m.I)/1000)*9.81*m.far[n]*(m.hw+m.hD[n])))**0.5))-(m.uhv[n]))*m.par[n]
        else:
            return pe.Constraint.Skip

    # Column construction

    # Column size
    
    m.Htotal = pe.Var(within=pe.NonNegativeReals, bounds=(6, 30), initialize=(((1+m.Sfactor)*sum(m.HS*pe.value(m.par[j]) for j in range(2, NT)))))   # Total column height [m]
    @m.Constraint()
    def Size(m):
        return m.Htotal == ((1+m.Sfactor)*sum(m.HS*m.par[j] for j in range(2, NT)))

    # Diameter-height relation
    @m.Constraint()
    def DtoLratio(m):
        return m.Htotal/m.D <= 20

    # ______________________________ Section 19 ______________________________
    # Objective function and infrastructure costs

    # Cost of energy for reboiler [$/kJ*h]
    m.CR_u = pe.Param(initialize=1*10**-5)
    # Cost of energy for condenser [$/kJ*h]
    m.CC_u = pe.Param(initialize=24.5/3600/6000)
    m.CEt = pe.Param(initialize=0.03)  # Ethanol cost [$/mol]

    # CEPCI Indexes
    m.I1 = pe.Param(initialize=397)
    m.I2 = pe.Param(initialize=575.4)

    # Aspen parameters
    m.sigma_d = pe.Param(initialize=50)  # Average superficial tension [dyn/cm]
    # Reference superficial tension [dyn/cm]
    m.sigma0 = pe.Param(initialize=20)
    # Condenser transfer coefficent [kJ/K*hr*m**2]
    m.Ucond = pe.Param(initialize=0.95*3600)
    m.lmtd_cond = pe.Param(initialize=35)  # Condenser average LMTD
    # Reboiler transfer coefficent [kJ/K*hr*m**2]
    m.Ureb = pe.Param(initialize=0.65*3600)
    m.lmtd_reb = pe.Param(initialize=35)  # Reboiler average LMTD

    # Luyben parameters
    m.c0_fair = pe.Param(initialize=439)
    m.c1_fair = pe.Param(initialize=2.5)
    m.c2_fair = pe.Param(initialize=1.2)
    m.an_ova = pe.Param(initialize=0.8)

    # Finance factors
    m.ophour = pe.Param(initialize=6000)   # Operation hours [hr/yr]
    m.anfact = pe.Param(initialize=0.25)    # Anualizing factor

    # Exchanger dimensions
    m.Area_cond = pe.Var(within=pe.NonNegativeReals,initialize=1.25*pe.value(m.Qc)/(m.Ucond*m.lmtd_cond))  # Condenser area [m**2]
    m.Area_reb = pe.Var(within=pe.NonNegativeReals,initialize=1.25*pe.value(m.Qr)/(m.Ureb*m.lmtd_reb))  # Reboiler area [m**2]

    @m.Constraint()
    def def_a_cond(m):
        return m.Area_cond*m.Ucond*m.lmtd_cond == 1.25*m.Qc

    @m.Constraint()
    def def_a_reb(m):
        return m.Area_reb*m.Ureb*m.lmtd_reb == 1.25*m.Qr

    # Column design

    dliqprom_init = sum(pe.value(m.par[n])*pe.value(m.rho[i, n])*pe.value(m.x[i, n]) for i in m.I for n in m.N)/sum(pe.value(m.par[n]) for n in m.N)/100
    m.dliqprom = pe.Var(within=pe.NonNegativeReals, initialize=dliqprom_init) # Average liquid density [kmol/m**3]
    @m.Constraint()
    def Eqdliqprom(m):
        return m.dliqprom == sum(m.par[n]*m.rho[i, n]*m.x[i, n] for i in m.I for n in m.N)/sum(m.par[n] for n in m.N)/100

    dvapprom_init = sum(pe.value(m.par[n])*pe.value(m.rhoV[n]) for n in m.N)/sum(pe.value(m.par[n]) for n in m.N)
    m.dvapprom = pe.Var(within=pe.NonNegativeReals, initialize=dvapprom_init)     # Average vapor density [kmol/m**3]
    @m.Constraint()
    def Eqdvapprom(m):
        return m.dvapprom == sum(m.par[n]*m.rhoV[n] for n in m.N)/sum(m.par[n] for n in m.N)

    mwpromliq_init = (sum(pe.value(m.par[n])*pe.value(m.MW[i])*pe.value(m.x[i, n])/100 for i in m.I for n in m.N)/1000)/sum(pe.value(m.par[n]) for n in m.N)
    m.mwpromliq = pe.Var(within=pe.NonNegativeReals, initialize=mwpromliq_init) # Average liquid molecular weight [kg/kmol]
    @m.Constraint()
    def Eqmwpromliq(m):
        return m.mwpromliq == (sum(m.par[n]*m.MW[i]*m.x[i, n]/100 for i in m.I for n in m.N)/1000)/sum(m.par[n] for n in m.N)

    mwpromvap_init =(sum(pe.value(m.par[n])*pe.value(m.MW[i])*pe.value(m.y[i, n])/100 for i in m.I for n in m.N)/1000)/sum(pe.value(m.par[n]) for n in m.N)
    m.mwpromvap = pe.Var(within=pe.NonNegativeReals, initialize=mwpromvap_init) # Average vapor molecular weight [kg/kmol]
    @m.Constraint()
    def Eqmwpromvap(m):
        return m.mwpromvap == (sum(m.par[n]*m.MW[i]*m.y[i, n]/100 for i in m.I for n in m.N)/1000)/sum(m.par[n] for n in m.N)

    V_prom_init = sum(pe.value(m.par[n])*pe.value(m.V[n]) for n in m.N)/sum(pe.value(m.par[n]) for n in m.N)
    m.V_prom = pe.Var(within=pe.NonNegativeReals, initialize=V_prom_init) # Average vapor flow [kmol/hr]
    @m.Constraint()
    def EqV_prom(m):
        return m.V_prom == sum(m.par[n]*m.V[n] for n in m.N)/sum(m.par[n] for n in m.N)

    m.FP_prom = pe.Var(within=pe.NonNegativeReals, initialize=0.07)  # Average flow factor
    @m.Constraint()
    def EqFP_prom(m):
        return m.FP_prom == 0.07

    cap_par_init = pe.value(m.c0_fair)/(1+m.c1_fair*(m.FP_prom**m.c2_fair))
    m.cap_par = pe.Var(within=pe.NonNegativeReals,initialize=cap_par_init)   # Capacity parameter
    @m.Constraint()
    def Eqcap_par(m):
        return m.cap_par*(1+m.c1_fair*(m.FP_prom**m.c2_fair)) == m.c0_fair

    m.f_flood = pe.Var(within=pe.NonNegativeReals, initialize=10)   # Flooding parameter
    m.scaleflood = pe.Param(initialize=1000)    # Scaling factor for flooding

    @m.Constraint()
    def Eqf_flood(m):
        return pe.sqrt(m.f_flood*m.scaleflood) == pe.sqrt(m.cap_par*m.an_ova*(m.sigma0/m.sigma_d)**0.2)*(m.dliqprom*m.mwpromliq-m.dvapprom*m.mwpromvap)

    @m.Constraint()
    def EqA_col(m):
        return pe.sqrt(m.A_col*0.6*m.f_flood*m.scaleflood)*m.dvapprom*m.mwpromvap == pe.sqrt(m.mwpromvap*m.V_prom)

    # Construction costs
    # Cost of infrastructure [1E5 $]
    m.tot_inf = pe.Var(within=pe.NonNegativeReals)
    m.scale_cost = pe.Param(initialize=100000)    # Scaling factor for costs

    # Cost of condenser [1E5 $]
    cost_cond_init = (m.I2/m.I1)*10**(3.7803+(0.8569*pe.log10(pe.value(m.Area_cond)))+0.0349*pe.sqrt(pe.log10(pe.value(m.Area_cond))))
    m.cost_cond = pe.Var(within=pe.NonNegativeReals, initialize=240913)
    @m.Constraint()
    def Eqcost_cond(m):
        return m.cost_cond == (m.I2/m.I1)*10**(3.7803+(0.8569*pe.log10(m.Area_cond))+0.0349*pe.sqrt(pe.log10(m.Area_cond)))
    #print(pe.value(m.cost_cond))

    # Cost of reboiler [1E5 $]
    cost_reb_init = (pe.value(m.I2)/pe.value(m.I1))*10**(4.4646-(0.5277*pe.log10(pe.value(m.Area_reb)))+0.3955*pe.sqrt(pe.log10(pe.value(m.Area_reb))))
    m.cost_reb = pe.Var(within=pe.NonNegativeReals, initialize=132826)
    @m.Constraint()
    def Eqcost_reb(m):
        return m.cost_reb == (m.I2/m.I1)*10**(4.4646-(0.5277*pe.log10(m.Area_reb))+0.3955*pe.sqrt(pe.log10(m.Area_reb)))

    #print(pe.value(m.cost_reb))

    cost_col_init = (m.I2/m.I1)*10**(3.4974+(0.4485*pe.log10(pe.value(m.Htotal)*pe.value(m.A_col)))+0.1074*pe.sqrt(pe.log10(pe.value(m.Htotal)*pe.value(m.A_col))))
    m.cost_col = pe.Var(within=pe.NonNegativeReals,initialize=12972)    # Cost of column [1E5 $]
    @m.Constraint()
    def Eqcost_col(m):
        return m.cost_col == (m.I2/m.I1)*10**(3.4974+(0.4485*pe.log10(m.Htotal*m.A_col))+0.1074*pe.sqrt(pe.log10(m.Htotal*m.A_col)))

    cost_sta_init = (m.I2/m.I1)*(sum(pe.value(m.par[n]) for n in m.N)-2)*10**(2.9949+(0.4465*pe.log10(pe.value(m.A_col)))+0.3961*pe.sqrt(pe.log10(pe.value(m.A_col))))
    m.cost_sta = pe.Var(within=pe.NonNegativeReals,initialize=15950)    # Cost of stages [1E5 $]
    @m.Constraint()
    def Eqcost_sta(m):
        return m.cost_sta == (m.I2/m.I1)*(sum(m.par[n] for n in m.N)-2)*10**(2.9949+(0.4465*pe.log10(m.A_col))+0.3961*pe.sqrt(pe.log10(m.A_col)))

    #print(pe.value(m.cost_sta))
    print(pe.value(m.Area_cond))
    print(pe.value(m.Htotal))

    @m.Constraint()
    def Eqtot_inf(m):
        return m.tot_inf == (m.cost_sta+m.cost_col+m.cost_cond+m.cost_reb)*m.anfact

    # Objective function
    m.CostMP = pe.Var(within=pe.NonNegativeReals)    # Raw material cost
    m.CostQr = pe.Var(within=pe.NonNegativeReals)    # Reboiler operation cost
    m.CostQc = pe.Var(within=pe.NonNegativeReals)    # Condenser operation cost
    m.GanEth = pe.Var(within=pe.NonNegativeReals)    # Ethanol utility

    @m.Constraint()
    def EqCostMP(m):
        return m.CostMP == m.CEt/2*m.FAz*1000*m.ophour

    @m.Constraint()
    def EqCostQr(m):
        return m.CostQr == m.Qr*m.CR_u*m.ophour

    @m.Constraint()
    def EqCostQc(m):
        return m.CostQc == m.Qc*m.CC_u*m.ophour

    @m.Constraint()
    def EqGanEth(m):
        return m.GanEth == m.CEt*m.V[1]*m.ophour

    def obj_rule(m):
        return (m.CostMP + m.CostQr + m.CostQc - m.GanEth) + (m.tot_inf)

    m.obj = pe.Objective(rule=obj_rule, sense=pe.minimize)

    # ______________________________ Section 20 ______________________________
    # Variable bounds were established in declaration

    # ______________________________ Section 21 ______________________________
    # Model solution
    dir_path = os.path.dirname(os.path.abspath(__file__))
    gams_path = os.path.join(dir_path, "gamsfiles/")
    if not(os.path.exists(gams_path)):
        print('Directory for automatically generated files ' +
              gams_path + ' does not exist. We will create it')
        os.makedirs(gams_path)
    solvername = 'gams'
    opt = SolverFactory(solvername, solver='sbb')
    results = opt.solve(m, tee=True,
                        # Uncomment the following lines if you want to save GAMS models
                        # keepfiles=True,
                        # tmpdir=gams_path,
                        # symbolic_solver_labels=True,

                        add_options=[
                            'option reslim = 4000;'
                            'option optcr = 0.0;'
                            'option nlp = conopt4;'
                            # Uncomment the following lines to setup IIS computation of BARON through option file
                            #'GAMS_MODEL.optfile = 1;'
                            # '\n'
                            #'$onecho > baron.opt \n'
                            #'CompIIS 1 \n'
                            # '$offecho'
                        ])

    print('Objective:', round(pe.value(m.obj), 3))

    return m


if __name__ == "__main__":
    NT = 30
    m = minlp_extractive_column(NT)
