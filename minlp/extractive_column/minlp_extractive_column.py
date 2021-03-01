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

    # ______________________________ Section 7 ______________________________
    # Calculation of superficial tension using critic DIPPR equation

    # Constants for DIPPR equation
    C1sig_init = {'Water':0.17766, 'Ethanol':0.0241005, 'Glycerol':0.0645335}
    m.C1sig = pe.Param(m.I, initialize=C1sig_init)
    C2sig_init = {'Water':2.567, 'Ethanol':-7.75658*10**-5, 'Glycerol':-5.38024}
    m.C2sig = pe.Param(m.I, initialize=C2sig_init)
    C3sig_init = {'Water':-3.3377, 'Ethanol':-1.025*10**-7, 'Glycerol':-2.1558*10**-7}
    m.C3sig = pe.Param(m.I, initialize=C3sig_init)
    C4sig_init = {'Water':1.9699, 'Ethanol':0, 'Glycerol':0}
    m.C4sig = pe.Param(m.I, initialize=C4sig_init)

    m.sigma = pe.Var(m.N, within=pe.NonNegativeReals, bounds=(0.005, 0.03)) # Liquid-vapor superficial tension [N/m]
    @m.Constraint(m.N)
    def Eqsigma(m,n):
        return m.sigma[n] == sum((m.x[i,n]/100)*m.C1sig[i]*(1-(m.Temp[n]/m.Tcritm[n]))**(m.C2sig[i]+m.C3sig[i]*(m.Temp[n]/m.Tcritm[n])+m.C4sig[i]*((m.Temp[n]/m.Tcritm[n]))**2) for i in m.I)

    # ______________________________ Sections 8-9 ______________________________
    # Calculation of activity coefficient using NRTL model

    a_nrtl_init= {('Water','Water'):0, ('Water','Ethanol'):3.4578, ('Water','Glycerol'):-1.2515,
                   ('Ethanol','Water'):-0.8009, ('Ethanol','Ethanol'):0, ('Ethanol','Glycerol'):0,
                   ('Glycerol','Water'):-0.7318, ('Glycerol','Ethanol'):0, ('Glycerol','Glycerol'):0}
    m.a_nrtl = pe.Param(m.I, m.I, initialize=a_nrtl_init)
    b_nrtl_init= {('Water','Water'):0, ('Water','Ethanol'):-586.0809, ('Water','Glycerol'):272.6075,
                   ('Ethanol','Water'):246.18, ('Ethanol','Ethanol'):0, ('Ethanol','Glycerol'):442.713,
                   ('Glycerol','Water'):170.9167, ('Glycerol','Ethanol'):36.139, ('Glycerol','Glycerol'):0}
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

    m.tao_nrtl = pe.Var(m.I, m.I, m.N, within=pe.Reals, bounds=(-5, 5))
    @m.Constraint(m.I, m.I, m.N)
    def Eq_tao_nrtl(m,i,i2,n):
        if i == i2:
            return m.tao_nrtl[i,i2,n] == 0
        else:
            return m.tao_nrtl[i,i2,n] == m.a_nrtl[i,i2] + (m.b_nrtl[i,i2]/m.Temp[n])
    
    m.g_nrtl = pe.Var(m.I, m.I, m.N, within=pe.Reals, bounds=(0, 2))
    @m.Constraint(m.I, m.I, m.N)
    def Eq_g_nrtl(m,i,i2,n):
        if i == i2:
            return m.g_nrtl[i,i2,n] == 1
        else:
            return m.g_nrtl[i,i2,n] == pe.exp(-m.alfa_nrtl[i,i2]*m.tao_nrtl[i,i2,n])
    
    m.gamma = pe.Var(m.I, m.N, within=pe.Reals, bounds=(0, 50))
    @m.Constraint(m.I, m.N)
    def Eqgamma(m,comp,n):
        return m.gamma[comp,n] == pe.exp(sum(m.x[comp1,n]*m.tao_nrtl[comp1,comp,n]*m.g_nrtl[comp1,comp,n] for comp1 in m.I)/sum(m.x[comp1,n]*m.g_nrtl[comp1,comp,n] for comp1 in m.I)+sum(m.x[comp1,n]*m.g_nrtl[comp,comp1,n]/sum(m.x[comp2,n]*m.g_nrtl[comp2,comp1,n] for comp2 in m.I)*(m.tao_nrtl[comp,comp1,n]-sum(m.x[comp2,n]*m.tao_nrtl[comp2,comp1,n]*m.g_nrtl[comp2,comp1,n] for comp2 in  m.I)/sum(m.x[comp3,n]*m.g_nrtl[comp3,comp1,n] for comp3 in m.I)) for comp1 in m.I))

    # ______________________________ Section 10 ______________________________
    # Phi calculation

    Omega_init = {'Water':0.344861, 'Ethanol':0.643558, 'Glycerol':0.51269}
    m.Omega = pe.Param(m.I, initialize=Omega_init)    # Acentric factor [*]
    m.phi = pe.Var(m.I, m.N, within=pe.NonNegativeReals, bounds=(0, 2))
    @m.Constraint(m.I, m.N)
    def EqPhi(m,i,n):
        return m.phi[i,n] == 0.985

    # ______________________________ Section 11 ______________________________
    # Entalphy calculation

    # Cp constants [kJ/mol*K]
    C1c_init = {'Water':35.86, 'Ethanol':6.638866920, 'Glycerol':12.02602320}
    m.C1c = pe.Param(m.I, initialize=C1c_init)
    C2c_init = {'Water':-0.0229, 'Ethanol':0.2249913420, 'Glycerol':0.4250839170}
    m.C2c = pe.Param(m.I, initialize=C2c_init)
    C3c_init = {'Water':6*10**-5, 'Ethanol':-1.0887479*10**-4, 'Glycerol':-2.8316157*10**-4}
    m.C3c = pe.Param(m.I, initialize=C3c_init)
    C4c_init = {'Water':-4*10**-8, 'Ethanol':1.83461616*10**-8, 'Glycerol':7.58362604*10**-8}
    m.C4c = pe.Param(m.I, initialize=C4c_init)
    m.Tref = pe.Param(initialize=298.15)    # Reference temperature [K]
    m.Hscale = pe.Param(initialize=1000)    # Scaling factor for entalphy

    m.HVi = pe.Var(m.I, m.N, within=pe.Reals, bounds=(-400, 1000))
    m.HV = pe.Var(m.N, within=pe.Reals, bounds=(-400, 1000))
    @m.Constraint(m.I, m.N)
    def EqHVi(m,i,n):
        return m.HVi[i,n] == ( (m.C1c[i]*(m.Temp[n]-m.Tref)) + ((m.C2c[i]/2)*((m.Temp[n]**2)-(m.Tref**2)))+ ((m.C3c[i]/3)*((m.Temp[n]**3)-(m.Tref**3))) + ((m.C4c[i]/4)*((m.Temp[n]**4)-(m.Tref**4))))/m.Hscale

    @m.Constraint(m.N)
    def EqHV(m,n):
        return m.HV[n]==sum(m.HVi[i,n]*m.y[i,n]/100 for i in m.I)

    # Vaporization entalphy constants [kJ/mol]
    C1v_init = {'Water':51546.000, 'Ethanol':55789.00, 'Glycerol':1.106700*10**5}
    m.C1v = pe.Param(m.I, initialize=C1v_init)
    C2v_init = {'Water':0.2840200, 'Ethanol':0.312450, 'Glycerol':0.4831900}
    m.C2v = pe.Param(m.I, initialize=C2v_init)
    C3v_init = {'Water':-0.1584300, 'Ethanol':0, 'Glycerol':0}
    m.C3v = pe.Param(m.I, initialize=C3v_init)
    C4v_init = {'Water':0.237500, 'Ethanol':0, 'Glycerol':0}
    m.C4v = pe.Param(m.I, initialize=C4v_init)
    Tb_init = {'Water':100, 'Ethanol':78.29, 'Glycerol':287.85}
    m.Tb = pe.Param(m.I, initialize=Tb_init)  # Component boiling temperature @1atm [C]

    def Tred_init(m,i):
        return m.Tb[i]/m.Tcrit[i]
    m.Tred = pe.Param(m.I, initialize=Tred_init)  # Reduced temperature [*]

    m.DHvap = pe.Var(m.I, m.N, within=pe.Reals) # Vaporization entalphy [kJ/mol]
    @m.Constraint(m.I, m.N)
    def EqdHvap(m,i,n):
        return m.DHvap[i,n] == (m.C1v[i]*((1-(m.Temp[n]/m.Tcrit[i])))**(m.C2v[i]+m.C3v[i]*(m.Temp[n]/m.Tcrit[i])+(m.C4v[i]*(m.Temp[n]/m.Tcrit[i]))))/m.Hscale

    # Liquid phase entalphy [kJ/mol]
    m.HLi = pe.Var(m.I, m.N, within=pe.Reals, bounds=(-1000, 1000))
    m.HL = pe.Var(m.N, within=pe.Reals, bounds=(-1000, 1000))
    @m.Constraint(m.I, m.N)
    def EqHLi(m,i,n):
        return m.HLi[i,n] == m.HVi[i,n]-m.DHvap[i,n]

    @m.Constraint(m.N)
    def EqHL(m,n):
        return m.HL[n] == sum(m.HLi[i,n]*m.x[i,n]/100 for i in m.I)

    # ______________________________ Section 12 ______________________________
    # Entalphy in feed calculation

    # Azeotropic feed
    def HV_Az_init(m,i):
        return ( (m.C1c[i]*(m.Tali2-m.Tref)) + ((m.C2c[i]/2)*((m.Tali2**2)-(m.Tref**2))) + ((m.C3c[i]/3)*((m.Tali2**3)-(m.Tref**3))) + ((m.C4c[i]/4)*((m.Tali2**4)-(m.Tref**4))))/m.Hscale

    m.HV_Az = pe.Param(m.I, initialize=HV_Az_init)  # Vapor enthalphy in feed [kJ/mol]

    def Tred_Az_init(m,i):
        return m.Tali2/m.Tcrit[i]

    m.Tred_Az = pe.Param(m.I, initialize=Tred_Az_init)  # Reduced temperature in feed [*]

    def DHVap_Az_init(m,i):
        return ( m.C1v[i]*( (1-m.Tred_Az[i])**( m.C2v[i] + (m.C3v[i]*m.Tred_Az[i]) + (m.C4v[i]*(m.Tred_Az[i]**2)) ) ) )/m.Hscale

    m.DHVap_Az = pe.Param(m.I, initialize=DHVap_Az_init)  # Vaporization enthalphy in feed [kJ/mol]

    def HL_Az_init(m,i):
        return m.HV_Az[i]-m.DHVap_Az[i]

    m.HL_Az = pe.Param(m.I, initialize=HL_Az_init)  # Liquid phase enthalphy in feed [kJ/mol]

    m.HFAz = pe.Var(m.N, within=pe.Reals) # Azeotropic feed enthalpy

    @m.Constraint(m.N)
    def EqHFAz(m,n):
        if n != NT and n != 1:
            return m.HFAz[n] == sum((m.zAz[i]/100)*(m.HL_Az[i]) for i in m.I)
        else:
            return pe.Constraint.Skip

    # Glycerol feed
    def HV_g_init(m,i):
        return ( (m.C1c[i]*(m.Tali1-m.Tref)) + ((m.C2c[i]/2)*((m.Tali1**2)-(m.Tref**2))) + ((m.C3c[i]/3)*((m.Tali1**3)-(m.Tref**3))) + ((m.C4c[i]/4)*((m.Tali1**4)-(m.Tref**4))))/m.Hscale
        
    m.HV_g = pe.Param(m.I, initialize=HV_g_init)  # Vapor enthalphy in feed [kJ/mol]

    def Tred_g_init(m,i):
        return m.Tali1/m.Tcrit[i]

    m.Tred_g = pe.Param(m.I, initialize=Tred_g_init)  # Reduced temperature in feed [*]

    def DHVap_g_init(m,i):
        return ( m.C1v[i]*( (1-m.Tred_g[i])**( m.C2v[i] + (m.C3v[i]*m.Tred_g[i]) + (m.C4v[i]*(m.Tred_g[i]**2)) ) ) )/m.Hscale

    m.DHVap_g = pe.Param(m.I, initialize=DHVap_g_init)  # Vaporization enthalphy in feed [kJ/mol]

    def HL_g_init(m,i):
        return m.HV_g[i]-m.DHVap_g[i]

    m.HL_g = pe.Param(m.I, initialize=HL_g_init)  # Liquid phase enthalphy in feed [kJ/mol]

    m.HFG = pe.Var(m.N, within=pe.Reals) # Azeotropic feed enthalpy

    @m.Constraint(m.N)
    def EqHFG(m,n):
        if n != NT and n != 1:
            return m.HFG[n] == sum((m.zg[i]/100)*(m.HL_g[i]) for i in m.I)
        else:
            return pe.Constraint.Skip

    # ______________________________ Section 13 ______________________________
    # Parameter, constraint and binary variable definition

    def yr_init(m,n):
        if n == 2:
            return 1
        else:
            return 0

    m.yr = pe.Param(m.N, initialize=yr_init)    # 1 if in stage n there is reflux, 0 otherwise
    m.yb = pe.Var(m.N, within=pe.Binary)    # 1 if in stage n there is boilup, 0 otherwise
    m.par = pe.Var(m.N, within=pe.Binary)    # 1 if stage n exists, 0 otherwise   (real?)
    m.yfG = pe.Var(m.N, within=pe.Binary)    # 1 if in stage n there is Glycerol feed, 0 otherwise
    m.yfAz = pe.Var(m.N, within=pe.Binary)    # 1 if in stage n there is Azeotrope feed, 0 otherwise



    return m




if __name__ == "__main__":
    NT = 30
    m = minlp_extractive_column(NT)
