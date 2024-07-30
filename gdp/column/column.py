"""
column.py
Fixed design solution of distillation column GDP model

This file defines an optimization model for the design and operation of a distillation column for benzene-toluene separation.
The objective is to minimize the operation cost (heat duties in condenser and reboiler) and fixed cost (number of trays in the column).
The constraints are the MESH equations (material balance, equilibrium, summation, and enthalpy balance) for each tray together with logical constraints that encode the existence of trays and position of reflux and boilup flows.
The continuous variables of this model are the flowrates of each component in liquid and vapor phase and temperatures at each tray, the reflux and boilup ratio, and the condenser and reboiler heat duties.
The logical variables are the exitence or non-existence of the trays, and the position of the reflux and boilup flows.
The complete model defines a Generalized Disjunctive Programming (GDP) problem.

After the model is defined, the boolean variables are reformulated into integer variables known as external variables.
The external variables are defined as user input and are used to activate or deactivate the trays and position the reflux and boilup flows.
Using the GDP problem, we activate and deactivate selectively corresponding constraints.
The model is initialized using data from an Excel sheet 'init.xlsx' and other model-related calculations.
The resulting problem is a Nonlinear Programming (NLP) problem, which is solved through GAMS solvers specified by the user.

References:
- Ghouse, Jaffer H., et al. "A comparative study between GDP and NLP formulations for conceptual design of distillation columns." Computer Aided Chemical Engineering. Vol. 44. Elsevier, 2018. 865-870.
- Bernal, David E., et al. "Process Superstructure Optimization through Discrete Steepest Descent Optimization: a GDP Analysis and Applications in Process Intensification." Computer Aided Chemical Engineering. Vol. 49. Elsevier, 2022. 1279-1284.
"""

from __future__ import division

import math
import os
import time

from pyomo.common.errors import InfeasibleConstraintException
from pyomo.contrib.fbbt.fbbt import fbbt
#from pyomo.contrib.gdpopt.data_class import MasterProblemResult
from pyomo.environ import (
    Block,
    BooleanVar,
    ConcreteModel,
    Constraint,
    NonNegativeReals,
    Objective,
    Param,
    RangeSet,
    Set,
    SolverFactory,
    TransformationFactory,
    Var,
    exactly,
    exp,
    land,
    log,
    lor,
    minimize,
    value,
)
from pyomo.gdp import Disjunct, Disjunction
from pyomo.opt import SolutionStatus, SolverResults
from pyomo.opt import TerminationCondition as tc
from pyomo.util.infeasible import log_infeasible_constraints
import pandas as pd


def initialize(m):
    """
    Initializes the values of the distillation model using provided data from an Excel sheet 'init.xlsx'
    and other model-related calculations.

    Args:
        m (pyomo.ConcreteModel): Model object representing the distillation column.

    Returns:
        None. The function modifies the model 'm' in-place.
    """
    m.reflux_frac.set_value(value(m.reflux_ratio / (1 + m.reflux_ratio)))
    m.boilup_frac.set_value(value(m.reboil_ratio / (1 + m.reboil_ratio)))

    # Import the Excel Sheet
    _excel_sheets = pd.read_excel('init.xlsx', sheet_name=None, engine='openpyxl')

    def set_value_if_not_fixed(var, val):
        """
        Set variable to the value if it is not fixed.
        Checks if the given variable is fixed. If not, sets its value to the provided value.

        Args:
            var (pyomo.Variable): A variable, typically part of a mathematical model.
            val (float): The value to be set for the variable if it's not fixed.

        Returns:
            None. The function modifies the variable 'var' in-place if it's not fixed.
        """
        if not var.fixed:
            var.set_value(val)

    # active_trays are the condenser, reboiler, feed, and those conditional trays whose indicator_var is 1
    active_trays = [
        t
        for t in m.trays
        if t not in m.conditional_trays
        or math.fabs(value(m.tray[t].indicator_var - 1)) <= 1e-3
    ]
    num_active_trays = len(active_trays)

    feed_tray = m.feed_tray

    # Process 'trays' sheet with temperature from Excel
    tray_indexed_data = _excel_sheets['trays']
    tray_indexed_data.sort_values(by=['tray'], inplace=True)
    tray_indexed_data.set_index('tray', inplace=True)

    # Process 'comps_and_trays' sheet with flows and compositions from Excel and set multi-index
    comp_and_tray_indexed_data = _excel_sheets['comps_and_trays']
    comp_and_tray_indexed_data.sort_values(by=['comp', 'tray'], inplace=True)
    comp_and_tray_indexed_data.set_index(['comp', 'tray'], inplace=True)
    comp_slices = {
        c: comp_and_tray_indexed_data.loc[c, :] for c in m.comps
    }  # Create dictionary of data slices per component

    num_data_trays = tray_indexed_data.index.size
    # If there are fewer active trays than data trays, average the data.
    if num_active_trays < num_data_trays:
        # Number of model trays is less than number of trays in data. Need to
        # do averaging
        new_indices = [1] + [
            1 + (num_data_trays - 1) / (num_active_trays - 1) * i
            for i in range(1, num_active_trays)
        ]
        # Combine data from adjacent trays based on linear combination.
        for tray in range(2, num_active_trays):
            indx = new_indices[tray - 1]
            lower = math.floor(indx)
            frac_above = indx - lower
            # Take linear combination of values
            tray_indexed_data.loc[tray] = (
                tray_indexed_data.loc[lower] * (1 - frac_above)
                + tray_indexed_data.loc[lower + 1] * frac_above
            )
            for c in m.comps:
                comp_slices[c].loc[tray] = (
                    comp_slices[c].loc[lower] * (1 - frac_above)
                    + comp_slices[c].loc[lower + 1] * frac_above
                )
        tray_indexed_data.loc[num_active_trays] = tray_indexed_data.loc[num_data_trays]
        # Trim data to match the number of active trays.
        tray_indexed_data = tray_indexed_data.head(num_active_trays)
        for c in m.comps:
            comp_slices[c].loc[num_active_trays] = comp_slices[c].loc[num_data_trays]
            comp_slices[c] = comp_slices[c].head(num_active_trays)
    # If there are more or equal active trays than data trays, interpolate the data.
    else:
        # Stretch the data out and do interpolation
        tray_indexed_data.index = pd.Index(
            [1]
            + [
                int(round(num_active_trays / num_data_trays * i))
                for i in range(2, num_data_trays + 1)
            ],
            name='tray',
        )
        tray_indexed_data = tray_indexed_data.reindex(
            [i for i in range(1, num_active_trays + 1)]
        ).interpolate()
        # Stretch component data, handle potential N/A values, and interpolate
        for c in m.comps:
            comp_slices[c].index = pd.Index(
                [1]
                + [
                    int(round(num_active_trays / num_data_trays * i))
                    for i in range(2, num_data_trays + 1)
                ],
                name='tray',
            )
            # special handling necessary for V near top of column and L
            # near column bottom. Do not want to interpolate with one end
            # being potentially 0. (ie. V from total condenser). Instead,
            # use back fill and forward fill.
            comp_slices[c] = comp_slices[c].reindex(
                [i for i in range(1, num_active_trays + 1)]
            )
            tray_below_condenser = sorted(active_trays, reverse=True)[1]
            if pd.isna(comp_slices[c]['V'][tray_below_condenser]):
                # V of the tray below the condenser is N/A. Find a valid
                # value lower down to use.
                val = next(
                    comp_slices[c]['V'][t]
                    for t in reversed(list(m.trays))
                    if pd.notna(comp_slices[c]['V'][t]) and not t == m.condens_tray
                )
                comp_slices[c]['V'][tray_below_condenser] = val
            if pd.isna(comp_slices[c]['L'][m.reboil_tray + 1]):
                # L of the tray above the reboiler is N/A. Find a valid
                # value higher up to use.
                val = next(
                    comp_slices[c]['L'][t]
                    for t in m.trays
                    if pd.notna(comp_slices[c]['L'][t]) and not t == m.reboil_tray
                )
                comp_slices[c]['L'][m.reboil_tray + 1] = val
            comp_slices[c] = comp_slices[c].interpolate()

    # Reindex the tray data and fill any missing values using back-fill method.
    tray_indexed_data.index = pd.Index(sorted(active_trays), name='tray')
    tray_indexed_data = tray_indexed_data.reindex(sorted(m.trays), method='bfill')

    # Set the temperature values for trays, if not fixed.
    for t in m.trays:
        set_value_if_not_fixed(m.T[t], tray_indexed_data['T [K]'][t])

    # Reindex and fill component data. Use back-fill for 'L' and 'x', and forward-fill for 'V' and 'y'
    for c in m.comps:
        comp_slices[c].index = pd.Index(sorted(active_trays), name='tray')
        comp_slices[c] = comp_slices[c].reindex(sorted(m.trays))
        comp_slices[c][['L', 'x']] = comp_slices[c][['L', 'x']].bfill()
        comp_slices[c][['V', 'y']] = comp_slices[c][['V', 'y']].ffill()

    comp_and_tray_indexed_data = pd.concat(comp_slices)

    # Set component values for each tray if they are not fixed.
    for c, t in m.comps * m.trays:
        set_value_if_not_fixed(m.L[c, t], comp_and_tray_indexed_data['L'][c, t])
        set_value_if_not_fixed(m.V[c, t], comp_and_tray_indexed_data['V'][c, t])
        set_value_if_not_fixed(m.x[c, t], comp_and_tray_indexed_data['x'][c, t])
        set_value_if_not_fixed(m.y[c, t], comp_and_tray_indexed_data['y'][c, t])

    # Set enthalpy specifications for each component in the feed.
    for c in m.comps:
        m.H_L_spec_feed[c].set_value(value(m.feed_liq_enthalpy_expr[c]))
        m.H_V_spec_feed[c].set_value(value(m.feed_vap_enthalpy_expr[c]))

    # Compute and set several values for each component in each tray.
    for t in m.trays:
        for c in m.comps:
            k = m.pvap_const[c]
            x = m.Pvap_X[c, t]

            x.set_value(value(1 - m.T[t] / k['Tc']))

            m.Pvap[c, t].set_value(
                value(
                    exp(
                        (
                            k['A'] * x
                            + k['B'] * x**1.5
                            + k['C'] * x**3
                            + k['D'] * x**6
                        )
                        / (1 - x)
                    )
                    * k['Pc']
                )
            )

            m.Kc[c, t].set_value(value(m.gamma[c, t] * m.Pvap[c, t] / m.P))

            m.H_L[c, t].set_value(value(m.liq_enthalpy_expr[t, c]))
            m.H_V[c, t].set_value(value(m.vap_enthalpy_expr[t, c]))

    # Setting initial values for distillate (D) and bottoms (B) for benzene and toluene.
    # NOTE: This can be improved by bringing up the spreadsheed data and using boil-up and reflux ratio.
    m.D['benzene'].set_value(42.3152714)
    m.D['toluene'].set_value(5.4446286)
    m.B['benzene'].set_value(7.67928)
    m.B['toluene'].set_value(44.56072)
    m.L['benzene', m.reboil_tray].set_value(7.67928)
    m.L['toluene', m.reboil_tray].set_value(44.56072)

    # Calculating and setting the vapor values (V) for benzene and toluene at the reboil tray.
    m.V['benzene', m.reboil_tray].set_value(
        value(m.L['benzene', m.reboil_tray + 1] - m.L['benzene', m.reboil_tray])
    )
    m.V['toluene', m.reboil_tray].set_value(
        value(m.L['toluene', m.reboil_tray + 1] - m.L['toluene', m.reboil_tray])
    )

    # Calculating and setting the liquid values (L) for benzene and toluene at the condensate tray.
    m.L['benzene', m.condens_tray].set_value(
        value(m.V['benzene', m.condens_tray - 1] - m.D['benzene'])
    )
    m.L['toluene', m.condens_tray].set_value(
        value(m.V['toluene', m.condens_tray - 1] - m.D['toluene'])
    )

    # Calculating and setting total liquid (liq) and vapor (vap) values for each tray
    for t in m.trays:
        m.liq[t].set_value(value(sum(m.L[c, t] for c in m.comps)))
        m.vap[t].set_value(value(sum(m.V[c, t] for c in m.comps)))

    # Setting the bottom and distillate values.
    m.bot.set_value(
        52.24
    )  # NOTE: This could be initialized as the sum of component flows in the Excel spreadsheet
    m.dis.set_value(
        47.7599
    )  # NOTE: This could be initialized as the sum of component flows in the Excel spreadsheet

    # Calculating and setting mole fraction values (x and y) for components in the reboil and condensate trays.
    for c in m.comps:
        m.x[c, m.reboil_tray].set_value(
            value(m.L[c, m.reboil_tray] / m.liq[m.reboil_tray])
        )
        m.y[c, m.reboil_tray].set_value(
            value(m.V[c, m.reboil_tray] / m.vap[m.reboil_tray])
        )
        m.x[c, m.condens_tray].set_value(
            value(m.L[c, m.condens_tray] / m.liq[m.condens_tray])
        )
        m.y[c, m.condens_tray].set_value(
            value(m.x[c, m.condens_tray] * m.Kc[c, m.condens_tray])
        )

    # Setting heat values for the boiler and condenser.
    m.Qb.set_value(2.307873115)
    m.Qc.set_value(3.62641882)


def build_column(
    min_trays,
    max_trays,
    xD,
    xB,
    x_input,
    nlp_solver,
    provide_init=False,
    init={},
    boolean_ref=False,
):
    t_start = time.process_time()
    """
    Builds the column model  which is a Generalized Disjunctive Program (GDP)
    References: Ghouse, Jaffer H., et al. "A comparative study between GDP and NLP formulations for conceptual design of distillation columns." Computer Aided Chemical Engineering. Vol. 44. Elsevier, 2018. 865-870.
    
    Args:
        min_trays (int): Minimum number of trays in the column
        max_trays (int): Maximum number of trays in the column
        xD (float): Distillate purity
        xB (float): Bottoms purity
        x_input (list): List of the external variable values with the reflex position and the boilup position in the column
        nlp_solver (str): Name of the NLP solver to use
        provide_init (bool): Whether to provide initialization values
        init (dict): Dictionary of initialization values
        boolean_ref (bool): Whether to use boolean reformulation
    Returns:
        m (pyomo.ConcreteModel): Pyomo model
    """
    m = ConcreteModel('benzene-toluene column')
    m.comps = Set(initialize=['benzene', 'toluene'])  # Initialize component set
    min_T, max_T = 300, 400  # [K] Define temperature range
    m.T_ref = 298.15  # Reference temperature [K]
    max_flow = 500  # Define maximum flow rate [mol/s]
    # Define number of trays, condenser and reboiler trays
    m.max_trays = max_trays
    m.condens_tray = max_trays  # Condenser at the top of the column
    m.feed_tray = math.ceil((max_trays / 2))  # Feed in the middle of the column
    m.reboil_tray = 1  # Reboiler at the bottom of the column
    m.distillate_purity = xD  # Purity of the distillate.
    m.bottoms_purity = xB  # Purity of the bottoms.
    # Vapor pressure constants for benzene and toluene
    m.pvap_const = {
        'benzene': {
            'A': -6.98273,
            'B': 1.33213,
            'C': -2.62863,
            'D': -3.33399,
            'Tc': 562.2,
            'Pc': 48.9,
        },
        'toluene': {
            'A': -7.28607,
            'B': 1.38091,
            'C': -2.83433,
            'D': -2.79168,
            'Tc': 591.8,
            'Pc': 41.0,
        },
    }
    # Heat capacity constants for vapor phase and liquid phase
    m.vap_Cp_const = {
        'benzene': {
            'A': -3.392e1,
            'B': 4.739e-1,
            'C': -3.017e-4,
            'D': 7.130e-8,
            'E': 0,
        },
        'toluene': {
            'A': -2.435e1,
            'B': 5.125e-1,
            'C': -2.765e-4,
            'D': 4.911e-8,
            'E': 0,
        },
    }
    m.liq_Cp_const = {
        'benzene': {'A': 1.29e5, 'B': -1.7e2, 'C': 6.48e-1, 'D': 0, 'E': 0},
        'toluene': {'A': 1.40e5, 'B': -1.52e2, 'C': 6.95e-1, 'D': 0, 'E': 0},
    }
    # Heat of vaporization for each component [J/mol]
    m.dH_vap = {'benzene': 33.770e3, 'toluene': 38.262e3}

    # Define set of potential trays
    m.trays = RangeSet(max_trays, doc='Set of potential trays')
    # Define set of trays that can be turned on and off
    m.conditional_trays = Set(
        initialize=m.trays - [m.condens_tray, m.feed_tray, m.reboil_tray],
        doc="Trays that may be turned on and off.",
    )
    # Disjunct for tray existence
    m.tray = Disjunct(m.conditional_trays, doc='Disjunct for tray existence')
    # Disjunct for tray absence
    m.no_tray = Disjunct(m.conditional_trays, doc='Disjunct for tray absence')

    # Disjunction statement defining whether a tray exists or not
    @m.Disjunction(m.conditional_trays, doc='Tray exists or does not')
    def tray_no_tray(b, t):  # NOTE: This function is not accesed
        """Disjunction statement defining whether a tray exists or not"""
        return [b.tray[t], b.no_tray[t]]

    # Constraint for minimum number of trays, adding 1 for feed tray
    m.minimum_num_trays = Constraint(
        expr=sum(m.tray[t].indicator_var for t in m.conditional_trays)
        + 1  # for feed tray
        >= min_trays
    )

    # If user provides initialization values, use them
    if provide_init:
        # Feed temperature variable [K]
        m.T_feed = Var(
            doc='Feed temperature [K]',
            domain=NonNegativeReals,
            bounds=(min_T, max_T),
            initialize=init['T_feed'],
        )

        # Vapor fraction of feed variable
        m.feed_vap_frac = Var(
            doc='Vapor fraction of feed',
            initialize=init['feed_vap_frac'],
            bounds=(0, 1),
        )

        # Component feed flow variable [mol/s]
        m.feed = Var(
            m.comps, doc='Total component feed flow [mol/s]', initialize=init['feed']
        )

        # Liquid mole fraction variable
        m.x = Var(
            m.comps,
            m.trays,
            doc='Liquid mole fraction',
            bounds=(0, 1),
            domain=NonNegativeReals,
            initialize=init['x'],
        )

        # Vapor mole fraction variable
        m.y = Var(
            m.comps,
            m.trays,
            doc='Vapor mole fraction',
            bounds=(0, 1),
            domain=NonNegativeReals,
            initialize=init['y'],
        )

        # Component liquid flows from tray variable [mol/s]
        m.L = Var(
            m.comps,
            m.trays,
            doc='component liquid flows from tray [mol/s]',
            domain=NonNegativeReals,
            bounds=(0, max_flow),
            initialize=init['L'],
        )

        # Component vapor flows from tray variable [mol/s]
        m.V = Var(
            m.comps,
            m.trays,
            doc='component vapor flows from tray [mol/s]',
            domain=NonNegativeReals,
            bounds=(0, max_flow),
            initialize=init['V'],
        )

        # Liquid flows from tray variable [mol/s]
        m.liq = Var(
            m.trays,
            domain=NonNegativeReals,
            doc='liquid flows from tray [mol/s]',
            initialize=init['liq'],
            bounds=(0, max_flow),
        )

        # Vapor flows from tray variable [mol/s]
        m.vap = Var(
            m.trays,
            domain=NonNegativeReals,
            doc='vapor flows from tray [mol/s]',
            initialize=init['vap'],
            bounds=(0, max_flow),
        )

        # Bottoms component flows variable [mol/s]
        m.B = Var(
            m.comps,
            domain=NonNegativeReals,
            doc='bottoms component flows [mol/s]',
            bounds=(0, max_flow),
            initialize=init['B'],
        )

        # Distillate component flows variable [mol/s]
        m.D = Var(
            m.comps,
            domain=NonNegativeReals,
            doc='distillate component flows [mol/s]',
            bounds=(0, max_flow),
            initialize=init['D'],
        )

        # Bottoms flow variable [mol/s]
        m.bot = Var(
            domain=NonNegativeReals,
            initialize=init['bot'],
            bounds=(0, 100),
            doc='bottoms flow [mol/s]',
        )

        # Distillate flow variable [mol/s]
        m.dis = Var(
            domain=NonNegativeReals,
            initialize=init['dis'],
            doc='distillate flow [mol/s]',
            bounds=(0, 100),
        )

        # Reflux ratio variable
        m.reflux_ratio = Var(
            domain=NonNegativeReals,
            bounds=(0.5, 4),
            doc='reflux ratio',
            initialize=init['reflux_ratio'],
        )

        # Reboil ratio variable
        m.reboil_ratio = Var(
            domain=NonNegativeReals,
            bounds=(1.3, 4),
            doc='reboil ratio',
            initialize=init['reboil_ratio'],
        )

        # Reflux fractions variable
        m.reflux_frac = Var(
            domain=NonNegativeReals,
            bounds=(0, 1 - 1e-6),
            doc='reflux fractions',
            initialize=init['reflux_frac'],
        )

        # Boilup fraction variable
        m.boilup_frac = Var(
            domain=NonNegativeReals,
            bounds=(0, 1 - 1e-6),
            doc='boilup fraction',
            initialize=init['boilup_frac'],
        )

        # Phase equilibrium constant variable
        m.Kc = Var(
            m.comps,
            m.trays,
            doc='Phase equilibrium constant',
            domain=NonNegativeReals,
            initialize=init['Kc'],
            bounds=(0, 1000),
        )

        # Temperature variable [K]
        m.T = Var(
            m.trays,
            doc='Temperature [K]',
            domain=NonNegativeReals,
            bounds=(min_T, max_T),
            initialize=init['T'],
        )

        # Pressure variable [bar]
        m.P = Var(doc='Pressure [bar]', bounds=(0, 5), initialize=init['P'])

        # Liquid activity coefficient variable
        m.gamma = Var(
            m.comps,
            m.trays,
            doc='liquid activity coefficent of component on tray',
            domain=NonNegativeReals,
            bounds=(0, 10),
            initialize=init['gamma'],
        )

        # Pure component vapor pressure variable [bar]
        m.Pvap = Var(
            m.comps,
            m.trays,
            doc='pure component vapor pressure of component on tray [bar]',
            domain=NonNegativeReals,
            bounds=(1e-3, 5),
            initialize=init['Pvap'],
        )

        # Variable related to fraction of critical temperature (1 - T/Tc)
        m.Pvap_X = Var(
            m.comps,
            m.trays,
            doc='Related to fraction of critical temperature (1 - T/Tc)',
            bounds=(0.25, 0.5),
            initialize=init['Pvap_X'],
        )

        # Liquid molar enthalpy variable [kJ/mol]
        m.H_L = Var(
            m.comps,
            m.trays,
            bounds=(0.1, 16),
            doc='Liquid molar enthalpy of component in tray [kJ/mol]',
            initialize=init['H_L'],
        )

        # Vapor molar enthalpy variable [kJ/mol]
        m.H_V = Var(
            m.comps,
            m.trays,
            bounds=(30, 16 + 40),
            doc='Vapor molar enthalpy of component in tray [kJ/mol]',
            initialize=init['H_V'],
        )

        # Component liquid molar enthalpy in feed variable [kJ/mol]
        m.H_L_spec_feed = Var(
            m.comps,
            doc='Component liquid molar enthalpy in feed [kJ/mol]',
            initialize=init['H_L_spec_feed'],
            bounds=(0.1, 16),
        )

        # Component vapor molar enthalpy in feed variable [kJ/mol]
        m.H_V_spec_feed = Var(
            m.comps,
            doc='Component vapor molar enthalpy in feed [kJ/mol]',
            initialize=init['H_V_spec_feed'],
            bounds=(30, 16 + 40),
        )

        # Reboiler duty variable [MJ/s]
        m.Qb = Var(
            domain=NonNegativeReals,
            doc='reboiler duty [MJ/s]',
            initialize=init['Qb'],
            bounds=(0, 8),
        )

        # Condenser duty variable [MJ/s]
        m.Qc = Var(
            domain=NonNegativeReals,
            doc='condenser duty [MJ/s]',
            initialize=init['Qc'],
            bounds=(0, 8),
        )

    else:
        # Feed temperature variable [K] within the minimum and maximum temperature
        m.T_feed = Var(
            doc='Feed temperature [K]',
            domain=NonNegativeReals,
            bounds=(min_T, max_T),
            initialize=368,  # Inlet temperature (95 C)
        )

        # Vapor fraction of the feed, value between 0 and 1
        m.feed_vap_frac = Var(doc='Vapor fraction of feed', initialize=0, bounds=(0, 1))

        # Total component feed flow [mol/s]
        m.feed = Var(m.comps, doc='Total component feed flow [mol/s]', initialize=50)

        # Liquid mole fraction variable for each component on each tray
        m.x = Var(
            m.comps,
            m.trays,
            doc='Liquid mole fraction',
            bounds=(0, 1),
            domain=NonNegativeReals,
            initialize=0.5,
        )

        # Vapor mole fraction variable for each component on each tray
        m.y = Var(
            m.comps,
            m.trays,
            doc='Vapor mole fraction',
            bounds=(0, 1),
            domain=NonNegativeReals,
            initialize=0.5,
        )

        # Component liquid flows from tray [mol/s]
        m.L = Var(
            m.comps,
            m.trays,
            doc='component liquid flows from tray in mol/s',
            domain=NonNegativeReals,
            bounds=(0, max_flow),
            initialize=50,
        )

        # Component vapor flows from tray [mol/s]
        m.V = Var(
            m.comps,
            m.trays,
            doc='component vapor flows from tray in mol/s',
            domain=NonNegativeReals,
            bounds=(0, max_flow),
            initialize=50,
        )

        # Liquid flows from each tray [mol/s]
        m.liq = Var(
            m.trays,
            domain=NonNegativeReals,
            doc='liquid flows from tray [mol/s]',
            initialize=100,
            bounds=(0, max_flow),
        )

        # Vapor flows from each tray [mol/s]
        m.vap = Var(
            m.trays,
            domain=NonNegativeReals,
            doc='vapor flows from tray [mol/s]',
            initialize=100,
            bounds=(0, max_flow),
        )

        # Bottoms component flows [mol/s]
        m.B = Var(
            m.comps,
            domain=NonNegativeReals,
            doc='bottoms component flows [mol/s]',
            bounds=(0, max_flow),
            initialize=50,
        )

        # Distillate component flows [mol/s]
        m.D = Var(
            m.comps,
            domain=NonNegativeReals,
            doc='distillate component flows [mol/s]',
            bounds=(0, max_flow),
            initialize=50,
        )

        # Bottoms flow [mol/s]
        m.bot = Var(
            domain=NonNegativeReals,
            initialize=50,
            bounds=(0, 100),
            doc='bottoms flow [mol/s]',
        )

        # Distillate flow [mol/s]
        m.dis = Var(
            domain=NonNegativeReals,
            initialize=50,
            doc='distillate flow [mol/s]',
            bounds=(0, 100),
        )

        # Reflux ratio variable
        m.reflux_ratio = Var(
            domain=NonNegativeReals, bounds=(0.5, 4), doc='reflux ratio', initialize=1.4
        )

        # Reboil ratio variable
        m.reboil_ratio = Var(
            domain=NonNegativeReals,
            bounds=(1.3, 4),
            doc='reboil ratio',
            initialize=0.9527,
        )

        # Reflux fractions variable
        m.reflux_frac = Var(
            domain=NonNegativeReals, bounds=(0, 1 - 1e-6), doc='reflux fractions'
        )

        # Boilup fraction variable
        m.boilup_frac = Var(
            domain=NonNegativeReals, bounds=(0, 1 - 1e-6), doc='boilup fraction'
        )

        # Phase equilibrium constant variable for each component on each tray
        m.Kc = Var(
            m.comps,
            m.trays,
            doc='Phase equilibrium constant',
            domain=NonNegativeReals,
            initialize=1,
            bounds=(0, 1000),
        )

        # Temperature variable for each tray [K]
        m.T = Var(
            m.trays,
            doc='Temperature [K]',
            domain=NonNegativeReals,
            bounds=(min_T, max_T),
        )

        # Pressure variable [bar]
        m.P = Var(doc='Pressure [bar]', bounds=(0, 5))

        # Liquid activity coefficient of component on tray variable
        m.gamma = Var(
            m.comps,
            m.trays,
            doc='liquid activity coefficent of component on tray',
            domain=NonNegativeReals,
            bounds=(0, 10),
            initialize=1,
        )

        # Pure component vapor pressure of component on tray [bar] variable
        m.Pvap = Var(
            m.comps,
            m.trays,
            doc='pure component vapor pressure of component on tray [bar]',
            domain=NonNegativeReals,
            bounds=(1e-3, 5),
            initialize=0.4,
        )

        # Variable related to fraction of critical temperature (1 - T/Tc)
        m.Pvap_X = Var(
            m.comps,
            m.trays,
            doc='Related to fraction of critical temperature (1 - T/Tc)',
            bounds=(0.25, 0.5),
            initialize=0.4,
        )

        # Liquid molar enthalpy of component in tray variable [kJ/mol]
        m.H_L = Var(
            m.comps,
            m.trays,
            bounds=(0.1, 16),
            doc='Liquid molar enthalpy of component in tray [kJ/mol]',
        )

        # Vapor molar enthalpy of component in tray variable [kJ/mol]
        m.H_V = Var(
            m.comps,
            m.trays,
            bounds=(30, 16 + 40),
            doc='Vapor molar enthalpy of component in tray [kJ/mol]',
        )

        # Component liquid molar enthalpy in feed variable [kJ/mol]
        m.H_L_spec_feed = Var(
            m.comps,
            doc='Component liquid molar enthalpy in feed [kJ/mol]',
            initialize=0,
            bounds=(0.1, 16),
        )

        # Component vapor molar enthalpy in feed variable [kJ/mol]
        m.H_V_spec_feed = Var(
            m.comps,
            doc='Component vapor molar enthalpy in feed [kJ/mol]',
            initialize=0,
            bounds=(30, 16 + 40),
        )

        # Reboiler duty variable [MJ/s]
        m.Qb = Var(
            domain=NonNegativeReals,
            doc='reboiler duty [MJ/s]',
            initialize=1,
            bounds=(0, 8),
        )

        # Condenser duty variable [MJ/s]
        m.Qc = Var(
            domain=NonNegativeReals,
            doc='condenser duty [MJ/s]',
            initialize=1,
            bounds=(0, 8),
        )

    m.partial_cond = Disjunct()  # Define a partial condenser disjunct
    m.total_cond = Disjunct()  # Define a total condenser disjunct
    m.condenser_choice = Disjunction(
        expr=[m.partial_cond, m.total_cond]
    )  # Condenser choice: partial or total condenser

    # Build mass balances for conditional trays, feed tray, condenser, and reboiler.
    for t in m.conditional_trays:
        _build_conditional_tray_mass_balance(m, t, m.tray[t], m.no_tray[t])
    _build_feed_tray_mass_balance(m)
    _build_condenser_mass_balance(m)
    _build_reboiler_mass_balance(m)

    @m.Constraint(m.comps, doc="Bottoms flow is equal to liquid leaving reboiler.")
    def bottoms_mass_balance(m, c):
        """Bottoms flow is equal to liquid leaving reboiler."""
        return m.B[c] == m.L[c, m.reboil_tray]

    @m.Constraint(
        doc="Boilup fraction is the ratio between the bottoms flow and the liquid leaving the reboiler."
    )
    def boilup_frac_defn(m):
        """Boilup fraction is the ratio between the bottoms flow and the liquid leaving the reboiler."""
        return m.bot == (1 - m.boilup_frac) * m.liq[m.reboil_tray + 1]

    @m.Constraint(
        doc="Reflux fraction is the ratio between the distillate flow and the difference in vapor flow in the condenser tray."
    )
    def reflux_frac_defn(m):
        """Reflux fraction is the ratio between the distillate flow and the difference in vapor flow in the condenser tray"""
        return m.dis == (1 - m.reflux_frac) * (
            m.vap[m.condens_tray - 1] - m.vap[m.condens_tray]
        )

    @m.Constraint(
        m.trays,
        doc="Total liquid flow on each tray is the sum of all component liquid flows on the tray.",
    )
    def liquid_sum(m, t):
        """Total liquid flow on each tray is the sum of all component liquid flows on the tray"""
        return sum(m.L[c, t] for c in m.comps) == m.liq[t]

    @m.Constraint(
        m.trays,
        doc="Total vapor flow on each tray is the sum of all component vapor flows on the tray.",
    )
    def vapor_sum(m, t):
        """Total vapor flow on each tray is the sum of all component vapor flows on the tray"""
        return sum(m.V[c, t] for c in m.comps) == m.vap[t]

    # Constraint to ensure the total bottoms flow is the sum of all component flows at the bottom
    m.bottoms_sum = Constraint(expr=sum(m.B[c] for c in m.comps) == m.bot)

    # Constraint to ensure the total distillate flow is the sum of all component flows at the top
    m.distil_sum = Constraint(expr=sum(m.D[c] for c in m.comps) == m.dis)

    @m.Constraint(
        m.trays,
        doc="Temperature decreases (or remains constant) from one tray to the next one down.",
    )
    def monotonoic_temperature(_, t):
        """Temperature decreases (or remains constant) from one tray to the next one down"""
        return m.T[t] >= m.T[t + 1] if t < max_trays else Constraint.Skip

    # Construct phase equilibrium relations for each conditional tray
    for t in m.conditional_trays:
        _build_tray_phase_equilibrium(m, t, m.tray[t])

    # Build blocks for the phase equilibrium relations on the feed tray, reboiler, and condenser
    m.feed_tray_phase_eq = Block()
    m.reboiler_phase_eq = Block()
    m.condenser_phase_eq = Block()

    # Construct phase equilibrium relations for the feed tray, reboiler, and condenser
    _build_tray_phase_equilibrium(m, m.feed_tray, m.feed_tray_phase_eq)
    _build_tray_phase_equilibrium(m, m.reboil_tray, m.reboiler_phase_eq)
    _build_tray_phase_equilibrium(m, m.condens_tray, m.condenser_phase_eq)

    # Construct heat relations for the column
    _build_column_heat_relations(m)

    @m.Constraint()
    def distillate_req(m):
        """Flow of benzene in the distillate meets the specified purity requirement"""
        return m.D['benzene'] >= m.distillate_purity * m.dis

    @m.Constraint()
    def bottoms_req(m):
        """Flow of toluene in the bottoms meets the specified purity requirement"""
        return m.B['toluene'] >= m.bottoms_purity * m.bot

    # Define the objective function for optimization.

    # The objective is to minimize the sum of condenser and reboiler duties, Qc and Qb, multiplied by 1E-3 [$/MJ/s] to convert from [MJ/s] to [$].
    # m.obj = Objective(expr=(m.Qc + m.Qb) * 1E-3, sense=minimize)

    # The objective is to minimize the sum of condenser and reboiler duties, Qc and Qb, multiplied by 1E3 to convert units,
    # and also the number of activated trays, which is obtained by summing up the indicator variables for the trays by 1E3 [$/No. of Trays].
    m.obj = Objective(
        expr=(m.Qc + m.Qb) * 1e3
        + (sum(m.tray[t].indicator_var for t in m.conditional_trays) + 1) * 1e3,
        sense=minimize,
    )

    # The objective is to minimize the number of activated trays, which is obtained by summing up the indicator variables for the trays.
    # m.obj = Objective(
    #     expr=sum(m.tray[t].indicator_var for t in m.conditional_trays) + 1)

    @m.Constraint()
    def reflux_ratio_calc(m):
        """Reflux ratio is the ratio between the distillate flow and the difference in vapor flow in the condenser tray."""
        return m.reflux_frac * (m.reflux_ratio + 1) == m.reflux_ratio

    @m.Constraint()
    def reboil_ratio_calc(m):
        """Reboil ratio is the ratio between the bottoms flow and the liquid leaving the reboiler."""
        return m.boilup_frac * (m.reboil_ratio + 1) == m.reboil_ratio

    @m.Constraint(m.conditional_trays)
    def tray_ordering(m, t):
        """Trays close to the feed should be activated first."""
        if t + 1 < m.condens_tray and t > m.feed_tray:
            return m.tray[t].indicator_var >= m.tray[t + 1].indicator_var
        elif t > m.reboil_tray and t + 1 < m.feed_tray:
            return m.tray[t + 1].indicator_var >= m.tray[t].indicator_var
        else:
            return Constraint.NoConstraint

    # Fix feed conditions
    m.feed['benzene'].fix(50)  # [mol/s]
    m.feed['toluene'].fix(50)  # [mol/s]
    m.T_feed.fix(368)  # [K]
    m.feed_vap_frac.fix(0.40395)  # [mol vapor / mol total]
    m.P.fix(1.01)  # [bar]
    # Fix to be total condenser
    m.partial_cond.deactivate()
    m.total_cond.indicator_var.fix(1)

    # -----------End of model declaration. These changes are required to run the DSDA--------------
    # FIX Indicator variables according to input
    ext_var_1 = x_input[0]  # Extract the first external variable from the input value.
    ext_var_2 = x_input[1]  # Extract the second external variable from the input value.

    if boolean_ref:
        # Boolean variables and intTrays set definition
        # Define the set of interior trays by excluding condensation and reboil trays from the total tray set.
        m.intTrays = Set(
            initialize=m.trays - [m.condens_tray, m.reboil_tray],
            doc='Interior trays of the column',
        )
        # Declare boolean variables for existence of boil-up and reflux flow in each stage,
        # as well as another boolean variable associated with tray presence and absence.
        m.YB = BooleanVar(m.intTrays, doc='Existence of boil-up flow in stage n')
        m.YR = BooleanVar(m.intTrays, doc='Existence of reflux flow in stage n')
        m.YP = BooleanVar(
            m.intTrays, doc='Boolean var associated with tray and no_tray'
        )
        m.YB_is_up = BooleanVar(
            doc='Ensure that there is no boil-up flow above the feed tray'
        )
        m.YR_is_down = BooleanVar(
            doc='Ensure that there is no reflux flow below the feed tray'
        )

        # Logical constraints
        @m.LogicalConstraint()
        def one_reflux(m):
            """Ensure that only one reflux flow exists."""
            return exactly(1, m.YR)

        @m.LogicalConstraint()
        def one_boilup(m):
            """Ensure that only one boil-up flow exists."""
            return exactly(1, m.YB)

        @m.LogicalConstraint()
        def boilup_fix(m):
            """Ensure that only one boil-up is happening."""
            return exactly(
                1, m.YB_is_up
            )  # NOTE: These are the one-dimensional Boolean variables, so it would have been enough to fix them to 1

        @m.LogicalConstraint()
        def reflux_fix(m):
            """Ensure that only one reflux is happening."""
            return exactly(
                1, m.YR_is_down
            )  # NOTE: These are the one-dimensional Boolean variables, so it would have been enough to fix them to 1

        @m.LogicalConstraint()
        def no_reflux_down(m):
            """Ensure no reflux flow exists below the feed tray."""
            return m.YR_is_down.equivalent_to(
                land(~m.YR[n] for n in range(m.reboil_tray + 1, m.feed_tray))
            )

        @m.LogicalConstraint()
        def no_boilup_up(m):
            """Ensure no boil-up flow exists above the feed tray."""
            return m.YB_is_up.equivalent_to(
                land(~m.YB[n] for n in range(m.feed_tray + 1, m.max_trays))
            )

        # Loop over internal trays
        for n in m.intTrays:
            # Fix the existence of reflux flow in stage n based on comparison with ext_var_1
            # ext_var_1 denotes the tray stage where reflux flow occurs.
            if n == ext_var_1:
                m.YR[n].fix(True)
            else:
                m.YR[n].fix(False)

            # Fix the existence of boil-up flow in stage n based on comparison with ext_var_2
            # ext_var_2 denotes the tray stage where boil-up flow occurs.
            if n == ext_var_2:
                m.YB[n].fix(True)
            else:
                m.YB[n].fix(False)

        # Check if the condition for all trays above the reboil tray holds,
        # and fix YR_is_down based on this
        temp = value(land(~m.YR[n] for n in range(m.reboil_tray + 1, m.feed_tray)))
        # temp is a temporary variable that stores the outcome of the logical condition for YR_is_down
        if temp == True:
            m.YR_is_down.fix(True)

        # Check if the condition for all trays above the feed tray holds,
        # and fix YB_is_up based on this
        temp = value(land(~m.YB[n] for n in range(m.feed_tray + 1, m.max_trays)))
        # temp is a temporary variable that stores the outcome of the logical condition for YB_is_up
        if temp == True:
            m.YB_is_up.fix(True)

        @m.LogicalConstraint(m.conditional_trays)
        def YP_or_notYP(m, n):
            """Define logical constraint YP_or_notYP based on conditions of reflux and boil-up flows above and including the tray"""
            return m.YP[n].equivalent_to(
                land(
                    lor(m.YR[j] for j in range(n, m.max_trays)),
                    lor(land(~m.YB[j] for j in range(n, m.max_trays)), m.YB[n]),
                )
            )

        # Loop over conditional trays
        # Associate YP variable with its corresponding binary variable
        for n in m.conditional_trays:
            m.YP[n].associate_binary_var(m.tray[n].indicator_var)

        # Loop over conditional trays again
        for n in m.conditional_trays:
            # Check the logical condition (similar to the one in YP_or_notYP)
            # temp is a temporary variable that stores the outcome of the logical condition for the current tray
            temp = value(
                land(
                    lor(m.YR[j] for j in range(n, m.max_trays)),
                    lor(land(~m.YB[j] for j in range(n, m.max_trays)), m.YB[n]),
                )
            )

            # Fix the binary variables based on the outcome of the logical condition
            if temp == True:
                m.tray[n].indicator_var.fix(True)
                m.no_tray[n].indicator_var.fix(False)
            else:
                m.tray[n].indicator_var.fix(False)
                m.no_tray[n].indicator_var.fix(True)

    else:
        # We are using a Boolean formulation because we know the Boolean variables and their relationships take longer for GAMS to write

        # Dictionaries to store the fixed values of YR and YB for each tray
        YR_fixed = {}
        YB_fixed = {}

        # Iterate over all trays except for the condenser and reboiler
        for n in m.trays - [m.condens_tray, m.reboil_tray]:
            # Fix the value of YR for the current tray based on comparison with ext_var_1
            if n == ext_var_1:
                YR_fixed[n] = 1
            else:
                YR_fixed[n] = 0

            # Fix the value of YB for the current tray based on comparison with ext_var_2
            if n == ext_var_2:
                YB_fixed[n] = 1
            else:
                YB_fixed[n] = 0

        # Iterate over all trays except for the condenser, reboiler, and feed trays
        for n in m.trays - [m.condens_tray, m.reboil_tray, m.feed_tray]:
            # Calculate the temporary variable based on the fixed values of YR and YB
            temp = (
                1
                - (
                    1
                    - sum(YR_fixed[j] for j in m.trays if j >= n and j <= max_trays - 1)
                )
                - (
                    sum(YB_fixed[j] for j in m.trays if j >= n and j <= max_trays - 1)
                    - YB_fixed[n]
                )
            )

            # Fix the indicator variables for the current tray based on the temporary variable, except for the feed tray
            if temp == 1 and n != m.feed_tray:
                m.tray[n].indicator_var.fix(True)
                m.no_tray[n].indicator_var.fix(False)
            elif temp == 0 and n != m.feed_tray:
                m.tray[n].indicator_var.fix(False)
                m.no_tray[n].indicator_var.fix(True)

    # Apply transformations to the model
    # Convert logical constraints to linear constraints
    TransformationFactory('core.logical_to_linear').apply_to(m)
    # Fix the disjuncts
    TransformationFactory('gdp.fix_disjuncts').apply_to(m)
    # Deactivate trivial constraints
    TransformationFactory('contrib.deactivate_trivial_constraints').apply_to(
        m, tmp=False, ignore_infeasible=True
    )

    # Set the status of the model and initialize an empty dictionary for the model's initialization
    m.dsda_status = 'Initialized'
    m.dsda_initialization = {}

    # Check equation feasibility
    try:
        fbbt(m)  # Apply feasibility-based bound tightening (FBBT) to the model

        # SOLVE
        if provide_init == False:
            initialize(m)  # Initialize the model if no initialization is provided

        # SOLVE
        # Get the path of the current file and the path where GAMS files will be stored
        dir_path = os.path.dirname(os.path.abspath(__file__))
        # Specify the path where GAMS files will be stored
        gams_path = os.path.join(dir_path, "gamsfiles/")
        if not (os.path.exists(gams_path)):
            print(
                'Directory for automatically generated files '
                + gams_path
                + ' does not exist. We will create it'
            )
            os.makedirs(gams_path)  # Create the directory if it doesn't exist

        # Specify the solver and its options
        solvername = 'gams'
        opt = SolverFactory(solvername, solver=nlp_solver)
        m.results = opt.solve(
            m,
            tee=True,
            # Uncomment the following lines if you want to save GAMS models
            # keepfiles=True,
            # tmpdir=gams_path,
            # symbolic_solver_labels=True,
            skip_trivial_constraints=True,
            add_options=[
                'option reslim = 120;'  # Limit the resource usage (time limit) [s]
                'option optcr = 0.0;'  # Set the relative optimality tolerance
                'option iterLim = 20000;'  # Set the maximum number of iterations
                # Uncomment the following lines to setup IIS computation of BARON through option file
                # 'GAMS_MODEL.optfile = 1;'
                # '\n'
                # '$onecho > baron.opt \n'
                # 'CompIIS 1 \n'
                # '$offecho'
                # 'display(execError);'
            ],
        )

        # Check the termination condition of the solver
        if (
            m.results.solver.termination_condition == 'locallyOptimal'
            or m.results.solver.termination_condition == 'optimal'
        ):
            m.dsda_status = 'Optimal'  # If the solver found an optimal or locally optimal solution, set the status to 'Optimal'
        elif m.results.solver.termination_condition == 'infeasible':
            m.dsda_status = 'Evaluated_Infeasible'  # If the solver found the problem to be infeasible, set the status to 'Evaluated_Infeasible'

        # Save results (for initialization)
        # Initialize empty dictionaries for storing the initialization values of various parameters
        (
            T_feed_init,
            feed_vap_frac_init,
            feed_init,
            x_init,
            y_init,
            L_init,
            V_init,
            liq_init,
            vap_init,
            B_init,
            D_init,
            bot_init,
            dis_init,
            reflux_ratio_init,
        ) = ({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {})
        (
            reboil_ratio_init,
            reflux_frac_init,
            boilup_frac_init,
            Kc_init,
            T_init,
            P_init,
            gamma_init,
        ) = ({}, {}, {}, {}, {}, {}, {})
        (
            Pvap_init,
            Pvap_X_init,
            H_L_init,
            H_V_init,
            H_L_spec_feed_init,
            H_V_spec_feed_init,
            Qb_init,
            Qc_init,
        ) = ({}, {}, {}, {}, {}, {}, {}, {})

        # Save the current value of different parameters in the model to the respective dictionaries
        T_feed_init = value(m.T_feed)
        feed_vap_frac_init = value(m.feed_vap_frac)
        bot_init = value(m.bot)
        dis_init = value(m.dis)
        reflux_ratio_init = value(m.reflux_ratio)
        reboil_ratio_init = value(m.reboil_ratio)
        reflux_frac_init = value(m.reflux_frac)
        boilup_frac_init = value(m.boilup_frac)
        Qb_init = value(m.Qb)
        Qc_init = value(m.Qc)
        P_init = value(m.P)

        # Loop over the components in the model, save the current value of feed, B, D, H_L_spec_feed, and H_V_spec_feed to the respective dictionaries
        for i in m.comps:
            feed_init[i] = value(m.feed[i])
            B_init[i] = value(m.B[i])
            D_init[i] = value(m.D[i])
            H_L_spec_feed_init[i] = value(m.H_L_spec_feed[i])
            H_V_spec_feed_init[i] = value(m.H_V_spec_feed[i])

        # Loop over the trays in the model, save the current value of liq, vap, and T to the respective dictionaries
        for n in m.trays:
            liq_init[n] = value(m.liq[n])
            vap_init[n] = value(m.vap[n])
            T_init[n] = value(m.T[n])

        # Loop over the components and trays in the model
        # Save the current value of parameters including mole fraction, liquid and vapor flow rates,
        # equilibrium constant, activity coefficient, vapor pressure, fugacity, and enthalpy of liquid and vapor phases to respective dictionaries
        for i in m.comps:
            for n in m.trays:
                x_init[i, n] = value(m.x[i, n])
                y_init[i, n] = value(m.y[i, n])
                L_init[i, n] = value(m.L[i, n])
                V_init[i, n] = value(m.V[i, n])
                Kc_init[i, n] = value(m.Kc[i, n])
                gamma_init[i, n] = value(m.gamma[i, n])
                Pvap_init[i, n] = value(m.Pvap[i, n])
                Pvap_X_init[i, n] = value(m.Pvap_X[i, n])
                H_L_init[i, n] = value(m.H_L[i, n])
                H_V_init[i, n] = value(m.H_V[i, n])

        # Save all the initialized dictionaries in the model for future initialization
        m.dsda_initialization = {
            'T_feed': T_feed_init,
            'feed_vap_frac': feed_vap_frac_init,
            'feed': feed_init,
            'x': x_init,
            'y': y_init,
            'L': L_init,
            'V': V_init,
            'liq': liq_init,
            'vap': vap_init,
            'B': B_init,
            'D': D_init,
            'bot': bot_init,
            'dis': dis_init,
            'reflux_ratio': reflux_ratio_init,
            'reboil_ratio': reboil_ratio_init,
            'reflux_frac': reflux_frac_init,
            'boilup_frac': boilup_frac_init,
            'Kc': Kc_init,
            'T': T_init,
            'P': P_init,
            'gamma': gamma_init,
            'Pvap': Pvap_init,
            'Pvap_X': Pvap_X_init,
            'H_L': H_L_init,
            'H_V': H_V_init,
            'H_L_spec_feed': H_L_spec_feed_init,
            'H_V_spec_feed': H_V_spec_feed_init,
            'Qb': Qb_init,
            'Qc': Qc_init,
        }

        # print('timer',time.process_time()-t_start)

        # print(m.results.solver.termination_condition)
    except InfeasibleConstraintException:
        m.dsda_status = 'FBBT_Infeasible'
        # print('Try an infeasible')

    return m


# ---------Other functions to define the model-------------------------------------------------
# Function for creating mass balances and flow rates for each component in a distillation column
def _build_conditional_tray_mass_balance(m, t, tray, no_tray):
    """
    Builds the constraints for mass balance, liquid and vapor composition for a given tray (t) in the distillation column.
    The constraints model the behavior of the mass balance for different components on a tray, accounting for feed, vapor, and liquid flows,
    as well as special conditions for the feed tray, condenser, and reboiler. Additional constraints define the liquid and vapor composition
    on the tray, as well as conditions for when the tray does not exist.

    Args:
        m (pyomo.ConcreteModel): The model object containing the relevant variables and parameters.
        t (int): Tray number for which the constraints are being defined (integer).
        tray (pyomo.Disjunct): Disjunct object representing the case when the tray exists in the column.
        no_tray (pyomo.Disjunct): Disjunct object representing the case when the tray is absent in the column.

    Return:
        None. The function adds constraints to the model but does not return a value.
    """

    @tray.Constraint(m.comps)
    def mass_balance(_, c):
        """Mass balance on each component on a tray."""

        # Compute mass balance components
        return (
            (
                m.feed[c] if t == m.feed_tray else 0
            )  # Total feed to the tray if it's a feed tray
            - m.V[c, t]  # Total vapor flow from the tray
            - (
                m.D[c] if t == m.condens_tray else 0
            )  # Total distillate if it's a condenser tray
            + (
                m.L[c, t + 1] if t < m.condens_tray else 0
            )  # Liquid flow from the tray above if it's not a condenser
            - (
                m.B[c] if t == m.reboil_tray else 0
            )  # Total bottoms if it's a reboiler tray
            - (
                m.L[c, t] if t > m.reboil_tray else 0
            )  # Liquid flow to the tray below if it's not a reboiler
            + (
                m.V[c, t - 1] if t > m.reboil_tray else 0
            )  # Vapor flow from the tray below if it's not a reboiler
        ) == 0

    @tray.Constraint(m.comps)
    def tray_liquid_composition(_, c):
        """Liquid composition constraint for the tray"""
        return m.L[c, t] == m.liq[t] * m.x[c, t]

    @tray.Constraint(m.comps)
    def tray_vapor_compositions(_, c):
        """Vapor composition constraint for the tray"""
        return m.V[c, t] == m.vap[t] * m.y[c, t]

    @no_tray.Constraint(m.comps)
    def liq_comp_pass_through(_, c):
        """If there is no tray, liquid composition from the tray above is passed directly to the tray below."""
        return m.x[c, t] == m.x[c, t + 1]

    # If there is no tray, liquid flow from the tray above is passed directly to the tray below
    @no_tray.Constraint(m.comps)
    def liq_flow_pass_through(_, c):
        """Liquid flow rate constraint for the case when the tray does not exist"""
        return m.L[c, t] == m.L[c, t + 1]

    # If there is no tray, vapor composition from the tray below is passed directly to the tray above
    @no_tray.Constraint(m.comps)
    def vap_comp_pass_through(_, c):
        """Vapor composition constraint for the case when the tray does not exist"""
        return m.y[c, t] == m.y[c, t - 1]

    # If there is no tray, vapor flow from the tray below is passed directly to the tray above
    @no_tray.Constraint(m.comps)
    def vap_flow_pass_through(_, c):
        """Vapor flow rate constraint for the case when the tray does not exist"""
        return m.V[c, t] == m.V[c, t - 1]


def _build_feed_tray_mass_balance(m):
    """
    Constructs the mass balance and composition constraints for the feed tray in the distillation column.

    This function defines the mass balance for each component in the feed tray, taking into account the feed,
    vapor, and liquid streams entering and leaving the tray. It also includes constraints for the liquid and vapor
    compositions on the feed tray.

    Args:
        m (pyomo.ConcreteModel): The model object containing variables and parameters related to the feed tray, components, and mass streams.

    Constraints:
        - feed_mass_balance: Ensures that the total mass in and out of the feed tray for each component is balanced.
        - feed_tray_liquid_composition: Defines the liquid composition for each component on the feed tray as a product of
          liquid flow rate and liquid mole fraction.
        - feed_tray_vapor_composition: Defines the vapor composition for each component on the feed tray as a product of
          vapor flow rate and vapor mole fraction.
    Returns:
        None: The function directly updates the model object, adding constraints to it.

    Example:
        Consider a model `dist_model` representing a distillation column with all the required sets, parameters,
        and variables defined. To add the feed tray mass balance and composition constraints to the model, call:

        _build_feed_tray_mass_balance(dist_model)

    Note:
        The function assumes that the feed enters only one specific tray in the column, known as the feed tray.
        The flow dynamics between the feed tray and its adjacent trays (above and below) play a crucial role in
        the distribution and separation of the components.
    """
    # Defining the tray as the feed tray for the mass balance
    t = m.feed_tray  # The feed tray number

    # Mass balance for each component on the feed tray
    @m.Constraint(m.comps)
    def feed_mass_balance(_, c):
        """Mass balance on each component on a tray."""
        return (
            m.feed[c]  # Feed into the tray
            - m.V[c, t]  # Vapor leaving from the tray
            + m.L[c, t + 1]  # Liquid coming from the tray above
            - m.L[c, t]  # Liquid flowing to the tray below
            + m.V[c, t - 1]  # Vapor coming from the tray below
            == 0
        )

    @m.Constraint(m.comps)
    def feed_tray_liquid_composition(_, c):
        """Liquid composition constraint for the feed tray"""
        return m.L[c, t] == m.liq[t] * m.x[c, t]

    @m.Constraint(m.comps)
    def feed_tray_vapor_composition(_, c):
        """Vapor composition on each component on a tray."""
        return m.V[c, t] == m.vap[t] * m.y[c, t]


def _build_condenser_mass_balance(m):
    """
    Constructs the mass balance equations for the condenser tray in a distillation column.
    This includes the definition of various constraints for mass balance, partial and total
    condensation, vapor and liquid composition, and distillate composition. It handles
    different scenarios, such as partial and total condensation, by defining corresponding
    constraints.

    Args:
        m (pyomo.ConcreteModel): A model object that includes the relevant variables, parameters,
        and expressions for the distillation process, such as component flows, liquid
        and vapor compositions, and energy expressions.

    Constraints:
        condenser_mass_balance: Ensures the conservation of mass for each component in the tray.
        condenser_liquid_composition: Relates liquid flow for each component to total liquid flow and mole fraction.
        condenser_vapor_composition: Relates vapor flow for each component to total vapor flow and mole fraction.
        no_vapor_flow: Ensures no vapor flow for each component in total condensation.
        no_total_vapor_flow: Ensures total vapor flow is zero in total condensation.
        liquid_fraction_pass_through: Ensures liquid mole fractions are equal to vapor phase mole fractions from the tray below in total condensation.
        condenser_distillate_composition: Defines the flow of each component in the distillate.

    Return:
        None: The function adds the constraints directly to the model object and does not
        return a value.
    """
    t = m.condens_tray  # The condenser tray number

    # Mass balance for each component in the condenser tray
    @m.Constraint(m.comps)
    def condenser_mass_balance(_, c):
        """Mass balance for each component in the condenser tray."""
        return (
            -m.V[c, t]  # Vapor leaving from the tray
            - m.D[c]  # Loss to distillate
            - m.L[c, t]  # Liquid flowing to the tray below
            + m.V[c, t - 1]  # Vapor coming from the tray below
            == 0
        )

    @m.partial_cond.Constraint(m.comps)
    def condenser_liquid_composition(_, c):
        """This is the composition of the liquid phase in the condenser tray for partial condensation"""
        return m.L[c, t] == m.liq[t] * m.x[c, t]

    @m.partial_cond.Constraint(m.comps)
    def condenser_vapor_composition(_, c):
        """This is the composition of the vapor phase in the condenser tray for partial condensation"""
        return m.V[c, t] == m.vap[t] * m.y[c, t]

    @m.total_cond.Constraint(m.comps)
    def no_vapor_flow(_, c):
        """This function ensures no vapor flow for each component in the case of total condensation."""
        return m.V[c, t] == 0

    @m.total_cond.Constraint()
    def no_total_vapor_flow(_):
        """This function ensures no total vapor flow in the case of total condensation."""
        return m.vap[t] == 0

    @m.total_cond.Constraint(m.comps)
    def liquid_fraction_pass_through(_, c):
        """This function ensures the mole fractions of the liquid phase are equal to those of the vapor phase from the tray below"""
        return m.x[c, t] == m.y[c, t - 1]

    @m.Constraint(m.comps)
    def condenser_distillate_composition(_, c):
        """This function defines the composition of the distillate in the condenser tray"""
        return m.D[c] == m.dis * m.x[c, t]


def _build_reboiler_mass_balance(m):
    """
    Constructs the mass balance equations for the reboiler tray in a distillation column.
    This includes the definition of constraints for mass balance, and liquid and vapor composition.
    The function considers the molar flow rates to establish the mass balance in the reboiler.

    Args:
        m (pyomo.ConcreteModel): A model object that includes the relevant variables, parameters,
        and expressions for the distillation process, such as component flows, liquid
        and vapor compositions.

    Constraints:
        reboiler_mass_balance: Ensures the conservation of mass for each component in the tray.
        reboiler_liquid_composition: Relates liquid flow for each component to total liquid flow and mole fraction.
        reboiler_vapor_composition: Relates vapor flow for each component to total vapor flow and mole fraction.

    Return:
        None: The function adds the constraints directly to the model object and does not
        return a value.
    """
    t = m.reboil_tray  # The reboiler tray number

    @m.Constraint(m.comps)
    def reboiler_mass_balance(_, c):
        """This function defines the mass balance for each component in the reboiler tray"""
        t = m.reboil_tray
        """The equation considers the mass balance in terms of molar flow rates"""
        return (
            -m.V[c, t]  # Vapor leaving the tray
            + m.L[c, t + 1]  # Liquid coming from the tray above
            - m.B[c]  # Loss to bottoms
            == 0
        )

    @m.Constraint(m.comps)
    def reboiler_liquid_composition(_, c):
        """This function defines the composition of the liquid phase in the reboiler tray"""
        return m.L[c, t] == m.liq[t] * m.x[c, t]

    @m.Constraint(m.comps)
    def reboiler_vapor_composition(_, c):
        """This function defines the composition of the vapor phase in the reboiler tray"""
        return m.V[c, t] == m.vap[t] * m.y[c, t]


def _build_tray_phase_equilibrium(m, t, tray):
    """
    Constructs the phase equilibrium equations for a given tray in a distillation column.
    The function encompasses defining Raoult's law for each component in the tray,
    the relationship between the phase equilibrium constant and activity coefficient,
    the relationship between vapor pressure and temperature for each component,
    and calculates the activity coefficient for each component.

    Args:
        m (pyomo.ConcreteModel): A model object containing the relevant variables, parameters,
            and expressions for the distillation process, such as vapor pressure constants,
            phase equilibrium constants, and activity coefficients.
        t (int): Tray index representing the specific tray within the column.
        tray (pyomo.Block): A Pyomo Block object representing the specific tray within the model,
            where constraints related to the tray are to be added.

    Constraints:
        raoults_law: Defines Raoult's law for each component in the tray.
        phase_equil_const: Defines the relationship between the phase equilibrium constant and the activity coefficient.
        vapor_pressure: Defines the relationship between vapor pressure and temperature for each component.
        activity_coefficient: Calculates the activity coefficient for each component.

    Return:
        None: The function adds constraints directly to the model object and does not
        return a value.

    Note:
        The function applies the Antoine equation to relate vapor pressure to composition,
        and in this particular case, assumes the activity coefficient to be 1 for simplicity.
    """

    @tray.Constraint(m.comps)
    def raoults_law(_, c):
        """Raoult's law for each component in a tray."""
        return m.y[c, t] == m.x[c, t] * m.Kc[c, t]

    @tray.Constraint(m.comps)
    def phase_equil_const(_, c):
        """This function defines the relationship between the phase equilibrium constant and the activity coefficient"""
        return m.Kc[c, t] * m.P == (m.gamma[c, t] * m.Pvap[c, t])

    # This function defines the relationship between vapor pressure and temperature for each component in a tray
    @tray.Constraint(m.comps)
    def Pvap_relation(_, c):
        """The equation uses the Antoine equation to relate the vapor pressure to one minus the reduced temperature (represented by the variable x)"""
        k = m.pvap_const[c]
        x = m.Pvap_X[c, t]
        return (log(m.Pvap[c, t]) - log(k['Pc'])) * (1 - x) == (
            k['A'] * x + k['B'] * x**1.5 + k['C'] * x**3 + k['D'] * x**6
        )

    @tray.Constraint(m.comps)
    def Pvap_X_defn(_, c):
        """Defines the relationship between the one minus the reduced temperature variable (Pvap_X) for each component in a tray,
        and the actual temperature of the tray, normalized by the critical temperature of the component (Tc).
        Pvap_X is defined as the difference between 1 and the ratio of the tray temperature to
        the critical temperature of the component, representing the deviation from the critical temperature.
        """
        k = m.pvap_const[c]
        return m.Pvap_X[c, t] == 1 - m.T[t] / k['Tc']

    # This function calculates the activity coefficient for each component in a tray
    @tray.Constraint(m.comps)
    def gamma_calc(_, c):
        """For simplicity, the activity coefficient is assumed to be 1 in this case"""
        return m.gamma[c, t] == 1


def _build_column_heat_relations(m):
    """
    This function calculates the enthalpy of the liquid phase for each component in each tray.
    Constructs the enthalpy relations for both liquid and vapor phases in each tray of a distillation column.
    It calculates the enthalpy of the liquid and vapor phases based on the heat capacity coefficients and the
    temperature difference from a reference temperature. Additionally, it builds energy balances for various
    trays including conditional trays, feed tray, condenser, and reboiler.

    Args:
        m (pyomo.ConcreteModel): A model object that includes relevant variables, parameters, and expressions for
        the distillation process, such as heat capacity coefficients, temperature references, and heat of
        vaporization for different components.

    Expressions:
        liq_enthalpy_expr: Calculates the enthalpy of the liquid phase for each component in each tray based
        on heat capacity coefficients and temperature difference [kJ/mol].
        vap_enthalpy_expr: Calculates the enthalpy of the vapor phase for each component in each tray based
        on heat of vaporization, heat capacity coefficients, and temperature difference [kJ/mol].

    Calls:
        _build_conditional_tray_energy_balance: Constructs energy balance for conditional trays.
        _build_feed_tray_energy_balance: Constructs energy balance for feed tray.
        _build_condenser_energy_balance: Constructs energy balance for the condenser.
        _build_reboiler_energy_balance: Constructs energy balance for the reboiler.

    Return:
        None: The function adds the expressions and calls the related functions directly to the model object,
        and does not return a value.
    """

    @m.Expression(m.trays, m.comps)
    def liq_enthalpy_expr(_, t, c):
        k = m.liq_Cp_const[c]
        """The equation calculates the enthalpy based on the heat capacity coefficients and the temperature difference from a reference temperature [kJ/mol]"""
        return (
            k['A'] * (m.T[t] - m.T_ref)
            + k['B'] * (m.T[t] ** 2 - m.T_ref**2) / 2
            + k['C'] * (m.T[t] ** 3 - m.T_ref**3) / 3
            + k['D'] * (m.T[t] ** 4 - m.T_ref**4) / 4
            + k['E'] * (m.T[t] ** 5 - m.T_ref**5) / 5
        ) * 1e-6  # Convert from [J/mol] to [MJ/mol]

    # This function calculates the enthalpy of the vapor phase for each component in each tray
    @m.Expression(m.trays, m.comps)
    def vap_enthalpy_expr(_, t, c):
        k = m.vap_Cp_const[c]
        """The equation calculates the enthalpy based on the heat of vaporization and the heat capacity coefficients, 
        as well as the temperature difference from a reference temperature [kJ/mol]"""
        return (
            m.dH_vap[c]
            + k['A'] * (m.T[t] - m.T_ref)
            + k['B'] * (m.T[t] ** 2 - m.T_ref**2) / 2
            + k['C'] * (m.T[t] ** 3 - m.T_ref**3) / 3
            + k['D'] * (m.T[t] ** 4 - m.T_ref**4) / 4
            + k['E'] * (m.T[t] ** 5 - m.T_ref**5) / 5
        ) * 1e-3  # Convert the results from [J/mol] to [kJ/mol]

    # Energy balance constraints for each tray
    for t in m.conditional_trays:
        _build_conditional_tray_energy_balance(m, t, m.tray[t], m.no_tray[t])
    _build_feed_tray_energy_balance(m)
    _build_condenser_energy_balance(m)
    _build_reboiler_energy_balance(m)


def _build_conditional_tray_energy_balance(m, t, tray, no_tray):
    """
    Constructs the energy balance constraints for a specific tray in a distillation column, considering both
    active and inactive (pass-through) scenarios. The function includes constraints for balancing the energy
    in the specified tray, accounting for liquid and vapor enthalpies, and also considers cases where the tray
    is bypassed (no liquid or vapor contact).

    Args:
        m (Model Object): A model object containing the relevant variables, parameters, and expressions
            for the distillation process, such as liquid and vapor enthalpy expressions, and components.

        t (integer): The index of the tray for which the energy balance is being constructed.

        tray (Block Object): A block object representing the active scenario where the tray is in operation.

        no_tray (Block Object): A block object representing the scenario where the tray is bypassed
            (pass-through without liquid or vapor contact).

    Constraints:
        energy_balance: Balances the energy in the tray, accounting for liquid and vapor enthalpies,
            and also considers cases where the tray is bypassed (no liquid or vapor contact).
        liquid_enthalpy: Calculates the enthalpy of the liquid phase for each component in each tray
            based on heat capacity coefficients and temperature difference [kJ/mol].
        vapor_enthalpy: Calculates the enthalpy of the vapor phase for each component in each tray
            based on heat of vaporization, heat capacity coefficients, and temperature difference [kJ/mol].


    Return:
        None: The function adds constraints directly to the model object and does not return a value.

    Note:
        The energy balance includes a scaling factor of 1E-3, to match units [kJ/mol].
        The active tray scenario involves energy balance equations for liquid and vapor leaving the tray.
        The pass-through scenario includes constraints that ensure the enthalpy of liquid and vapor is
        unchanged across the inactive tray.
    """

    @tray.Constraint()
    def energy_balance(_):
        """Ensuring net heat in tray is zero (equilibrium)"""
        return (
            sum(
                m.L[c, t + 1] * m.H_L[c, t + 1]  # heat of liquid from tray above
                - m.L[c, t] * m.H_L[c, t]  # heat of liquid to tray below
                + m.V[c, t - 1] * m.H_V[c, t - 1]  # heat of vapor from tray below
                - m.V[c, t] * m.H_V[c, t]  # heat of vapor to tray above
                for c in m.comps
            )
            * 1e-3  # Convert the result from [J/mol] to [kJ/mol]
            == 0
        )

    @tray.Constraint(m.comps)
    def liq_enthalpy_calc(_, c):
        """Liquid enthalpy as the function of Temperature"""
        return m.H_L[c, t] == m.liq_enthalpy_expr[t, c]

    @tray.Constraint(m.comps)
    def vap_enthalpy_calc(_, c):
        """Enthalpy of vapor leaving a tray"""
        return m.H_V[c, t] == m.vap_enthalpy_expr[t, c]

    # In case the tray does not exist, pass the enthalpy values through to the next tray
    @no_tray.Constraint(m.comps)
    def liq_enthalpy_pass_through(_, c):
        """Pass through liquid enthalpy"""
        return m.H_L[c, t] == m.H_L[c, t + 1]

    @no_tray.Constraint(m.comps)
    def vap_enthalpy_pass_through(_, c):
        """Pass through vapor enthalpy"""
        return m.H_V[c, t] == m.H_V[c, t - 1]


def _build_feed_tray_energy_balance(m):
    """
    Constructs the energy balance constraints for the feed tray in a distillation column. The function includes
    constraints for balancing the energy in the feed tray, accounting for both liquid and vapor enthalpies.
    It defines how the enthalpies for liquid and vapor feed are calculated, considering the temperature of the feed
    and the corresponding heat capacity coefficients.

    Args:
        m (Model Object): A model object containing the relevant variables, parameters, and expressions
            for the distillation process, such as liquid and vapor enthalpy expressions, feed vapor fraction,
            components, feed temperature, and feed tray identifier.

    Constraints:
        - feed_tray_energy_balance: Ensures that the net heat in the feed tray is zero, considering heat from liquid and vapor streams.
        - feed_tray_liq_enthalpy_calc: Calculates liquid enthalpy in the feed tray based on temperature.
        - feed_tray_vap_enthalpy_calc: Calculates vapor enthalpy in the feed tray based on temperature.
        - feed_liq_enthalpy_calc: Constraint for feed liquid enthalpy based on feed temperature.
        - feed_vap_enthalpy_calc: Constraint for feed vapor enthalpy based on feed temperature.

    Expressions:
        - feed_liq_enthalpy_expr: Defines the feed liquid enthalpy as a function of feed temperature.
        - feed_vap_enthalpy_expr: Defines the feed vapor enthalpy as a function of feed temperature.

    Return:
        None: The function adds constraints and expressions directly to the model object and does not return a value.

    Note:
        The energy balance includes a scaling factor of 1E-3, to match units [kJ/mol].
        Specific enthalpy expressions for liquid and vapor feed are created within this function.
    """
    t = m.feed_tray

    @m.Constraint()
    def feed_tray_energy_balance(_):
        """Energy balance for the feed tray"""
        return (
            sum(
                m.feed[c]
                * (
                    m.H_L_spec_feed[c] * (1 - m.feed_vap_frac)
                    + m.H_V_spec_feed[c] * m.feed_vap_frac  # Heat of feed liquid
                )  # Heat of feed vapor
                for c in m.comps
            )
            + sum(
                m.L[c, t + 1] * m.H_L[c, t + 1]  # Heat of liquid from tray above
                - m.L[c, t] * m.H_L[c, t]  # Heat of liquid to tray below
                + m.V[c, t - 1] * m.H_V[c, t - 1]  # Heat of vapor from tray below
                - m.V[c, t] * m.H_V[c, t]  # Heat of vapor to tray above
                for c in m.comps
            )
        ) * (
            1e-3  # Convert the result from [kJ/mol] to [MJ/mol]
        ) == 0

    @m.Constraint(m.comps)
    def feed_tray_liq_enthalpy_calc(_, c):
        """Enthalpy of liquid from feed tray"""
        return m.H_L[c, t] == m.liq_enthalpy_expr[t, c]

    @m.Constraint(m.comps)
    def feed_tray_vap_enthalpy_calc(_, c):
        """Enthalpy of vapor from feed tray"""
        return m.H_V[c, t] == m.vap_enthalpy_expr[t, c]

    @m.Expression(m.comps)
    def feed_liq_enthalpy_expr(_, c):
        """Enthalpy of liquid feed"""
        k = m.liq_Cp_const[c]
        return (
            k['A'] * (m.T_feed - m.T_ref)
            + k['B'] * (m.T_feed**2 - m.T_ref**2) / 2
            + k['C'] * (m.T_feed**3 - m.T_ref**3) / 3
            + k['D'] * (m.T_feed**4 - m.T_ref**4) / 4
            + k['E'] * (m.T_feed**5 - m.T_ref**5) / 5
        ) * 1e-6  # Convert the result from [J/mol] to [MJ/mol]

    @m.Constraint(m.comps)
    def feed_liq_enthalpy_calc(_, c):
        """Enthalpy of liquid feed"""
        return m.H_L_spec_feed[c] == m.feed_liq_enthalpy_expr[c]

    @m.Expression(m.comps)
    def feed_vap_enthalpy_expr(_, c):
        """Enthalpy of vapor feed"""
        k = m.vap_Cp_const[c]
        return (
            m.dH_vap[c]
            + k['A'] * (m.T_feed - m.T_ref)
            + k['B'] * (m.T_feed**2 - m.T_ref**2) / 2
            + k['C'] * (m.T_feed**3 - m.T_ref**3) / 3
            + k['D'] * (m.T_feed**4 - m.T_ref**4) / 4
            + k['E'] * (m.T_feed**5 - m.T_ref**5) / 5
        ) * 1e-3  # convert [J/mol] into [kJ/mol]

    @m.Constraint(m.comps)
    def feed_vap_enthalpy_calc(_, c):
        """The enthalpy of the vapor feed is calculated using the feed temperature and the vapor heat capacity"""
        return m.H_V_spec_feed[c] == m.feed_vap_enthalpy_expr[c]


def _build_condenser_energy_balance(m):
    """
    Constructs the energy balance for the condenser tray in a distillation column, accommodating both partial
    and total condensation scenarios. The function includes the energy balances that conserve energy by balancing
    input and output enthalpies for both cases. It also contains constraints for calculating liquid and vapor
    enthalpies on the condenser tray based on given expressions.

    Args:
        m (Model Object): A model object containing the relevant variables, parameters, and expressions
            for the distillation process, such as liquid and vapor enthalpy expressions, components,
            and the condenser tray identifier.

    Constraints:
        - partial_condenser_energy_balance: Ensures that the net heat in the partial condenser is zero, considering
          the heat contributions of liquid distillate, liquid to the tray below, vapor from the tray below, and vapor
          from the partial condenser.
        - total_condenser_energy_balance: Ensures that the net heat in the total condenser is zero, considering the
          heat contributions of liquid distillate, liquid to the tray below, and vapor from the tray below.
        - condenser_liq_enthalpy_calc: Calculates liquid enthalpy in the condenser based on temperature.
        - vap_enthalpy_calc: Calculates vapor enthalpy in the condenser based on temperature (only in the case of
          a partial condenser).

    Return:
        None: The function adds constraints directly to the model object and does not return a value.

    Note:
        The condenser energy balances include a scaling factor of 1E-3, to match units [kJ/mol].
        The constraints vary depending on whether it is a partial or total condenser.
    """
    t = m.condens_tray

    @m.partial_cond.Constraint()
    def partial_condenser_energy_balance(_):
        """The partial condenser energy balance is calculated using the enthalpy of the liquid distillate."""
        return (
            -m.Qc
            + sum(
                -m.D[c] * m.H_L[c, t]  # Enthalpy of liquid distillate
                - m.L[c, t] * m.H_L[c, t]  # Enthalpy of liquid to tray below
                + m.V[c, t - 1] * m.H_V[c, t - 1]  # Enthalpy of vapor from tray below
                - m.V[c, t] * m.H_V[c, t]  # Enthalpy of vapor from partial condenser
                for c in m.comps
            )
            * 1e-3  # Convert [kJ/s] into [MJ/s]
            == 0
        )  # Ensuring net heat in partial condenser is zero (equilibrium)

    @m.total_cond.Constraint()
    def total_condenser_energy_balance(_):
        """The total condenser energy balance is calculated using the enthalpy of the liquid distillate."""
        return (
            -m.Qc
            + sum(
                -m.D[c] * m.H_L[c, t]  # Enthalpy of liquid distillate
                - m.L[c, t] * m.H_L[c, t]  # Enthalpy of liquid to tray below
                + m.V[c, t - 1] * m.H_V[c, t - 1]  # Enthalpy of vapor from tray below
                for c in m.comps
            )
            * 1e-3  # Convert [kJ/s] into [MJ/s]
            == 0
        )  # Ensuring net heat in total condenser is zero (equilibrium)

    @m.Constraint(m.comps)
    def condenser_liq_enthalpy_calc(_, c):
        """The enthalpy of liquid from the partial condenser is calculated using the liquid Cp constants."""
        return m.H_L[c, t] == m.liq_enthalpy_expr[t, c]

    @m.partial_cond.Constraint(m.comps)
    def vap_enthalpy_calc(_, c):
        """The enthalpy of vapor from the partial condenser is calculated using the vapor Cp constants."""
        return m.H_V[c, t] == m.vap_enthalpy_expr[t, c]


def _build_reboiler_energy_balance(m):
    """
    Constructs the energy balance and enthalpy calculations for the reboiler tray in a distillation column.
    The function includes the reboiler energy balance that conserves energy by balancing input and output enthalpies.
    It also contains constraints for calculating liquid and vapor enthalpies on the reboiler tray based on given expressions.

    Args:
        m (Model Object): A model object containing the relevant variables, parameters, and expressions
            for the distillation process, such as liquid and vapor enthalpy expressions, components,
            and the reboiler tray identifier.

    Constraints:
        reboiler_energy_balance: The reboiler energy balance involves conserving energy by balancing input and output enthalpies.
        reboiler_liq_enthalpy_calc: The enthalpy of liquid from the reboiler is calculated using the liquid Cp constants.
        reboiler_vap_enthalpy_calc: The enthalpy of vapor from the reboiler is calculated using the vapor Cp constants.

    Return:
        None: The function adds constraints directly to the model object and does not return a value.

    Note:
        The reboiler energy balance includes a scaling factor of 1E-3, to match units [kJ/mol].
    """
    t = m.reboil_tray

    @m.Constraint()
    def reboiler_energy_balance(_):
        """Ensures that the net heat in the reboiler is zero (equilibrium), considering the heat contributions from the
        liquid from the tray above, liquid bottoms, and vapor to the tray above."""
        return (
            m.Qb
            + sum(
                m.L[c, t + 1] * m.H_L[c, t + 1]  # Heat of liquid from tray above
                - m.B[c] * m.H_L[c, t]  # Heat of liquid bottoms on reboiler
                - m.V[c, t] * m.H_V[c, t]  # Heat of vapor to tray above
                for c in m.comps
            )
            * 1e-3  # Convert the result from [kJ/s] to [MJ/s]
            == 0
        )  # Ensuring net heat in reboiler is zero (equilibrium)

    @m.Constraint(m.comps)
    def reboiler_liq_enthalpy_calc(_, c):
        """This constraint sets the liquid enthalpy on the reboiler tray equal to the calculated liquid enthalpy from a given expression."""
        return m.H_L[c, t] == m.liq_enthalpy_expr[t, c]

    @m.Constraint(m.comps)
    def reboiler_vap_enthalpy_calc(_, c):
        """Constraint sets the vapor enthalpy on the reboiler tray equal to the calculated vapor enthalpy from a given expression."""
        return m.H_V[c, t] == m.vap_enthalpy_expr[t, c]


# Conditional checks if the script is being executed as the main program.
# If so, it runs the build_column function with the specified arguments.
if __name__ == "__main__":
    # Inputs
    NT = 17  # Total number of trays
    model_args = {
        'min_trays': 8,
        'max_trays': NT,
        'xD': 0.95,
        'xB': 0.95,
        'x_input': [16, 7],  # Reflux position  # Boilup position
        'nlp_solver': 'ipopth',
    }  # Model arguments
    m = build_column(
        **model_args
    )  # Calculate the fixed value of the external variables
    # NOTE: rename function to solve_fixed_externalvars_column
