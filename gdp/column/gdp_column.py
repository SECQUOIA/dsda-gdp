import math  # Provides functions for mathematical operations.
import os  # Provides functions for interacting with the operating system.

# Imports from the Pyomo library for building and solving optimization problems.
from pyomo.common.errors import InfeasibleConstraintException
from pyomo.contrib.fbbt.fbbt import fbbt
from pyomo.contrib.gdpopt.data_class import MasterProblemResult
from pyomo.core.base.misc import display
from pyomo.core.plugins.transform.logical_to_linear import \
    update_boolean_vars_from_binary  # Transforms logical constraints into binary constraints.
from pyomo.environ import (Block, BooleanVar, ConcreteModel, Constraint,
                           NonNegativeReals, Objective, Param, RangeSet, Set,
                           SolverFactory, TransformationFactory, Var, exactly,
                           land, log, lor, minimize, value)  # Core components of Pyomo's modeling language.
from pyomo.gdp import Disjunct, Disjunction  # Imports for disjunctive programming.
from pyomo.opt import SolutionStatus, SolverResults  # Classes to represent solver results.
from pyomo.opt import TerminationCondition as tc  # Enum of possible termination conditions for solvers.
from pyomo.opt.base.solvers import SolverFactory # Base class for solver factories.


def build_column(min_trays, max_trays, xD, xB):
    """Builds the column model."""
    # References: A comparative study between GDP and NLP formulations for conceptual design of distillation columns (Ghouse et al., 2018)

    # This function defines a model for a benzene-toluene distillation column.
    # min_trays: Minimum number of trays in the column.
    # max_trays: Maximum number of trays in the column.
    # xD: Required distillate purity.
    # xB: Required bottoms purity.
    m = ConcreteModel('benzene-toluene column')
    m.comps = Set(initialize=['benzene', 'toluene']) # The components in the feed [mol/s]
    min_T, max_T = 300, 400
    m.T_ref = 298.15
    max_flow = 500
    m.max_trays = max_trays
    m.condens_tray = max_trays
    m.feed_tray = math.ceil((max_trays / 2))
    m.reboil_tray = 1
    m.distillate_purity = xD
    m.bottoms_purity = xB
    m.pvap_const = {
        'benzene': {'A': -6.98273, 'B': 1.33213, 'C': -2.62863,
                    'D': -3.33399, 'Tc': 562.2, 'Pc': 48.9},
        'toluene': {'A': -7.28607, 'B': 1.38091, 'C': -2.83433,
                    'D': -2.79168, 'Tc': 591.8, 'Pc': 41.0}}
    m.vap_Cp_const = {
        'benzene': {'A': -3.392E1, 'B': 4.739E-1, 'C': -3.017E-4,
                    'D': 7.130E-8, 'E': 0},
        'toluene': {'A': -2.435E1, 'B': 5.125E-1, 'C': -2.765E-4,
                    'D': 4.911E-8, 'E': 0}}
    m.liq_Cp_const = {
        'benzene': {'A': 1.29E5, 'B': -1.7E2, 'C': 6.48E-1,
                    'D': 0, 'E': 0},
        'toluene': {'A': 1.40E5, 'B': -1.52E2, 'C': 6.95E-1,
                    'D': 0, 'E': 0}}
    m.dH_vap = {'benzene': 33.770E3, 'toluene': 38.262E3}  # J/mol

    m.trays = RangeSet(max_trays, doc='Set of potential trays') # Define a set of trays in the column
    m.conditional_trays = Set(initialize=m.trays - [m.condens_tray, m.feed_tray, m.reboil_tray], doc="Trays that may be turned on and off.") # Define a set of trays that can be turned on and off
    m.tray = Disjunct(m.conditional_trays, doc='Disjunct for tray existence') # Define a disjunction for tray existence
    m.no_tray = Disjunct(m.conditional_trays, doc='Disjunct for tray absence') # Define a disjunction for tray absence

    # Define a disjunction function that can toggle between tray existence and absence
    @m.Disjunction(m.conditional_trays, doc='Tray exists or does not')
    def tray_no_tray(b, t):
        return [b.tray[t], b.no_tray[t]]
    m.minimum_num_trays = Constraint(expr=sum(m.tray[t].indicator_var for t in m.conditional_trays) + 1 >= min_trays) # Ensure minimum number of trays

    # Define variables
    m.T_feed = Var(doc='Feed temperature [K]', domain=NonNegativeReals, bounds=(min_T, max_T), initialize=368) # Feed temperature
    m.feed_vap_frac = Var(doc='Vapor fraction of feed', initialize=0, bounds=(0, 1)) # Vapor fraction of feed
    m.feed = Var(m.comps, doc='Total component feed flow [mol/s]', initialize=50) # Total component feed flow
    m.x = Var(m.comps, m.trays, doc='Liquid mole fraction', bounds=(0, 1), domain=NonNegativeReals, initialize=0.5) # Liquid mole fraction
    m.y = Var(m.comps, m.trays, doc='Vapor mole fraction', bounds=(0, 1), domain=NonNegativeReals, initialize=0.5) # Vapor mole fraction
    m.L = Var(m.comps, m.trays, doc='component liquid flows from tray in kmol', domain=NonNegativeReals, bounds=(0, max_flow), initialize=50) # Component liquid flows from tray
    m.V = Var(m.comps, m.trays, doc='component vapor flows from tray in kmol', domain=NonNegativeReals, bounds=(0, max_flow), initialize=50) # Component vapor flows from tray
    m.liq = Var(m.trays, domain=NonNegativeReals, doc='liquid flows from tray in kmol', initialize=100, bounds=(0, max_flow)) # Liquid flows from tray
    m.vap = Var(m.trays, domain=NonNegativeReals, doc='vapor flows from tray in kmol', initialize=100, bounds=(0, max_flow)) # Vapor flows from tray
    m.B = Var(m.comps, domain=NonNegativeReals, doc='bottoms component flows in kmol', bounds=(0, max_flow), initialize=50) # Bottoms component flows
    m.D = Var(m.comps, domain=NonNegativeReals, doc='distillate component flows in kmol', bounds=(0, max_flow), initialize=50) # Distillate component flows
    m.bot = Var(domain=NonNegativeReals, initialize=50, bounds=(0, 100), doc='bottoms flow in kmol') # Bottoms flow
    m.dis = Var(domain=NonNegativeReals, initialize=50, doc='distillate flow in kmol', bounds=(0, 100)) # Distillate flow
    m.reflux_ratio = Var(domain=NonNegativeReals, bounds=(0.5, 4), doc='reflux ratio', initialize=1.4) # Reflux ratio
    m.reboil_ratio = Var(domain=NonNegativeReals, bounds=(1.3, 4), doc='reboil ratio', initialize=0.9527) # Reboil ratio
    m.reflux_frac = Var(domain=NonNegativeReals, bounds=(0, 1 - 1E-6), doc='reflux fractions') # Reflux fractions
    m.boilup_frac = Var(domain=NonNegativeReals, bounds=(0, 1 - 1E-6), doc='boilup fraction') # Boilup fraction
    m.Kc = Var(m.comps, m.trays, doc='Phase equilibrium constant', domain=NonNegativeReals, initialize=1, bounds=(0, 1000)) # Phase equilibrium constant
    m.T = Var(m.trays, doc='Temperature [K]', domain=NonNegativeReals, bounds=(min_T, max_T)) # Tray temperature
    m.P = Var(doc='Pressure [bar]', bounds=(0, 5)) # Pressure
    m.gamma = Var(m.comps, m.trays, doc='liquid activity coefficent of component on tray', domain=NonNegativeReals, bounds=(0, 10), initialize=1) # Liquid activity coefficient
    m.Pvap = Var(m.comps, m.trays, doc='pure component vapor pressure of component on tray in bar', domain=NonNegativeReals, bounds=(1E-3, 5), initialize=0.4) # Pure component vapor pressure
    m.Pvap_rel = Var(m.comps, m.trays, doc='pure component relative vapor pressure of component on tray in bar (to avoid numerical problems)', domain=NonNegativeReals, bounds=(0, 5), initialize=0.4) # Pure component relative vapor pressure
    m.Pvap_X = Var(m.comps, m.trays, doc='Related to fraction of critical temperature (1 - T/Tc)', bounds=(0.25, 0.5), initialize=0.4) # Related to fraction of critical temperature
    m.H_L = Var(m.comps, m.trays, bounds=(0.1, 16), doc='Liquid molar enthalpy of component in tray (kJ/mol)') # Liquid molar enthalpy of component in tray
    m.H_V = Var(m.comps, m.trays, bounds=(30, 16 + 40), doc='Vapor molar enthalpy of component in tray (kJ/mol)') # Vapor molar enthalpy of component in tray
    m.H_L_spec_feed = Var(m.comps, doc='Component liquid molar enthalpy in feed [kJ/mol]', initialize=0, bounds=(0.1, 16)) # Component liquid molar enthalpy in feed
    m.H_V_spec_feed = Var(m.comps, doc='Component vapor molar enthalpy in feed [kJ/mol]', initialize=0, bounds=(30, 16 + 40)) # Component vapor molar enthalpy in feed
    m.Qb = Var(domain=NonNegativeReals, doc='reboiler duty (MJ/s)', initialize=1, bounds=(0, 8)) # Reboiler duty
    m.Qc = Var(domain=NonNegativeReals, doc='condenser duty (MJ/s)', initialize=1, bounds=(0, 8)) # Condenser duty


    m.partial_cond = Disjunct()  # Define a partial condenser disjunct
    m.total_cond = Disjunct()  # Define a total condenser disjunct
    m.condenser_choice = Disjunction(expr=[m.partial_cond, m.total_cond])  # Condenser choice: partial or total condenser

    # Build mass balance for conditional trays
    for t in m.conditional_trays:
        _build_conditional_tray_mass_balance(m, t, m.tray[t], m.no_tray[t])

    # Build mass balance for feed tray
    _build_feed_tray_mass_balance(m)

    # Build mass balance for condenser
    _build_condenser_mass_balance(m)

    # Build mass balance for reboiler
    _build_reboiler_mass_balance(m)

    # Constraint to ensure the bottoms flow is equal to the liquid leaving the reboiler
    @m.Constraint(m.comps, doc="Bottoms flow is equal to liquid leaving reboiler.")
    def bottoms_mass_balance(m, c):
        return m.B[c] == m.L[c, m.reboil_tray]

    # Constraint to define boilup fraction
    @m.Constraint()
    def boilup_frac_defn(m):
        return m.bot == (1 - m.boilup_frac) * m.liq[m.reboil_tray + 1]

    # Constraint to define reflux fraction
    @m.Constraint()
    def reflux_frac_defn(m):
        return m.dis == (1 - m.reflux_frac) * (m.vap[m.condens_tray - 1] - m.vap[m.condens_tray])

    # Constraint to ensure the sum of component liquid flows from each tray equals the total liquid flow from that tray
    @m.Constraint(m.trays)
    def liquid_sum(m, t):
        return sum(m.L[c, t] for c in m.comps) == m.liq[t]

    # Constraint to ensure the sum of component vapor flows from each tray equals the total vapor flow from that tray
    @m.Constraint(m.trays)
    def vapor_sum(m, t):
        return sum(m.V[c, t] for c in m.comps) == m.vap[t]

    # Constraint to ensure the sum of component bottoms flows equals the total bottoms flow
    m.bottoms_sum = Constraint(expr=sum(m.B[c] for c in m.comps) == m.bot)

    # Constraint to ensure the sum of component distillate flows equals the total distillate flow
    m.distil_sum = Constraint(expr=sum(m.D[c] for c in m.comps) == m.dis)

    # Constraint to ensure that the temperature on each tray is greater or equal to the temperature on the tray below
    @m.Constraint(m.trays)
    def monotonoic_temperature(_, t):
        return m.T[t] >= m.T[t + 1] if t < max_trays else Constraint.Skip


    # Building phase equilibrium for each conditional tray
    for t in m.conditional_trays:
        _build_tray_phase_equilibrium(m, t, m.tray[t])

    # Defining a block for feed tray, reboiler and condenser phase equilibrium
    m.feed_tray_phase_eq = Block()
    m.reboiler_phase_eq = Block()
    m.condenser_phase_eq = Block()

    # Building phase equilibrium for feed tray, reboiler and condenser
    _build_tray_phase_equilibrium(m, m.feed_tray, m.feed_tray_phase_eq)
    _build_tray_phase_equilibrium(m, m.reboil_tray, m.reboiler_phase_eq)
    _build_tray_phase_equilibrium(m, m.condens_tray, m.condenser_phase_eq)

    # Building heat relations for the column
    _build_column_heat_relations(m)

    # Constraint to ensure the distillate contains at least a certain purity of benzene
    @m.Constraint()
    def distillate_req(m):
        return m.D['benzene'] >= m.distillate_purity * m.dis

    # Constraint to ensure the bottoms contains at least a certain purity of toluene
    @m.Constraint()
    def bottoms_req(m):
        return m.B['toluene'] >= m.bottoms_purity * m.bot

    # Define the objective function as the sum of reboiler and condenser duty plus an indicator for tray activation 
    m.obj = Objective(expr=(m.Qc + m.Qb) * 1E3 + 1E3 * (sum(m.tray[t].indicator_var for t in m.conditional_trays) + 1), sense=minimize)

    # Constraint to calculate the reflux ratio
    @m.Constraint()
    def reflux_ratio_calc(m):
        return m.reflux_frac * (m.reflux_ratio + 1) == m.reflux_ratio

    # Constraint to calculate the reboil ratio
    @m.Constraint()
    def reboil_ratio_calc(m):
        return m.boilup_frac * (m.reboil_ratio + 1) == m.reboil_ratio

    # Constraint to ensure a specific tray order: trays close to the feed should be activated first
    @m.Constraint(m.conditional_trays)
    def tray_ordering(m, t):
        """Trays close to the feed should be activated first."""
        if t + 1 < m.condens_tray and t > m.feed_tray:
            return m.tray[t].indicator_var >= m.tray[t + 1].indicator_var
        elif t > m.reboil_tray and t + 1 < m.feed_tray:
            return m.tray[t + 1].indicator_var >= m.tray[t].indicator_var
        else:
            return Constraint.NoConstraint


    # Defining set of interior trays in the column (excluding condenser and reboiler trays)
    m.intTrays = Set(initialize=m.trays - [m.condens_tray, m.reboil_tray], doc='Interior trays of the column')

    # Defining Boolean variables to denote existence of boil-up flow and reflux flow in each interior tray
    m.YB = BooleanVar(m.intTrays, initialize=False, doc='Existence of boil-up flow in stage n')
    m.YR = BooleanVar(m.intTrays, initialize=False, doc='Existence of reflux flow in stage n')

    # Initializing at least one reflux and boilup tray to avoid errors in Mixed-Integer Nonlinear Programming (MINLP) solvers
    m.YB[m.reboil_tray+1].set_value(True)
    m.YR[m.max_trays-1].set_value(True)

    # Defining additional Boolean variables for logical constraints
    m.YP = BooleanVar(m.intTrays, doc='Boolean var associated with tray and no_tray')
    m.YB_is_up = BooleanVar(doc='Boolean var for intermediate sum determining if Boilup is above the feed')
    m.YR_is_down = BooleanVar(doc='Boolean var for intermediate sum determining if Reflux is below the feed')

    # Defining logical constraints to ensure only one reflux and one boilup
    @m.LogicalConstraint()
    def one_reflux(m):
        return exactly(1, m.YR)

    @m.LogicalConstraint()
    def one_boilup(m):
        return exactly(1, m.YB)

    # Defining logical constraint for YP Boolean variable
    @m.LogicalConstraint(m.conditional_trays)
    def YP_or_notYP(m, n):
        return m.YP[n].equivalent_to(land(lor(m.YR[j] for j in range(n, m.max_trays)), lor(land(~m.YB[j] for j in range(n, m.max_trays)), m.YB[n])))

    # Associating YP Boolean variable with tray activation
    for n in m.conditional_trays:
        m.YP[n].associate_binary_var(m.tray[n].indicator_var)

    # Fixing feed conditions
    m.feed['benzene'].fix(50)  # Fixing benzene flow in the feed at 50 mol/s
    m.feed['toluene'].fix(50)  # Fixing toluene flow in the feed at 50 mol/s
    m.T_feed.fix(368)  # Fixing feed temperature at 368 K
    m.feed_vap_frac.fix(0.40395)  # Fixing feed vapor fraction
    m.P.fix(1.01)  # Fixing pressure at 1.01 bar

    # Fixing the system to be a total condenser
    m.partial_cond.deactivate()  # Deactivating partial condenser
    m.total_cond.indicator_var.fix(1)  # Activating total condenser

    # Fixing auxiliary Boolean variables for logical position of boilup and reflux
    m.YB_is_up.fix(True)
    m.YR_is_down.fix(True)

    # Returning the model
    return m


# This function builds the constraints for mass balance, liquid and vapor composition for a given tray (t) in the column
def _build_conditional_tray_mass_balance(m, t, tray, no_tray):
    """
    t = tray number
    tray = tray exists disjunct
    no_tray = tray absent disjunct
    """
    
    # Define mass balance constraint for each component in the tray
    @tray.Constraint(m.comps)
    def mass_balance(_, c):
        return (
            # Include feed flow rate if the current tray is the feed tray
            (m.feed[c] if t == m.feed_tray else 0)
            # Subtract vapor flow rate leaving the current tray
            - m.V[c, t]
            # Subtract distillate flow rate if current tray is the condenser
            - (m.D[c] if t == m.condens_tray else 0)
            # Include liquid flow rate from the tray above if current tray is not the condenser
            + (m.L[c, t + 1] if t < m.condens_tray else 0)
            # Subtract bottoms flow rate if current tray is the reboiler
            - (m.B[c] if t == m.reboil_tray else 0)
            # Subtract liquid flow rate to the tray below if current tray is not the reboiler
            - (m.L[c, t] if t > m.reboil_tray else 0)
            # Include vapor flow rate from the tray below if current tray is not the reboiler
            + (m.V[c, t - 1] if t > m.reboil_tray else 0) == 0)

    # Define constraints for liquid composition on the tray
    @tray.Constraint(m.comps)
    def tray_liquid_composition(_, c):
        return m.L[c, t] == m.liq[t] * m.x[c, t]

    # Define constraints for vapor composition on the tray
    @tray.Constraint(m.comps)
    def tray_vapor_compositions(_, c):
        return m.V[c, t] == m.vap[t] * m.y[c, t]

    # If the tray does not exist, the liquid composition should be equal to that of the tray above
    @no_tray.Constraint(m.comps)
    def liq_comp_pass_through(_, c):
        return m.x[c, t] == m.x[c, t + 1]

    # If the tray does not exist, the liquid flow rate should be equal to that of the tray above
    @no_tray.Constraint(m.comps)
    def liq_flow_pass_through(_, c):
        return m.L[c, t] == m.L[c, t + 1]

    # If the tray does not exist, the vapor composition should be equal to that of the tray below
    @no_tray.Constraint(m.comps)
    def vap_comp_pass_through(_, c):
        return m.y[c, t] == m.y[c, t - 1]

    # If the tray does not exist, the vapor flow rate should be equal to that of the tray below
    @no_tray.Constraint(m.comps)
    def vap_flow_pass_through(_, c):
        return m.V[c, t] == m.V[c, t - 1]


# This function constructs the mass balance and composition constraints for the feed tray
def _build_feed_tray_mass_balance(m):
    t = m.feed_tray  # The feed tray number

    # Mass balance for each component on the feed tray
    @m.Constraint(m.comps)
    def feed_mass_balance(_, c):
        return (
            m.feed[c]        # Feed to the tray
            - m.V[c, t]      # Vapor from the tray
            + m.L[c, t + 1]  # Liquid from the tray above
            - m.L[c, t]      # Liquid to the tray below
            + m.V[c, t - 1]  # Vapor from the tray below
            == 0)

    # Liquid composition constraint for each component on the feed tray
    @m.Constraint(m.comps)
    def feed_tray_liquid_composition(_, c):
        return m.L[c, t] == m.liq[t] * m.x[c, t]

    # Vapor composition constraint for each component on the feed tray
    @m.Constraint(m.comps)
    def feed_tray_vapor_composition(_, c):
        return m.V[c, t] == m.vap[t] * m.y[c, t]

# This function constructs the mass balance and composition constraints for the condenser
def _build_condenser_mass_balance(m):
    t = m.condens_tray  # The condenser tray number

    # Mass balance for each component in the condenser
    @m.Constraint(m.comps)
    def condenser_mass_balance(_, c):
        return (
            - m.V[c, t]      # Vapor from the tray
            - m.D[c]         # Loss to distillate
            - m.L[c, t]      # Liquid to the tray below
            + m.V[c, t - 1]  # Vapor from the tray below
            == 0)

    # For a partial condenser:
    @m.partial_cond.Constraint(m.comps)
    def condenser_liquid_composition(_, c):
        return m.L[c, t] == m.liq[t] * m.x[c, t]

    @m.partial_cond.Constraint(m.comps)
    def condenser_vapor_composition(_, c):
        return m.V[c, t] == m.vap[t] * m.y[c, t]

    # For a total condenser:
    @m.total_cond.Constraint(m.comps)
    def no_vapor_flow(_, c):
        return m.V[c, t] == 0  # No vapor flow for a total condenser

    @m.total_cond.Constraint()
    def no_total_vapor_flow(_):
        return m.vap[t] == 0  # Total vapor flow is zero for a total condenser

    @m.total_cond.Constraint(m.comps)
    def liquid_fraction_pass_through(_, c):
        return m.x[c, t] == m.y[c, t - 1]  # Liquid composition is the same as vapor composition from the tray below

    @m.Constraint(m.comps)
    def condenser_distillate_composition(_, c):
        return m.D[c] == m.dis * m.x[c, t]  # Define distillate composition

# This function constructs the mass balance and composition constraints for the reboiler
def _build_reboiler_mass_balance(m):
    t = m.reboil_tray  # The reboiler tray number

    # Mass balance for each component in the reboiler
    @m.Constraint(m.comps)
    def reboiler_mass_balance(_, c):
        t = m.reboil_tray
        return (
            - m.V[c, t]      # Vapor from the tray
            + m.L[c, t + 1]  # Liquid from the tray above
            - m.B[c]         # Loss to bottoms
            == 0)

    # Liquid composition constraint for each component in the reboiler
    @m.Constraint(m.comps)
    def reboiler_liquid_composition(_, c):
        return m.L[c, t] == m.liq[t] * m.x[c, t]

    # Vapor composition constraint for each component in the reboiler
    @m.Constraint(m.comps)
    def reboiler_vapor_composition(_, c):
        return m.V[c, t] == m.vap[t] * m.y[c, t]


# This function constructs the phase equilibrium constraints for a given tray
def _build_tray_phase_equilibrium(m, t, tray):
    # Raoult's law for each component on the tray
    @tray.Constraint(m.comps)
    def raoults_law(_, c):
        return m.y[c, t] == m.x[c, t] * m.Kc[c, t]

    # Definition of the phase equilibrium constant for each component on the tray
    @tray.Constraint(m.comps)
    def phase_equil_const(_, c):
        return m.Kc[c, t] * m.P == m.gamma[c, t] * m.Pvap[c, t]

    # Definition of the relative vapor pressure for each component on the tray
    @tray.Constraint(m.comps)
    def Pvap_relative(_, c):
        return m.Pvap_rel[c, t] == m.Pvap[c, t] - m.Pvap[c, t].lb

    # Relation between the relative vapor pressure and temperature for each component on the tray
    @tray.Constraint(m.comps)
    def Pvap_relation(_, c):
        k = m.pvap_const[c]
        x = m.Pvap_X[c, t]
        return (log(m.Pvap_rel[c, t] + m.Pvap[c, t].lb) - log(k['Pc'])) * (1 - x) == (
            k['A'] * x + k['B'] * x ** 1.5 + k['C'] * x ** 3 + k['D'] * x ** 6)

    # Definition of the temperature-dependent part of the relative vapor pressure for each component on the tray
    @tray.Constraint(m.comps)
    def Pvap_X_defn(_, c):
        k = m.pvap_const[c]
        return m.Pvap_X[c, t] == 1 - m.T[t] / k['Tc']

    # Assumption of ideal solution (gamma = 1) for each component on the tray
    @tray.Constraint(m.comps)
    def gamma_calc(_, c):
        return m.gamma[c, t] == 1

# This function constructs the expressions for liquid and vapor enthalpy and the energy balance constraints for each tray
def _build_column_heat_relations(m):
    # Liquid enthalpy expression for each component on each tray
    @m.Expression(m.trays, m.comps)
    def liq_enthalpy_expr(_, t, c):
        k = m.liq_Cp_const[c]
        return (
            k['A'] * (m.T[t] - m.T_ref) + k['B'] * (m.T[t] ** 2 - m.T_ref ** 2) / 2 +
            k['C'] * (m.T[t] ** 3 - m.T_ref ** 3) / 3 + k['D'] * (m.T[t] ** 4 - m.T_ref ** 4) / 4 +
            k['E'] * (m.T[t] ** 5 - m.T_ref ** 5) / 5) * 1E-6

    # Vapor enthalpy expression for each component on each tray
    @m.Expression(m.trays, m.comps)
    def vap_enthalpy_expr(_, t, c):
        k = m.vap_Cp_const[c]
        return (
            m.dH_vap[c] + k['A'] * (m.T[t] - m.T_ref) +
            k['B'] * (m.T[t] ** 2 - m.T_ref ** 2) / 2 + k['C'] * (m.T[t] ** 3 - m.T_ref ** 3) / 3 +
            k['D'] * (m.T[t] ** 4 - m.T_ref ** 4) / 4 + k['E'] * (m.T[t] ** 5 - m.T_ref ** 5) / 5) * 1E-3

    # Energy balance constraints for each tray
    for t in m.conditional_trays:
        _build_conditional_tray_energy_balance(m, t, m.tray[t], m.no_tray[t])
    _build_feed_tray_energy_balance(m)
    _build_condenser_energy_balance(m)
    _build_reboiler_energy_balance(m)


def _build_conditional_tray_energy_balance(m, t, tray, no_tray):
    # Energy balance constraint for each tray
    @tray.Constraint()
    def energy_balance(_):
        return sum(
            m.L[c, t + 1] * m.H_L[c, t + 1]  # heat of liquid from tray above
            - m.L[c, t] * m.H_L[c, t]  # heat of liquid to tray below
            + m.V[c, t - 1] * m.H_V[c, t - 1]  # heat of vapor from tray below
            - m.V[c, t] * m.H_V[c, t]  # heat of vapor to tray above
            for c in m.comps) * 1E-3 == 0  # Ensuring net heat in tray is zero (equilibrium)

    # Constraints to calculate enthalpy for liquid and vapor based on temperature for each tray
    @tray.Constraint(m.comps)
    def liq_enthalpy_calc(_, c):
        return m.H_L[c, t] == m.liq_enthalpy_expr[t, c]  # Liquid enthalpy = f(T)

    @tray.Constraint(m.comps)
    def vap_enthalpy_calc(_, c):
        return m.H_V[c, t] == m.vap_enthalpy_expr[t, c]  # Vapor enthalpy = f(T)

    # In case the tray does not exist, pass the enthalpy values through to the next tray
    @no_tray.Constraint(m.comps)
    def liq_enthalpy_pass_through(_, c):
        return m.H_L[c, t] == m.H_L[c, t + 1]  # Pass through liquid enthalpy

    @no_tray.Constraint(m.comps)
    def vap_enthalpy_pass_through(_, c):
        return m.H_V[c, t] == m.H_V[c, t - 1]  # Pass through vapor enthalpy


def _build_feed_tray_energy_balance(m):
    t = m.feed_tray

    # Energy balance constraint for the feed tray
    @m.Constraint()
    def feed_tray_energy_balance(_):
        return (
            sum(m.feed[c] * (
                m.H_L_spec_feed[c] * (1 - m.feed_vap_frac) +  # Heat of feed liquid
                m.H_V_spec_feed[c] * m.feed_vap_frac)  # Heat of feed vapor
                for c in m.comps) +
            sum(
                m.L[c, t + 1] * m.H_L[c, t + 1]  # Heat of liquid from tray above
                - m.L[c, t] * m.H_L[c, t]  # Heat of liquid to tray below
                + m.V[c, t - 1] * m.H_V[c, t - 1]  # Heat of vapor from tray below
                - m.V[c, t] * m.H_V[c, t]  # Heat of vapor to tray above
                for c in m.comps)) * 1E-3 == 0  # Ensuring net heat in feed tray is zero (equilibrium)

    # Constraints to calculate enthalpy for liquid and vapor based on temperature for feed tray
    @m.Constraint(m.comps)
    def feed_tray_liq_enthalpy_calc(_, c):
        return m.H_L[c, t] == m.liq_enthalpy_expr[t, c]  # Liquid enthalpy = f(T)

    @m.Constraint(m.comps)
    def feed_tray_vap_enthalpy_calc(_, c):
        return m.H_V[c, t] == m.vap_enthalpy_expr[t, c]  # Vapor enthalpy = f(T)

    # Expressions to calculate feed liquid and vapor enthalpy based on feed temperature
    @m.Expression(m.comps)
    def feed_liq_enthalpy_expr(_, c):
        k = m.liq_Cp_const[c]
        return (
            k['A'] * (m.T_feed - m.T_ref) +
            k['B'] * (m.T_feed ** 2 - m.T_ref ** 2) / 2 +
            k['C'] * (m.T_feed ** 3 - m.T_ref ** 3) / 3 +
            k['D'] * (m.T_feed ** 4 - m.T_ref ** 4) / 4 +
            k['E'] * (m.T_feed ** 5 - m.T_ref ** 5) / 5) * 1E-6  # Feed liquid enthalpy = f(T_feed)

    @m.Constraint(m.comps)
    def feed_liq_enthalpy_calc(_, c):
        return m.H_L_spec_feed[c] == m.feed_liq_enthalpy_expr[c]  # Constraint for feed liquid enthalpy

    @m.Expression(m.comps)
    def feed_vap_enthalpy_expr(_, c):
        k = m.vap_Cp_const[c]
        return (
            m.dH_vap[c] +
            k['A'] * (m.T_feed - m.T_ref) +
            k['B'] * (m.T_feed ** 2 - m.T_ref ** 2) / 2 +
            k['C'] * (m.T_feed ** 3 - m.T_ref ** 3) / 3 +
            k['D'] * (m.T_feed ** 4 - m.T_ref ** 4) / 4 +
            k['E'] * (m.T_feed ** 5 - m.T_ref ** 5) / 5) * 1E-3  # m.H_V_spec_feed[c]

    @m.Constraint(m.comps)
    def feed_vap_enthalpy_calc(_, c):
        return m.H_V_spec_feed[c] == m.feed_vap_enthalpy_expr[c]  # Constraint for feed vapor enthalpy



def _build_condenser_energy_balance(m):
    t = m.condens_tray

    # Energy balance for partial condenser
    @m.partial_cond.Constraint()
    def partial_condenser_energy_balance(_):
        return -m.Qc + sum(
            - m.D[c] * m.H_L[c, t]  # heat of liquid distillate
            - m.L[c, t] * m.H_L[c, t]  # heat of liquid to tray below
            + m.V[c, t - 1] * m.H_V[c, t - 1]  # heat of vapor from tray below
            - m.V[c, t] * m.H_V[c, t]  # heat of vapor from partial condenser
            for c in m.comps) * 1E-3 == 0  # Ensuring net heat in partial condenser is zero (equilibrium)

    # Energy balance for total condenser
    @m.total_cond.Constraint()
    def total_condenser_energy_balance(_):
        return -m.Qc + sum(
            - m.D[c] * m.H_L[c, t]  # heat of liquid distillate
            - m.L[c, t] * m.H_L[c, t]  # heat of liquid to tray below
            + m.V[c, t - 1] * m.H_V[c, t - 1]  # heat of vapor from tray below
            for c in m.comps) * 1E-3 == 0  # Ensuring net heat in total condenser is zero (equilibrium)

    # Constraints to calculate enthalpy for liquid and vapor based on temperature for condenser
    @m.Constraint(m.comps)
    def condenser_liq_enthalpy_calc(_, c):
        return m.H_L[c, t] == m.liq_enthalpy_expr[t, c]  # Liquid enthalpy = f(T)

    @m.partial_cond.Constraint(m.comps)
    def vap_enthalpy_calc(_, c):
        return m.H_V[c, t] == m.vap_enthalpy_expr[t, c]  # Vapor enthalpy = f(T)


def _build_reboiler_energy_balance(m):
    t = m.reboil_tray

    # Energy balance for reboiler
    @m.Constraint()
    def reboiler_energy_balance(_):
        return m.Qb + sum(
            m.L[c, t + 1] * m.H_L[c, t + 1]  # Heat of liquid from tray above
            - m.B[c] * m.H_L[c, t]  # heat of liquid bottoms if reboiler
            - m.V[c, t] * m.H_V[c, t]  # heat of vapor to tray above
            for c in m.comps) * 1E-3 == 0  # Ensuring net heat in reboiler is zero (equilibrium)

    # Constraints to calculate enthalpy for liquid and vapor based on temperature for reboiler
    @m.Constraint(m.comps)
    def reboiler_liq_enthalpy_calc(_, c):
        return m.H_L[c, t] == m.liq_enthalpy_expr[t, c]  # Liquid enthalpy = f(T)

    @m.Constraint(m.comps)
    def reboiler_vap_enthalpy_calc(_, c):
        return m.H_V[c, t] == m.vap_enthalpy_expr[t, c]  # Vapor enthalpy = f(T)


if __name__ == "__main__":
    # Inputs
    NT = 17  # Total number of trays
    model_args = {'min_trays': 8, 'max_trays': NT, 'xD': 0.95, 'xB': 0.95}  # Model arguments
    m = build_column(**model_args)  # Building the column model

