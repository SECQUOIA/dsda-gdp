"""
Distillation column model for 2018 PSE conference formulated into GDP models.
References:
- Ghouse, Jaffer H., et al. "A comparative study between GDP and NLP formulations for conceptual design of distillation columns." Computer Aided Chemical Engineering. Vol. 44. Elsevier, 2018. 865-870.
"""
# The gdp_column.py formulates the build column model, state the energy and the mass balances for every part of the column.
# The model use the default value for the initial guess.
# The model does not have the iteration
# The model states the Boolean variable, the LD-SDA is done on main_column.py

import math  # Provides functions for mathematical operations.
import os  # Provides functions for interacting with the operating system.

# Imports from the Pyomo library for building and solving optimization problems.
from pyomo.common.errors import InfeasibleConstraintException
from pyomo.contrib.fbbt.fbbt import fbbt
from pyomo.contrib.gdpopt.data_class import MasterProblemResult
from pyomo.core.base.misc import display
from pyomo.core.plugins.transform.logical_to_linear import (
    update_boolean_vars_from_binary,
)  # Transforms logical constraints into binary constraints.
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
    land,
    log,
    lor,
    minimize,
    value,
)  # Core components of Pyomo's modeling language.
from pyomo.gdp import Disjunct, Disjunction  # Imports for disjunctive programming.
from pyomo.opt import (
    SolutionStatus,
    SolverResults,
)  # Classes to represent solver results.
from pyomo.opt import (
    TerminationCondition as tc,
)  # Enum of possible termination conditions for solvers.
from pyomo.opt.base.solvers import SolverFactory  # Base class for solver factories.


def build_column(min_trays, max_trays, xD, xB):
    """
    Builds the column model.
    References: Ghouse, Jaffer H., et al. "A comparative study between GDP and NLP formulations for conceptual design of distillation columns." Computer Aided Chemical Engineering. Vol. 44. Elsevier, 2018. 865-870.

    Args:
        min_trays (int): Minimum number of trays in the column
        max_trays (int): Maximum number of trays in the column
        xD (float): Distillate purity
        xB (float): Bottoms purity

    Returns:
        m (ConcreteModel): Pyomo model
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
    def tray_no_tray(b, t):
        """Disjunction statement defining whether a tray exists or not"""
        return [b.tray[t], b.no_tray[t]]

    # Constraint for minimum number of trays, adding 1 for feed tray
    m.minimum_num_trays = Constraint(
        expr=sum(m.tray[t].indicator_var for t in m.conditional_trays)
        + 1  # for feed tray
        >= min_trays
    )

    # Define variables
    m.T_feed = Var(
        doc='Feed temperature [K]',
        domain=NonNegativeReals,
        bounds=(min_T, max_T),
        initialize=368,
    )  # Feed temperature [K]
    m.feed_vap_frac = Var(
        doc='Vapor fraction of feed', initialize=0, bounds=(0, 1)
    )  # Vapor fraction of feed
    m.feed = Var(
        m.comps, doc='Total component feed flow [mol/s]', initialize=50
    )  # Total component feed flow [mol/s]
    m.x = Var(
        m.comps,
        m.trays,
        doc='Liquid mole fraction',
        bounds=(0, 1),
        domain=NonNegativeReals,
        initialize=0.5,
    )  # Liquid mole fraction
    m.y = Var(
        m.comps,
        m.trays,
        doc='Vapor mole fraction',
        bounds=(0, 1),
        domain=NonNegativeReals,
        initialize=0.5,
    )  # Vapor mole fraction
    m.L = Var(
        m.comps,
        m.trays,
        doc='component liquid flows from tray [mol/s]',
        domain=NonNegativeReals,
        bounds=(0, max_flow),
        initialize=50,
    )  # Component liquid flows from tray [mol/s]
    m.V = Var(
        m.comps,
        m.trays,
        doc='component vapor flows from tray [mol/s]',
        domain=NonNegativeReals,
        bounds=(0, max_flow),
        initialize=50,
    )  # Component vapor flows from tray [mol/s]
    m.liq = Var(
        m.trays,
        domain=NonNegativeReals,
        doc='liquid flows from tray [mol/s]',
        initialize=100,
        bounds=(0, max_flow),
    )  # Liquid flows from tray [mol/s]
    m.vap = Var(
        m.trays,
        domain=NonNegativeReals,
        doc='vapor flows from tray [mol/s]',
        initialize=100,
        bounds=(0, max_flow),
    )  # Vapor flows from tray [mol/s]
    m.B = Var(
        m.comps,
        domain=NonNegativeReals,
        doc='bottoms component flows in mol/s',
        bounds=(0, max_flow),
        initialize=50,
    )  # Bottoms component flows [mol/s]
    m.D = Var(
        m.comps,
        domain=NonNegativeReals,
        doc='distillate component flows [mol/s]',
        bounds=(0, max_flow),
        initialize=50,
    )  # Distillate component flows [mol/s]
    m.bot = Var(
        domain=NonNegativeReals,
        initialize=50,
        bounds=(0, 100),
        doc='bottoms flow [mol/s]',
    )  # Bottoms flow [mol/s]
    m.dis = Var(
        domain=NonNegativeReals,
        initialize=50,
        doc='distillate flow [mol/s]',
        bounds=(0, 100),
    )  # Distillate flow [mol/s]
    m.reflux_ratio = Var(
        domain=NonNegativeReals, bounds=(0.5, 4), doc='reflux ratio', initialize=1.4
    )  # Reflux ratio
    m.reboil_ratio = Var(
        domain=NonNegativeReals, bounds=(1.3, 4), doc='reboil ratio', initialize=0.9527
    )  # Reboil ratio
    m.reflux_frac = Var(
        domain=NonNegativeReals, bounds=(0, 1 - 1e-6), doc='reflux fractions'
    )  # Reflux fractions
    m.boilup_frac = Var(
        domain=NonNegativeReals, bounds=(0, 1 - 1e-6), doc='boilup fraction'
    )  # Boilup fraction
    m.Kc = Var(
        m.comps,
        m.trays,
        doc='Phase equilibrium constant',
        domain=NonNegativeReals,
        initialize=1,
        bounds=(0, 1000),
    )  # Phase equilibrium constant
    m.T = Var(
        m.trays, doc='Temperature [K]', domain=NonNegativeReals, bounds=(min_T, max_T)
    )  # Tray temperature [K]
    m.P = Var(doc='Pressure [bar]', bounds=(0, 5))  # Pressure [bar]
    m.gamma = Var(
        m.comps,
        m.trays,
        doc='liquid activity coefficent of component on tray',
        domain=NonNegativeReals,
        bounds=(0, 10),
        initialize=1,
    )  # Liquid activity coefficient
    m.Pvap = Var(
        m.comps,
        m.trays,
        doc='pure component vapor pressure of component on tray [bar]',
        domain=NonNegativeReals,
        bounds=(1e-3, 5),
        initialize=0.4,
    )  # Pure component vapor pressure [bar]
    m.Pvap_rel = Var(
        m.comps,
        m.trays,
        doc='pure component relative vapor pressure of component on tray in bar (to avoid numerical problems)',
        domain=NonNegativeReals,
        bounds=(0, 5),
        initialize=0.4,
    )  # Pure component relative vapor pressure [bar]
    m.Pvap_X = Var(
        m.comps,
        m.trays,
        doc='Related to fraction of critical temperature (1 - T/Tc)',
        bounds=(0.25, 0.5),
        initialize=0.4,
    )  # Related to fraction of critical temperature [kJ/mol]
    m.H_L = Var(
        m.comps,
        m.trays,
        bounds=(0.1, 16),
        doc='Liquid molar enthalpy of component in tray [kJ/mol]',
    )  # Liquid molar enthalpy of component in tray [kJ/mol]
    m.H_V = Var(
        m.comps,
        m.trays,
        bounds=(30, 16 + 40),
        doc='Vapor molar enthalpy of component in tray [kJ/mol]',
    )  # Vapor molar enthalpy of component in tray [kJ/mol]
    m.H_L_spec_feed = Var(
        m.comps,
        doc='Component liquid molar enthalpy in feed [kJ/mol]',
        initialize=0,
        bounds=(0.1, 16),
    )  # Component liquid molar enthalpy in feed [kJ/mol]
    m.H_V_spec_feed = Var(
        m.comps,
        doc='Component vapor molar enthalpy in feed [kJ/mol]',
        initialize=0,
        bounds=(30, 16 + 40),
    )  # Component vapor molar enthalpy in feed [kJ/mol]
    m.Qb = Var(
        domain=NonNegativeReals, doc='reboiler duty [MJ/s]', initialize=1, bounds=(0, 8)
    )  # Reboiler duty [MJ/s]
    m.Qc = Var(
        domain=NonNegativeReals,
        doc='condenser duty [MJ/s]',
        initialize=1,
        bounds=(0, 8),
    )  # Condenser duty [MJ/s]

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

    # Define the objective function as the sum of reboiler and condenser duty plus an indicator for tray activation
    # The objective is to minimize the sum of condenser and reboiler duties, Qc and Qb, multiplied by 1E3 to convert units,
    # and also the number of activated trays, which is obtained by summing up the indicator variables for the trays by 1E3 [$/No. of Trays].
    m.obj = Objective(
        expr=(m.Qc + m.Qb) * 1e3
        + (sum(m.tray[t].indicator_var for t in m.conditional_trays) + 1) * 1e3,
        sense=minimize,
    )

    # Constraint to calculate the reflux ratio
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

    # Defining set of interior trays in the column (excluding condenser and reboiler trays)
    m.intTrays = Set(
        initialize=m.trays - [m.condens_tray, m.reboil_tray],
        doc='Interior trays of the column',
    )

    # Defining Boolean variables to denote existence of boil-up flow and reflux flow in each interior tray
    m.YB = BooleanVar(
        m.intTrays, initialize=False, doc='Existence of boil-up flow in stage n'
    )
    m.YR = BooleanVar(
        m.intTrays, initialize=False, doc='Existence of reflux flow in stage n'
    )

    # Initializing at least one reflux and boilup tray to avoid errors in Mixed-Integer Nonlinear Programming (MINLP) solvers
    m.YB[m.reboil_tray + 1].set_value(True)
    m.YR[m.max_trays - 1].set_value(True)

    # Defining additional Boolean variables for logical constraints
    m.YP = BooleanVar(m.intTrays, doc='Boolean var associated with tray and no_tray')
    m.YB_is_up = BooleanVar(
        doc='Boolean var for intermediate sum determining if Boilup is above the feed'
    )
    m.YR_is_down = BooleanVar(
        doc='Boolean var for intermediate sum determining if Reflux is below the feed'
    )

    # Defining logical constraints to ensure only one reflux and one boilup
    @m.LogicalConstraint()
    def one_reflux(m):
        """Ensure that only one reflux flow exists."""
        return exactly(1, m.YR)

    @m.LogicalConstraint()
    def one_boilup(m):
        """Ensure that only one boil-up flow exists."""
        return exactly(1, m.YB)

    # Defining logical constraint for YP Boolean variable
    @m.LogicalConstraint(m.conditional_trays)
    def YP_or_notYP(m, n):
        return m.YP[n].equivalent_to(
            land(
                lor(m.YR[j] for j in range(n, m.max_trays)),
                lor(land(~m.YB[j] for j in range(n, m.max_trays)), m.YB[n]),
            )
        )

    # Associating YP Boolean variable with tray activation
    for n in m.conditional_trays:
        m.YP[n].associate_binary_var(m.tray[n].indicator_var)

    # Fixing feed conditions
    m.feed['benzene'].fix(50)  # Fixing benzene flow in the feed at 50 [mol/s]
    m.feed['toluene'].fix(50)  # Fixing toluene flow in the feed at 50 [mol/s]
    m.T_feed.fix(368)  # Fixing feed temperature at 368 [K]
    m.feed_vap_frac.fix(0.40395)  # Fixing feed vapor fraction
    m.P.fix(1.01)  # Fixing pressure at 1.01 [bar]

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
    Builds the constraints for mass balance, liquid and vapor composition for a given tray (t) in the distillation column.
    The constraints model the behavior of the mass balance for different components on a tray, accounting for feed, vapor, and liquid flows,
    as well as special conditions for the feed tray, condenser, and reboiler. Additional constraints define the liquid and vapor composition
    on the tray, as well as conditions for when the tray does not exist.

    Args:
        m: The model object containing the relevant variables and parameters.
        t: Tray number for which the constraints are being defined (integer).
        tray: Disjunct object representing the case when the tray exists in the column.
        no_tray: Disjunct object representing the case when the tray is absent in the column.

    Return:
        None. The function adds constraints to the model but does not return a value.
    """

    @tray.Constraint(m.comps)
    def mass_balance(_, c):
        """Mass balance on each component on a tray."""
        return (
            m.feed[c] if t == m.feed_tray else 0
        ) - m.V[  # Total feed to the tray if it's a feed tray
            c, t
        ] - (  # Total vapor flow from the tray
            m.D[c] if t == m.condens_tray else 0
        ) + (  # Total distillate if it's a condenser tray
            m.L[c, t + 1] if t < m.condens_tray else 0
        ) - (  # Liquid flow from the tray above if it's not a condenser
            m.B[c] if t == m.reboil_tray else 0
        ) - (  # Total bottoms if it's a reboiler tray
            m.L[c, t] if t > m.reboil_tray else 0
        ) + (  # Liquid flow to the tray below if it's not a reboiler
            m.V[c, t - 1] if t > m.reboil_tray else 0
        ) == 0  # Vapor flow from the tray below if it's not a reboiler

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
        m: The model object containing variables and parameters related to the feed tray, components, and mass streams.

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
        m (Model Object): A model object that includes the relevant variables, parameters,
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
        m (Model Object): A model object that includes the relevant variables, parameters,
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
    Constructs the phase equilibrium constraints for a given tray in the distillation column.
    This function models the equilibrium relationships between the vapor and liquid phases on a tray,
    based on Raoult's law, the phase equilibrium constant, the relative vapor pressure, and the
    Antoine equation.

    Args:
        m: The model object containing the relevant variables and parameters such as vapor and liquid compositions,
           phase equilibrium constants, relative vapor pressure, and temperature-dependent factors.
        t: The specific tray number for which the phase equilibrium is being modeled.
        tray: A container object within the model representing the tray for which the constraints are being built.

    Constraints:
        - raoults_law: Models the relationship between vapor and liquid composition for each component on the tray
          using Raoult's law.
        - phase_equil_const: Defines the phase equilibrium constant for each component on the tray.
        - Pvap_relative: Defines the relative vapor pressure for each component on the tray.
        - Pvap_relation: Establishes the relationship between the relative vapor pressure and temperature for each
          component on the tray using the Antoine equation.
        - Pvap_X_defn: Defines the temperature-dependent part of the relative vapor pressure for each component on the
          tray.
        - gamma_calc: Assumes an ideal solution (activity coefficient, gamma = 1) for each component on the tray.

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

    # Definition of the relative vapor pressure for each component on the tray
    @tray.Constraint(m.comps)
    def Pvap_relative(_, c):
        return m.Pvap_rel[c, t] == m.Pvap[c, t] - m.Pvap[c, t].lb

    # Relation between the relative vapor pressure and temperature for each component on the tray
    @tray.Constraint(m.comps)
    def Pvap_relation(_, c):
        """The equation uses the Antoine equation to relate the vapor pressure to one minus the reduced temperature (represented by the variable x)"""
        k = m.pvap_const[c]
        x = m.Pvap_X[c, t]
        return (log(m.Pvap_rel[c, t] + m.Pvap[c, t].lb) - log(k['Pc'])) * (1 - x) == (
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
        m (Model Object): A model object that includes relevant variables, parameters, and expressions for
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


if __name__ == "__main__":
    # Inputs
    NT = 17  # Total number of trays
    model_args = {
        'min_trays': 8,
        'max_trays': NT,
        'xD': 0.95,
        'xB': 0.95,
    }  # Model arguments
    m = build_column(**model_args)  # Building the column model
