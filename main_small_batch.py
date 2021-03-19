from math import ceil, fabs
import pyomo.environ as pe
from pyomo.gdp import Disjunct, Disjunction
from pyomo.util.infeasible import log_infeasible_constraints

from gdp.small_batch.gdp_small_batch import build_small_batch, external_ref
from gdp.dsda.dsda_functions import (
    generate_initialization, initialize_model, solve_subproblem, solve_with_dsda,
    solve_with_gdpopt, solve_with_minlp, visualize_dsda)


if __name__ == "__main__":
    # Inputs
    timelimit = 10
    model_args = {}

    # MINLP solution
    m = build_small_batch()
    m_init = initialize_model(m, from_feasible=True, feasible_model='small_batch')
    m_solved = solve_with_minlp(m_init, transformation='bigm', minlp='baron', timelimit=timelimit, gams_output=False)
    print(m_solved.results)

    # GDPopt method
    m = build_small_batch()
    m_init = initialize_model(m, from_feasible=True, feasible_model='small_batch')
    m_solved = solve_with_gdpopt(m_init, mip='cplex', nlp='conopt',
                                 timelimit=timelimit, strategy='LOA', mip_output=False, nlp_output=False)
    print(m_solved.results)

    # # D-SDA
    k = 'Infinity'
    starting_point = [3, 3, 3]
    min_allowed = {i: 1 for i in range(1, len(starting_point)+1)}
    max_allowed = {i: 3 for i in range(1, len(starting_point)+1)}

    m_solved, route = solve_with_dsda(model_function=build_small_batch, model_args={}, starting_point=starting_point, reformulation_function=external_ref, k=k,
                                      provide_starting_initialization=True, feasible_model='small_batch', subproblem_solver='msnlp', min_allowed=min_allowed, max_allowed=max_allowed, iter_timelimit=10, timelimit=3600)
    print(route)
    print(m_solved.results)
