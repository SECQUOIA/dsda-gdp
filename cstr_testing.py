import csv
import logging
import os
from math import ceil, fabs
import time

import pyomo.environ as pe
from pyomo.environ import SolverFactory, Suffix, value
from pyomo.gdp import Disjunct, Disjunction
from gdp.cstr.gdp_reactor import build_cstrs

if __name__ == "__main__":

    # Results
    NTs = range(5, 31, 1)
    # NTs = [5]
    time_limit = 900
    nlps = ['knitro', 'baron']
    starting_point = [1, 1]

    globaltee = True
   
    
    transformations = ['bigm', 'hull']
    ks = ['L2', 'Linf']
    strategies = ['LOA', 'GLOA', 'LBB']

    # Opening the CSV file here
    with open('cstr_ldsda.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        # Write headers
        writer.writerow(['NT', 'NLP_Solver', 'k', 'Objective', 'Time (s)'])
    
        # LD-SDA
        for NT in NTs:
            for k in ks:
                    for solver in nlps:
                        m = build_cstrs(NT)
                        new_results = {}
                        start_time = time.time()
                        result = pe.SolverFactory('gdpopt.ldsda').solve(
                            m,
                            minlp_solver='gams',
                            nlp_solver=solver,
                            tee=True,
                            starting_point=starting_point,
                            logical_constraint_list=[m.one_unreacted_feed,
                            m.one_recycle
                            ],
                            direction_norm=k,
                            time_limit=time_limit,
                        )

                        end_time = time.time()  # End timing
                        elapsed_time = end_time - start_time  # Calculate elapsed time
                        new_results = {
                                'NT': NT,
                                'Solver': solver, 
                                'k': str(k), 
                                'Objective': pe.value(m.obj), 
                                'Time': elapsed_time
                            }
                        writer.writerow([new_results['NT'], new_results['Solver'], new_results['k'], new_results['Objective'], new_results['Time']])
                        print(new_results)

    #     # D-SDA
    #     m = build_cstrs(NT)
    #     ext_ref = {m.YF: m.N, m.YR: m.N}
    #     get_external_information(m, ext_ref, tee=False)

    #     for solver in nlps:
    #         for k in ks:
    #             for transformation in ['hull','bigm']:
    #                 new_result = {}
    #                 m_solved, _, _ = solve_with_dsda(
    #                     model_function=build_cstrs,
    #                     model_args={'NT': NT},
    #                     starting_point=starting_point,
    #                     ext_dict=ext_ref,
    #                     ext_logic=problem_logic_cstr,
    #                     mip_transformation=True,
    #                     transformation=transformation,
    #                     k=k,
    #                     provide_starting_initialization=True,
    #                     feasible_model='cstr_' + str(NT),
    #                     subproblem_solver=solver,
    #                     subproblem_solver_options=nlp_opts[solver],
    #                     iter_timelimit=timelimit,
    #                     timelimit=timelimit,
    #                     gams_output=False,
    #                     tee=False,
    #                     global_tee=False,
    #                 )
    #                 new_result = {'Method': str('D-SDA_MIP_'+transformation), 'Approach': str('k='+k), 'Solver': solver, 'Objective': pe.value(
    #                     m_solved.obj), 'Time': m_solved.dsda_time, 'Status': m_solved.dsda_status, 'User_time': m_solved.dsda_usertime, 'NT': NT}
    #                 dict_data.append(new_result)
    #                 print(new_result)


