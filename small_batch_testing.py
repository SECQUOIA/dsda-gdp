import csv
import logging
import os
from math import ceil, fabs
import time

import pyomo.environ as pe
from pyomo.gdp import Disjunct, Disjunction
from pyomo.util.infeasible import log_infeasible_constraints

from gdp.small_batch.gdp_small_batch import build_small_batch

if __name__ == "__main__":
    # Inputs
    time_limit = 900
    model_args = {}
    starting_point = [3, 3, 3]

    globaltee = True

    nlps = ['knitro', 'baron']

    transformations = ['bigm', 'hull']
    ks = ['L2', 'Linf']

    m = build_small_batch()

    with open('small_batch_ldsda.csv', mode='w', newline='') as file:
        writer = csv.writer(file)
        # Write headers
        writer.writerow(['NLP_Solver', 'k', 'Objective', 'Time (s)'])

        for solver in nlps:
            for k in ks:
                new_results={}
                start_time = time.time()
                result = pe.SolverFactory('gdpopt.ldsda').solve(
                    m,
                    minlp_solver='gams',
                    nlp_solver=solver,
                    tee=True,
                    starting_point=[3, 3, 3],
                    logical_constraint_list=[
                        m.lim['mixer'].name,
                        m.lim['reactor'].name,
                        m.lim['centrifuge'].name,
                    ],
                    direction_norm=k,
                    time_limit=time_limit,
                )
                end_time = time.time()
                elapsed_time = end_time - start_time

                print(result)
                new_results={'Solver': solver, 'k': str(k), 'Objective': pe.value(m.obj), 'Time': elapsed_time}

                writer.writerow([new_results['Solver'], new_results['k'], new_results['Objective'], new_results['Time']])

                print(new_results)

                print("objective: ", pe.value(m.obj))
                for k in m.k:
                    for j in m.j:
                        print(f"Y[{k},{j}] = {pe.value(m.Y[k, j])}")
