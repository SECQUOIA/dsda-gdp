from __future__ import division

# Import various modules and functions needed for the script
import csv  # To handle CSV files
import logging  # To keep logs for tracking
import os  # To access the OS functionalities for file and directory handling
from math import ceil, fabs  # Importing ceil and fabs functions from math
import time  # To keep track of time taken to solve the optimization problem
import pyomo.environ as pe  # To create and solve optimization models

# Importing specific classes and functions from pyomo.environ
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
    Suffix,
    TransformationFactory,
    Var,
    exactly,
    land,
    log,
    lor,
    minimize,
    value,
)

# Importing Disjunct and Disjunction classes from pyomo.gdp for creating generalized disjunctive programming models
from pyomo.gdp import Disjunct, Disjunction

# Importing utility function to log infeasible constraints
from pyomo.util.infeasible import log_infeasible_constraints

# Importing build_column function from gdp.column.gdp_column
from gdp.column.gdp_column import build_column

if __name__ == "__main__":
    # This part of the script initializes the variables that are going to be used throughout the script.
    # NT: Number of trays in the distillation column
    # timelimit: time limit for solving the optimization problem
    # model_args: dictionary containing arguments for the distillation column model
    # starting_point: List containing initial guesses for the number of trays and reflux.
    # globaltee: Boolean indicating whether to display solver output.
    # logging: Configuring the logging level to ERROR. This will avoid printing out warning messages.

    NT = 17
    time_limit = 900  # [s]
    model_args = {'min_trays': 8, 'max_trays': NT, 'xD': 0.95, 'xB': 0.95}
    starting_point = [NT - 2, 1]
    globaltee = True
    logging.basicConfig(level=logging.ERROR)

    nlps = ['knitro', 'baron']
    ks = ['Linf', 'L2']

    m = build_column(**model_args)

    results = []

    for solver in nlps:
        for k in ks:
            new_results={}
            start_time = time.time()
            result = pe.SolverFactory('gdpopt.ldsda').solve(
                m,
                minlp_solver='gams',
                nlp_solver=solver,
                tee=True,
                starting_point=starting_point,
                logical_constraint_list=[m.one_reflux,
                m.one_boilup
                ],
                direction_norm=k,
                time_limit=time_limit,
            )
            end_time = time.time()
            elapsed_time = end_time - start_time
            new_results={"Solver": solver, "k": str(k), "Objective": pe.value(m.obj), "Time": elapsed_time}
            results.append(new_results)
            print(new_results)

        # Write the results to a CSV file
        csv_file_path = 'column_testing.csv'
        with open(csv_file_path, mode='w', newline='') as csv_file:
            fieldnames = ['Solver', 'k', 'Objective', 'Time']
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

            writer.writeheader()
            for result in results:
                writer.writerow(result)

        print(f"Results have been written to {csv_file_path}")