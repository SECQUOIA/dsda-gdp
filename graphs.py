import csv
import json
import os
import matplotlib
from math import ceil, fabs, log10

import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pyomo.environ as pe
from pyomo.environ import SolverFactory, Suffix, value
from pyomo.gdp import Disjunct, Disjunction
from pyomo.util.infeasible import log_infeasible_constraints

from gdp.cstr.gdp_reactor import build_cstrs
from gdp.column.gdp_column import build_column
from gdp.small_batch.gdp_small_batch import build_small_batch
from gdp.dsda.dsda_functions import solve_with_dsda, visualize_dsda, get_external_information
from main_cstr import problem_logic_cstr
from main_column import problem_logic_column


if __name__ == "__main__":

    # # Column Graphs
    # ks = ['Infinity','2']
    # NT = 17
    # starting_point = [NT-2,1]
    # model_args = {'min_trays': 8, 'max_trays': NT, 'xD': 0.95, 'xB': 0.95}
    # routes = []

    # for k in ks:
    #     column = pd.read_csv("compl_enum_column_17_knitro_.csv") 
    #     column.head()
    #     column["Scaled_Objective"] = ( (column["Objective"])/ (19346) )

    #     m = build_column(**model_args)
    #     ext_ref = {m.YB: m.intTrays, m.YR: m.intTrays}
    #     get_external_information(m, ext_ref, tee=False)

    #     _, route, _ = solve_with_dsda(
    #         model_function=build_column,
    #         model_args=model_args,
    #         starting_point=starting_point,
    #         ext_dict=ext_ref,
    #         ext_logic=problem_logic_column,
    #         k=k,
    #         provide_starting_initialization=True,
    #         feasible_model='column_' + str(NT),
    #         subproblem_solver='knitro',
    #         subproblem_solver_options={},
    #         iter_timelimit=15,
    #         timelimit=3600,
    #         gams_output=False,
    #         tee=False,
    #         global_tee=False,
    #     )
    #     routes.append(route)

  

    # def drawArrow(A, B, color, shape):
    #     plt.arrow(A[0], A[1], B[0] - A[0], B[1] - A[1], width=0.00005,
    #             head_width=0.15, head_length=0.1, color=color, shape=shape)

    # shape = 'full'
    # color = ''
    # for j in range(len(routes)):
    #     if j == 0:
    #         color = 'black'
    #     else:
    #         color = 'orangered'
    #     for i in range(len(routes[j])-1):
    #         drawArrow(routes[j][i], routes[j][i+1], color, shape)

    # x_inf, y_inf, x_feas, y_feas, obj_feas = [], [], [], [], []
    # for i in range(len(column)):
    #     if column['Status'][i]=='Evaluated_Infeasible':
    #         x_inf.append(column['x'][i])
    #         y_inf.append(column['y'][i])
    #     else:
    #         x_feas.append(column['x'][i])
    #         y_feas.append(column['y'][i])
    #         obj_feas.append(column['Objective'][i])

    # sizer = [i / 300 for i in obj_feas]
    # cm = plt.cm.get_cmap('viridis_r')
    # sc1 = plt.scatter(x_inf, y_inf, s=60, c='white', marker='o', edgecolors='black')
    # sc = plt.scatter(x_feas, y_feas, s=sizer, c=obj_feas, cmap=cm)
    # cbar = plt.colorbar(sc)
    # cbar.set_label('Objective', rotation=90)
    # title_string = 'D-SDA for the Column Superstructure NT=17'
    # plt.title(title_string)
    # plt.xlabel('Reflux Position XR')
    # plt.ylabel('Boil-up Position XB')
    # plt.show()

    # _______________________________________________________________________________

    # # CSTR D-SDA Graphs
    # ks = ['Infinity','2']
    # starting_point = [1,1]
    # NT = 25
    # routes = []
    # for k in ks:
    #     cstr = pd.read_csv("compl_enum_cstr_25_baron_.csv") 
    #     cstr.head()
    #     cstr["Scaled_Objective"] = ( (cstr["Objective"] - 2.71)/ (9.89 - 2.71) )

    #     m = build_cstrs(NT)
    #     ext_ref = {m.YF: m.N, m.YR: m.N}
    #     get_external_information(m, ext_ref, tee=False)

        
    #     _, route, _ = solve_with_dsda(
    #         model_function=build_cstrs,
    #         model_args={'NT': NT},
    #         starting_point=starting_point,
    #         ext_dict=ext_ref,
    #         ext_logic=problem_logic_cstr,
    #         k=k,
    #         provide_starting_initialization=True,
    #         feasible_model='cstr_' + str(NT),
    #         subproblem_solver='knitro',
    #         subproblem_solver_options={},
    #         iter_timelimit=15,
    #         timelimit=3600,
    #         gams_output=False,
    #         tee=False,
    #         global_tee=False,
    #     )
    #     routes.append(route)
    # # route = []
    # cm = plt.cm.get_cmap('viridis_r')

    # def drawArrow(A, B, color, shape):
    #     plt.arrow(A[0], A[1], B[0] - A[0], B[1] - A[1], width=0.00005,
    #             head_width=0.45, head_length=0.25, color=color, shape=shape)

    # shape = 'full'
    # color = ''
    # for j in range(len(routes)):
    #     if j == 0:
    #         color = 'black'
    #     else:
    #         color = 'orangered'
    #     for i in range(len(routes[j])-1):
    #         drawArrow(routes[j][i], routes[j][i+1], color, shape)

    # #drawArrow([5,15],[6,15],'black',shape)
    # #drawArrow([5,23],[6,23],'orangered',shape)

    # sc = plt.scatter(cstr.x, cstr.y, s=10*cstr.Objective, c=cstr.Scaled_Objective, cmap=cm, norm=matplotlib.colors.LogNorm())
    # cbar = plt.colorbar(sc)
    # cbar.set_label('Scaled Objective', rotation=90)
    # title_string = 'D-SDA for the CSTR Superstructure NT=25'
    # plt.title(title_string)
    # plt.xlabel('Number of Reactors XF')
    # #plt.xlabel('k = Infinity       k = 1')
    # plt.ylabel('Reflux Position XR')
    # plt.show()
