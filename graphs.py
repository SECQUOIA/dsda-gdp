import csv
import json
import os
import matplotlib
import ast
from math import ceil, fabs, log10

import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pyomo.environ as pe
import mpl_toolkits.mplot3d
from pyomo.environ import SolverFactory, Suffix, value
from pyomo.gdp import Disjunct, Disjunction
from pyomo.util.infeasible import log_infeasible_constraints
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.axes_grid1 import Divider, Size
from matplotlib.ticker import MaxNLocator
from gdp.cstr.gdp_reactor import build_cstrs
from gdp.column.gdp_column import build_column
from gdp.small_batch.gdp_small_batch import build_small_batch
from gdp.dsda.dsda_functions import solve_with_dsda, visualize_dsda, get_external_information
from main_cstr import problem_logic_cstr
from main_column import problem_logic_column

class Arrow3D(FancyArrowPatch):
    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0,0), (0,0), *args, **kwargs)
        self._xyz = (x,y,z)
        self._dxdydz = (dx,dy,dz)

    def draw(self, renderer):
        x1,y1,z1 = self._xyz
        dx,dy,dz = self._dxdydz
        x2,y2,z2 = (x1+dx,y1+dy,z1+dz)

        xs, ys, zs = proj_transform((x1,x2),(y1,y2),(z1,z2), renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        super().draw(renderer)

    def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
        '''Add an 3d arrow to an `Axes3D` instance.'''

        arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
        ax.add_artist(arrow)

    setattr(Axes3D,'arrow3D',_arrow3D)

def set_axes_equal(ax: plt.Axes):
    """Set 3D plot axes to equal scale.

    Make axes of 3D plot have equal scale so that spheres appear as
    spheres and cubes as cubes.  Required since `ax.axis('equal')`
    and `ax.set_aspect('equal')` don't work on 3D.
    """
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    _set_axes_radius(ax, origin, radius)

def _set_axes_radius(ax, origin, radius):
    x, y, z = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    ax.set_zlim3d([z - radius, z + radius])


if __name__ == "__main__":

    # # Small Batch Graph
    batch = pd.read_csv("results/compl_enum_small_batch_baron.csv") 
    # batch.head()


    # fig = plt.figure(figsize = (9, 7))
    # cm = plt.cm.get_cmap('viridis_r')
    # ax = plt.axes(projection ="3d")

    # sctt = ax.scatter3D(batch.x, batch.y, batch.z, s=120, c = batch['Objective'], cmap=cm, alpha=1)

    # # k = Infinity
    # ax.arrow3D(3,3,3,-1,-1,-1,mutation_scale=25,arrowstyle="-|>",ec ='black',fc='black')
    # ax.arrow3D(2,2,2,0,0,-1,mutation_scale=25,arrowstyle="-|>",ec ='black',fc='black')

    # # k = 2
    # ax.arrow3D(3,3,3,0,0,-1,mutation_scale=25,arrowstyle="-|>",ec ='orangered',fc='orangered')
    # ax.arrow3D(3,3,2,0,0,-1,mutation_scale=25,arrowstyle="-|>",ec ='orangered',fc='orangered')
    # ax.arrow3D(3,3,1,-1,0,0,mutation_scale=25,arrowstyle="-|>",ec ='orangered',fc='orangered')
    # ax.arrow3D(2,3,1,0,-1,0,mutation_scale=25,arrowstyle="-|>",ec ='orangered',fc='orangered')

    # cbar = plt.colorbar(sctt, shrink=0.8, aspect = 20)
    # cbar.set_label('Objective', rotation=90)

    # ax.invert_xaxis()

    # ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    # ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    # ax.zaxis.set_major_locator(MaxNLocator(integer=True))

    # ax.grid(b = True, color ='grey',
    #     linestyle ='-.', linewidth = 0.3,
    #     alpha = 0.2)

    # ax.set_box_aspect([1,1,2]) 
    # ax.set_xlabel('Number of Mixers '+ r'$X_{mixer}$' , fontsize=10)
    # ax.set_ylabel('Number of Reactors '+ r'$X_{reactor}$', fontsize=10)
    # ax.set_zlabel('Number of Centrifuges '+ r'$X_{centrifuge}$', fontsize=10)

    # plt.show()

    # # _______________________________________________________________________________
    # # CSTR Gap vs NT Graph
    # cstr = pd.read_csv("results/cstr_results.csv") 
    # cstr.head()
    # enumb = pd.read_csv("results/compl_enum_cstr_25_baron.csv")
    # enumk = pd.read_csv("results/compl_enum_cstr_25_knitro.csv")
    # enumm = pd.read_csv("results/compl_enum_cstr_25_msnlp.csv")

    
    # colors={'bigm':'blue','hull':'turquoise','LOA':'purple','GLOA':'red','LBB':'green','k=2':'slategrey','k=Infinity':'black','total':'gold'}
    # shapes={'antigone':'2','baron':'o','scip':'s','dicopt':'*','sbb':'+','msnlp':'D','knitro':'^'}

    # absolutes = {}
    # sumb, sumk, summ = {}, {}, {}
    # objb, objk, objm = {}, {}, {}

    # for i in range(len(enumb)):
    #     if enumb['x'][i] == enumb['y'][i]:
    #         absolutes[enumb['x'][i]] = enumb['Objective'][i]
    #         sumb[enumb['x'][i]] = enumb['Global_Time'][i]
    #         objb[enumb['x'][i]] = enumb['Objective'][i]

    # for i in range(len(enumk)):
    #     if enumk['x'][i] == enumk['y'][i]:
    #         sumk[enumk['x'][i]] = enumk['Global_Time'][i]
    #         objk[enumk['x'][i]] = enumk['Objective'][i]

    # for i in range(len(enumm)):
    #     if enumm['x'][i] == enumm['y'][i]:
    #         summ[enumm['x'][i]] = enumm['Global_Time'][i]
    #         objm[enumm['x'][i]] = enumm['Objective'][i]  

    # print(cstr)

    # new_row = {}
    # for i in sumb.keys():
    #     new_row = {'Method':'Enumeration', 'Approach':'total', 'Solver':'baron', 'Objective':objb[i],'Time':sumb[i],'Status':'NA','NT':i}
    #     cstr = cstr.append(new_row, ignore_index=True)

    # new_row = {}
    # for i in sumk.keys():
    #     new_row = {'Method':'Enumeration', 'Approach':'total', 'Solver':'knitro', 'Objective':objk[i],'Time':sumk[i],'Status':'NA','NT':i}
    #     cstr = cstr.append(new_row, ignore_index=True)

    # new_row = {}
    # for i in summ.keys():
    #     new_row = {'Method':'Enumeration', 'Approach':'total', 'Solver':'msnlp', 'Objective':objm[i],'Time':summ[i],'Status':'NA','NT':i}
    #     cstr = cstr.append(new_row, ignore_index=True)


    # cstr['Gap'] = 101

    # cstr['Method'] = pd.Categorical(cstr['Method'], ["GDPopt", "MINLP", "Enumeration","D-SDA"])
    # cstr['Solver'] = pd.Categorical(cstr['Solver'], ["baron", "antigone", "msnlp", "sbb",'knitro','dicopt','scip'])
    # cstr['Approach'] = pd.Categorical(cstr['Approach'], ["bigm", "hull", "LOA","GLOA",'LBB','total','k=2','k=Infinity'])
    # #cstr = cstr.sort_values(by='Approach')
    # cstr = cstr.sort_values(by=['Method','Approach', 'Solver'])
    # cstr = cstr.reset_index(drop=True)

    # print(cstr)

    # artists = []
    # legends = []
    # max_gaps = [10, 0]
    # for max_gap in max_gaps:
    #     fig = plt.figure()
    #     ax = plt.gca()
    #     for i in range(len(cstr)):
    #         for j in absolutes.keys():
    #             if cstr['NT'][i] == j:
    #                 cstr['Gap'][i] = float(100*(cstr['Objective'][i] - absolutes[j])/absolutes[j])

    #         if cstr['Time'][i] < 1000 and cstr['NT'][i] >= 5:
    #             if cstr['Gap'][i] <= max_gap + 1e-3:
    #                 art = plt.scatter(cstr['NT'][i], cstr['Time'][i], marker=shapes[cstr['Solver'][i]],  c=colors[cstr['Approach'][i]])
    #                 leg = str(cstr['Approach'][i]+'-'+cstr['Solver'][i])
    #                 artists.append(art)
    #                 legends.append(leg)
    #             else:
    #                 art = plt.scatter([10], 10, marker='.',  c='white')
    #                 leg = "" #str(cstr['Approach'][i]+'-'+cstr['Solver'][i])
    #                 artists.append(art)
    #                 legends.append(leg)

  
    #     artists2, legends2 = [], []

    #     for i in range(len(legends)):
    #         if legends[i] not in legends2:
    #             legends2.append(legends[i])
    #             artists2.append(artists[i])

    #     order_legends, order_artists = [], []
    #     order = [22,25,19,5,1,3,8,11,13,17,23,26,20,6,2,4,0,9,14,15,24,27,21,7,0,0,12,10,18,16]
    #     for o in order:
    #         order_legends.append(legends2[o])
    #         order_artists.append(artists2[o])


    #     ax.legend(order_artists, order_legends, loc=8, prop={'size': 8}, scatterpoints=1, 
    #     ncol=3, framealpha=1, fancybox=False, edgecolor='black', bbox_to_anchor=(0.5, -0.75))
    #     ax.set_yscale('log')
    #     #fig.subplots_adjust(bottom=0.01)
    #     xs = [5,25]
    #     xint = range(min(xs), ceil(max(xs))+1,5)
    #     gap_str = ''
    #     if max_gap == 0:
    #         gap_str = 'Gap = 0%)'
    #     else:
    #         gap_str = 'Gap < 10%)'
    #     title_string = 'CSTR Superstructure Method Comparison (Achieved '+gap_str
    #     plt.title(title_string)
    #     plt.xlabel('Superstructure Size NT')
    #     plt.ylabel('Execution Time [s]')
    #     matplotlib.pyplot.xticks(xint)
    #     plt.show()




    # # _______________________________________________________________________________
    # # Column D-SDA Graphs
    # ks = ['Infinity','2']
    # NT = 17
    # starting_point = [NT-2,1]
    # model_args = {'min_trays': 8, 'max_trays': NT, 'xD': 0.95, 'xB': 0.95}
    # routes = []

    # for k in ks:
    #     column = pd.read_csv("results/compl_enum_column_17_knitro.csv") 
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
    # x_ns, y_ns = [], []
    # for i in range(len(column)):
    #     if column['Status'][i]=='Evaluated_Infeasible':
    #         x_inf.append(column['x'][i])
    #         y_inf.append(column['y'][i])
    #     elif column['Status'][i]=='No_Sol_Found':
    #         x_ns.append(column['x'][i])
    #         y_ns.append(column['y'][i])
    #     else:
    #         x_feas.append(column['x'][i])
    #         y_feas.append(column['y'][i])
    #         obj_feas.append(column['Objective'][i])

    # sizer = [i / 300 for i in obj_feas]
    # cm = plt.cm.get_cmap('viridis_r')
    # sc0 = plt.scatter(x_ns, y_ns, s=70, c='white', marker='o', edgecolors='black')
    # sc1 = plt.scatter(x_inf, y_inf, s=70, c='white', marker='^', edgecolors='black')
    # sc = plt.scatter(x_feas, y_feas, s=sizer, c=obj_feas, cmap=cm)
    # cbar = plt.colorbar(sc)
    # cbar.set_label('Objective', rotation=90)
    # title_string = 'D-SDA for the Column Superstructure NT=17'
    # plt.title(title_string)
    # #plt.xlabel('Reflux Position XR')
    # plt.xlabel('k = 2')
    # plt.ylabel('Boil-up Position XB')
    # plt.show()

    # # _______________________________________________________________________________

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