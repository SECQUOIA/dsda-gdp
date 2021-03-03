import pyomo.environ as pe
from pyomo.gdp import (Disjunct, Disjunction)
import numpy as np
import time as time
import itertools as it
import matplotlib.pyplot as plt
from pyomo.core.base.misc import display
from pyomo.opt.base.solvers import SolverFactory
from pyomo.common.errors import InfeasibleConstraintException
from pyomo.contrib.fbbt.fbbt import fbbt
from pyomo.contrib.gdpopt.data_class import MasterProblemResult
from pyomo.opt import SolutionStatus
from pyomo.opt import TerminationCondition as tc, SolverResults
import os

from column import build_column

def list_generator(NT):
    X1, X2, aux, aux2, x = [], [], [], 1, {}

    for i in range(1, NT+1):
        X1.append(i)
        aux.append(i)
        X2.append(aux2)

    for i in range(NT-1):
        aux.pop(0)
        aux2 += 1
        for j in aux:
            X1.append(j)
            X2.append(aux2)
    return X1, X2

def complete_enumeration(NT):
    X1, X2 = list_generator(NT)

    print('=============================')
    print('%6s %6s %12s' % ('x1', 'x2', 'Objective'))
    print('-----------------------------')

    for i in range(len(X1)):
        x = [X1[i], X2[i]]
        m, _, _ = build_column(min_trays=8,max_trays=NT,xD=0.95,xB=0.95,x_input=x,provide_init=False,init={})
        print('%6s %6s %12s' % (X1[i], X2[i], round(pe.value(m.obj), 5)))
    print('=============================')


def visualization(NT, points):
    X1, X2 = list_generator(NT)

    def drawArrow(A, B):
        plt.arrow(A[0], A[1], B[0] - A[0], B[1] - A[1], width=0.00005,
                  head_width=0.15, head_length=0.05, color='black', shape='full')

    for i in range(len(points)-1):
        drawArrow(points[i], points[i+1])

    plt.scatter(X1, X2, color='gray', marker='o')
    plt.show()

# Creates all posible directions with k=2 for num_ext variables

# num_ext is type int and means number of external variables
# directions is type dictionary starting from 1 with all posible directions


def neighborhood_k_eq_2(num_ext):
    num_neigh = 2*num_ext
    neighbors = np.concatenate((np.eye(num_ext), -np.eye(num_ext)), axis=1)
    directions = {}
    for i in range(num_neigh):
        direct = []
        directions[i+1] = direct
        for j in range(num_ext):
            direct.append(neighbors[j, i])
    return directions

def neighborhood_k_eq_inf(num_ext): # number
    num_neigh = 3*num_ext-1
    neighbors = list(it.product([-1,0,1], repeat=num_ext))
    directions = {}
    for i in range(len(neighbors)):
        directions[i+1]=list(neighbors[i])
    temp = directions.copy()
    for i in directions.keys():
        if temp[i] == [0]*num_ext:
            temp.pop(i,None)
    return temp
# Creates neighbor of a given point
# Optimize option will discard out of bounds points given by min and max allowed
# Cheating will discard infeasible points for CSTR problem (X2 - X1 > 0)
# Cheating will only be used while solving GAMS issue

# start is type list and stands for the actual point
# neighborhood is type dict and is the output of a k-Neighborhood function
# newbors or new_newbors is type  dict starting in 0 with neighbor of a given point
# Neighbor 0 is the actual point
def my_neighbors(start, neighborhood, optimize=True, min_allowed={}, max_allowed={}, cheating=False):
    neighbors = {0: start}
    for i in neighborhood.keys():
        neighbors[i] = list(map(sum, zip(start, list(neighborhood[i]))))

    if optimize:
        new_neighbors = {}
        num_vars = len(neighbors[0])
        for i in neighbors.keys():
            checked = 0
            for j in range(num_vars):
                if neighbors[i][j] >= min_allowed[j+1] and neighbors[i][j] <= max_allowed[j+1]:
                    checked += 1
            if checked == num_vars:
                new_neighbors[i] = neighbors[i]

        if cheating:
            temp = new_neighbors
            for i in list(temp.keys()):
                if new_neighbors[i][1] - new_neighbors[i][0] > 0:
                    new_neighbors.pop(i, None)

        return new_neighbors
    return neighbors


# Evaluates a group of given points and returns the best
# ext_vars is dict with given points where 0 is actual point
# init is type dict and contains solved variables for the actual point
# fmin is type int and stands for objective at actual point
# tol is type int and stands for numerical tolereance for equallity
# fmin as return is type int and gives the best neighbor's objective
# best_var is type list and gives the best neighbor
# best_dir is type int and is the steepest direction (key in neighborhood)
# best_init is type dict and contains solved variables for the best point
# improve is type bool and shows if an improvement was made while looking for neighbors
def evaluate_neighbors(ext_vars, init, fmin, tol=0.000001):
    improve = False
    best_var = ext_vars[0]
    here = ext_vars[0]
    best_dir = 0
    best_init = init
    temp = ext_vars
    temp.pop(0, None)
    objectives = {}
    feasibles = {}
    initials = {}
    for i in temp.keys():
        m, status, new_init = build_column(min_trays=8,max_trays=NT,xD=0.95,xB=0.95,x_input=temp[i],provide_init=True,init=init)
                            #fnlp_gdp(NT, temp[i], provide_init=True, init=init)

        if status == pe.SolverStatus.ok:
            objectives[i] = pe.value(m.obj)
            feasibles[i] = temp[i]
            initials[i] = new_init

    key_min = min(objectives.keys(), key=(lambda k: objectives[k]))
    min_obj = objectives[key_min]
    mins = 0
    min_points = {}

    for i in objectives.keys():
        if abs(objectives[i] - min_obj) < tol:
            min_points[i] = feasibles[i]
            mins += 1

    if mins > 1:
        ssums = {}
        for i in min_points.keys():
            ssum = 0
            for j in range(len(best_var)):
                ssum += (min_points[i][j] - here[j])**2
            ssums[i] = ssum
        
        key_max = max(ssums.keys(), key=(lambda k: ssums[k]))

        if objectives[key_max] + tol < fmin:
            fmin = objectives[key_max]
            best_var = ext_vars[key_max]
            best_dir = key_max
            best_init = initials[key_max]
            improve = True
    else:
        if objectives[key_min] + tol < fmin:
            fmin = objectives[key_min]
            best_var = ext_vars[key_min]
            best_dir = key_min
            best_init = initials[key_min]
            improve = True

    return fmin, best_var, best_dir, best_init, improve


# Moves from a certain start in a given direction and evaluates it
# start is type list with the actual point
# init is type dict and contains solved variables for the actual point
# fmin is type int and stands for objective at actual point
# direction is type int and is the moving direction direction (key in neighborhood)
# optimize option will discard out of bounds points given by min and max allowed
# tol is type int and stands for numerical tolereance for equallity
# fmin as return is type int and gives the best point objective (between moved and actual)
# best_var is type list and gives the best point (between moved and actual)
# move is type bool and shows if an improvement was made while looking for neighbors
# best_init is type dict and contains solved variables for the best point
def move_and_evaluate(start, init, fmin, direction, optimize=True, min_allowed={}, max_allowed={}, tol=0.000001):
    best_var = start
    best_init = init
    moved = False

    moved_point = list(map(sum, zip(list(start), list(direction))))

    if optimize:
        checked = 0
        for j in range(len(moved_point)):
            if moved_point[j] >= min_allowed[j+1] and moved_point[j] <= max_allowed[j+1]:
                checked += 1
        if checked == len(moved_point):
            m, status, new_init = build_column(min_trays=8,max_trays=NT,xD=0.95,xB=0.95,x_input=moved_point,provide_init=True,init=init)
            #                   fnlp_gdp(NT, moved_point, provide_init=True, init=init)
            if status == pe.SolverStatus.ok:
                act_obj = pe.value(m.obj)
                if act_obj + tol < fmin:
                    fmin = act_obj
                    best_var = moved_point
                    best_init = new_init
                    moved = True
    else:
        m, status, new_init = build_column(min_trays=8,max_trays=NT,xD=0.95,xB=0.95,x_input=moved_point,provide_init=True,init=init)
        #                      fnlp_gdp(NT, moved_point, provide_init=True, init=init)
        if status == pe.SolverStatus.ok:
            act_obj = pe.value(m.obj)
            if act_obj + tol < fmin:
                fmin = act_obj
                best_var = moved_point
                best_init = new_init
                moved = True

    return fmin, best_var, moved, best_init


def dsda(NT, k='inf', visualize=False):
    print('\n Starting D-SDA with k =',k)
    # Initialize
    t_start = time.process_time()
    route = []
    ext_var = [16,2]
    route.append(ext_var)
    m, _, init = build_column(min_trays=8,max_trays=NT,xD=0.95,xB=0.95,x_input=ext_var,provide_init=False,init={})
    fmin = pe.value(m.obj)
    min_allowed = {i: 2 for i in range(1, len(ext_var)+1)}
    max_allowed = {i: NT-1 for i in range(1, len(ext_var)+1)}

    # Define neighborhood
    if k == '2':
        neighborhood = neighborhood_k_eq_2(len(ext_var))
    elif k == 'inf':
        neighborhood = neighborhood_k_eq_inf(len(ext_var))
    elif k == 'l_flat':
        neighborhood = {1: [1, 1], 2: [-1, -1], 3: [1, 0], 4: [-1, 0], 5: [0, 1], 6: [0, -1]}
    elif k == 'm_flat':
        neighborhood = {1: [1, -1], 2: [1, 0], 3: [-1, 1], 4: [0, 1], 5: [-1, 0], 6: [0, -1]}
    else: 
        return "Enter a valid neighborhood ('inf', '2', 'l_flat' or 'm_flat')"

    looking_in_neighbors = True

    # Look in neighbors (outter cycle)
    while looking_in_neighbors:

        # Find neighbors of the actual point
        neighbors = my_neighbors(ext_var, neighborhood, optimize=True,
                                 min_allowed=min_allowed, max_allowed=max_allowed, cheating=False)

        # Evaluate neighbors of the actual point
        fmin, best_var, best_dir, best_init, improve = evaluate_neighbors(
            neighbors, init, fmin)

        # Stopping condition in case there is no improvement amongst neighbors
        if improve == True:
            line_searching = True
            route.append(best_var)

            # If improvement was made start line search (inner cycle)
            while line_searching:

                # Move in given direction and evaluate
                fmin, best_var, moved, best_init = move_and_evaluate(
                    best_var, best_init, fmin, neighborhood[best_dir], optimize=True, min_allowed=min_allowed, max_allowed=max_allowed)

                # Stopping condition in case no movement was done
                if moved == True:
                    route.append(best_var)
                else:
                    ext_var = best_var
                    line_searching = False

        else:
            looking_in_neighbors = False

    t_end = round(time.process_time() - t_start, 2)

    if visualize:
        visualization(NT, route)

    # Return visited points / final point / objective at that point / execution time
    return route[-1], round(fmin, 5), t_end


if __name__ == "__main__":
    NT = 17
    k = 'inf'  # or k = '2'
    #complete_enumeration(NT)
    print(dsda(NT, k, visualize=False))
    #m, _, _, = build_column(min_trays=8,max_trays=NT,xD=0.95,xB=0.95,x_input=[13,4],provide_init=False,init={})
    #print('Objective',round(pe.value(m.obj),2))


