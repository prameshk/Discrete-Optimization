# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 22:26:17 2020

@author: Pramesh Kumar
"""
import numpy as np
import time
from gurobipy import *
from math import floor, fabs


import cplex 

from cplex.callbacks import MIPInfoCallback
import cplex.callbacks as CPX_CB


def generateFacilityLocationData(C, F):
    # Unbounded ray instance seed 159
    np.random.seed(15645)
    p =  np.random.randint(1000, size=(C, F))
    f = np.random.randint(1000, size=(F))

    for j in range(F):
        for i in range(C):
            f[j] += float(round(0.05*p[i,j]))

    return C, F, p, f


# Step 1: Initialize variables
C = 5
F = 5
bigM = 1



# Step 2: Start clock
ts = time.time()
# Step 3: Generate instance
C, F, p, f = generateFacilityLocationData(C,F)

class TimeLimitCallback(MIPInfoCallback):

    def __call__(self):
        print(self.get_time())
        if not self.aborted and self.has_incumbent():
            gap = 100.0 * self.get_MIP_relative_gap()
            timeused = self.get_time() - self.starttime
            if timeused > self.timelimit and gap < self.acceptablegap:
                print("Good enough solution at", timeused, "sec., gap =",
                      gap, "%, quitting.")
                self.aborted = True
                self.abort()


class MyBranch(CPX_CB.BranchCallback):

    def __call__(self):
        self.times_called += 1
        br_type = self.get_branch_type()
        if (br_type == self.branch_type.SOS1 or
                br_type == self.branch_type.SOS2):
            return

        x = self.get_values()

        objval = self.get_objective_value()
        obj = self.get_objective_coefficients()
        feas = self.get_feasibilities()

        maxobj = -CPX.infinity
        maxinf = -CPX.infinity
        bestj = -1

        for j in range(len(x)):
            if feas[j] == self.feasibility_status.infeasible:
                xj_inf = x[j] - floor(x[j])
                if xj_inf > 0.5:
                    xj_inf = 1.0 - xj_inf

                if (xj_inf >= maxinf and
                        (xj_inf > maxinf or fabs(obj[j]) >= maxobj)):
                    bestj = j
                    maxinf = xj_inf
                    maxobj = fabs(obj[j])

        if bestj < 0:
            return

        xj_lo = floor(x[bestj])
        # the (bestj, xj_lo, direction) triple can be any python object to
        # associate with a node
        self.make_branch(objval, variables=[(bestj, "L", xj_lo + 1)],
                         node_data=(bestj, xj_lo, "UP"))
        self.make_branch(objval, variables=[(bestj, "U", xj_lo)],
                         node_data=(bestj, xj_lo, "DOWN"))
        # equivalent to
        # self.make_branch(
        #     objval,
        #     constraints=[([[bestj], [1.0]], "G", float(xj_lo + 1))],
        #     node_data=(bestj, xj_lo, "UP"))
        # self.make_branch(
        #     objval,
        #     constraints=[([[bestj], [1.0]], "L", float(xj_lo))],
        #     node_data=(bestj, xj_lo, "DOWN"))



class UserCallback(MIPInfoCallback):
    
    def __call__(self):
        print("I'm here")
        print(self.has_incumbent())


def solveModelCPLEX():
    m = cplex.Cplex()
    facilities = [str(i) for i in range(F)]




    y = list(m.variables.add(obj=(f*-1).tolist(), types=["B"]*len(facilities), names = facilities))

    x = {}
    for i in range(C):
        assignments = [str(i)+","+str(j) for j in range(F)]
        x[i] = m.variables.add(obj = p[i].tolist(), names = assignments)

    assignmentConstrs = {}
    for i in range(C):
        assignmentConstrs[i] = m.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=x[i], val = [1.0]*F)], senses = ["E"], rhs=[1.0])


    indicatorConstrs = {}
    for i in range(C):
        for j in range(F):
            indicatorConstrs[i, j] = m.linear_constraints.add(lin_expr=[cplex.SparsePair(ind = [x[i][0] + j, y[j]], val=[1.0, -bigM])], senses= ["L"], rhs=[0.0])

    m.objective.set_sense(m.objective.sense.maximize)
    #m.set_error_stream(None)
    #m.set_warning_stream(None)
    #m.set_results_stream(None)
    #m.parameters.benders.strategy.set(-1) # For no benders

    #m.parameters.benders.strategy.set(3) Full Benders
    
    """  mastervalue = m.long_annotations.benders_mastervalue
    idx = m.long_annotations.add(
        name=m.long_annotations.benders_annotation,
        defval=mastervalue + 1)
    objtype = m.long_annotations.object_type.variable




    m.long_annotations.set_values(idx, objtype,
                                    [(y[i], mastervalue)
                                        for i in range(F)])
     """

    #m.set_callback(UserCallback)
    branch_instance = m.register_callback(MyBranch)
    branch_instance.times_called = 0

    m.solve()
    print("Solution status =", m.solution.get_status_string())
    print("Optimal value using CPLEX:", m.solution.get_objective_value())
    tol = m.parameters.mip.tolerances.integrality.get()
    values = m.solution.get_values()
    #print(values[:len(y)])


def solveModelGurobi():
    m2 = Model()
    x = {j: m2.addVar(lb=0, vtype=GRB.BINARY) for j in range(F)}
    y = {(i, j): m2.addVar(lb=0, vtype=GRB.BINARY) for i in range(C) for j in range(F)}
    for i in range(C):
        m2.addConstr(sum([y[i, j] for j in range(F)]) == 1)
    for j in range(F):
        for i in range(C):
            m2.addConstr(y[i, j] <= bigM*x[j])
    obj = 0
    for j in range(F):
        obj = obj -f[j] * x[j]
        for i in range(C):
            obj += p[i, j] * y[i, j]
    m2.setObjective(obj, sense=GRB.MAXIMIZE)
    m2.update()
    m2.Params.OutputFlag = 0
    m2.setParam('OutputFlag', False)
    m2.optimize()
    xVal= [x[j].x for j in range(F)]
    yVal =[y[i, j].x for i in range(C) for j in range(F)]
    print("Objective using Gurobipy", m2.objVal)
    

    return m2.objVal, xVal, yVal

start = time.time()
#obg, xg, yg = solveModelGurobi()
print("Gurobi took...", round(time.time() - start, 2), "seconds")

start = time.time()
solveModelCPLEX()
print("CPLEX took...", round(time.time() - start, 2), "seconds")


##################################################################################################################################


'''
x_initial = np.zeros(F)
x_initial[1] = 1
x_initial[2] = 1
start = time.time()

#xb, yb, obb = solveUFLBenders(100, x_initial, 1000, 1)
print("Benders took...", round(time.time() - start, 2), "seconds")

start = time.time()
obc, xc, yc = solveModelCplex()
print("CPLEX took...", round(time.time() - start, 2), "seconds")
start = time.time()
obcb, xcb, ycb = solveModelCplexBenders()
print("CPLEX Benders took...", round(time.time() - start, 2), "seconds")
#checkGurobiBendersSimilarity(xb, yb, xg, yg)

'''
