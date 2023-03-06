'''
This script is written by Pramesh Kumar for outer approximation algorithm

Problem: This example is taken from https://optimization.cbe.cornell.edu/index.php?title=Outer-approximation_(OA)

    
MINIMIZE      y_1 + y_2 + x_1^2 + x_2^2
s.t.          (x_1-2)^2 - x_2 \le 0
              x_1 - 2y_1 \ge 0
              x_1 - x_2 - 3(1-y_1) \ge 0
              x_1 + y_1 - 1 \ge 0
              x_2 - y_2 \ge 0
              x_1 + x_2 \ge 3y_1
              y_1 + y_2 \ge 1
              0 \le x_1 \le 4
              0 \le x_2 \le 4
              y_1 binary, y_2 binary
'''

import sys
import numpy as np
from gurobipy import *
import time

############################################################################################################################
##############################################################################################################################

    


def solveModelGurobi():
    m2 = Model()
    y_1 = m2.addVar(lb=0, vtype=GRB.BINARY)
    y_2 = m2.addVar(lb=0, vtype=GRB.BINARY)
    x_1 = m2.addVar(lb=0, vtype=GRB.CONTINUOUS)
    x_2 = m2.addVar(lb=0, vtype=GRB.CONTINUOUS)
    
    m2.addConstr((x_1 - 2)*(x_1 - 2) - x_2 <= 0)
    m2.addConstr((x_1 - 2*y_1) >= 0)
    m2.addConstr((x_1 - x_2) -3*(1-y_1) >= 0)
    m2.addConstr(x_1 + y_1 >= 1)
    m2.addConstr(x_1 - y_2 >= 0)
    m2.addConstr((x_1 + x_2) >= 3*y_1)
    m2.addConstr(y_1 + y_2 >= 1)
    m2.addConstr(x_1 <= 4)
    m2.addConstr(x_2 <= 4)
    obj = y_1 + y_2 + x_1*x_1 + x_2 * x_2
    m2.setObjective(obj, sense=GRB.MINIMIZE)
    m2.update()
    m2.Params.OutputFlag = 0
    m2.setParam('OutputFlag', False)
    m2.optimize()
   
    if m2.status == 2:
        return (m2.objVal, y_1.x, y_2.x, x_1.x, x_2.x)

def prepareMasterModel():
    m2 = Model()
    m2._y_1 = m2.addVar(lb=0, vtype=GRB.BINARY)
    m2._y_2 = m2.addVar(lb=0, vtype=GRB.BINARY)
    m2._x_1 = m2.addVar(lb=0, vtype=GRB.CONTINUOUS)
    m2._x_2 = m2.addVar(lb=0, vtype=GRB.CONTINUOUS)
    m2._eta = m2.addVar(lb=0, vtype=GRB.CONTINUOUS)
    
    m2.addConstr((m2._x_1 - 2*m2._y_1) >= 0)
    m2.addConstr((m2._x_1 - m2._x_2) -3*(1-m2._y_1) >= 0)
    m2.addConstr(m2._x_1 + m2._y_1 >= 1)
    m2.addConstr(m2._x_1 - m2._y_2 >= 0)
    m2.addConstr((m2._x_1 + m2._x_2) >= 3*m2._y_1)
    m2.addConstr(m2._y_1 + m2._y_2 >= 1)
    m2.addConstr(m2._x_1 <= 4)
    m2.addConstr(m2._x_2 <= 4)
    obj = m2._y_1 + m2._y_2 + m2._eta
    m2.setObjective(obj, sense=GRB.MINIMIZE)
    m2.update()
    m2._cuts = 0
    m2.Params.LazyConstraints = 1
    return m2
    
def subproblem(y1, y2):
    ms = Model()
    x_1 = ms.addVar(lb=0, vtype=GRB.CONTINUOUS)
    x_2 = ms.addVar(lb=0, vtype=GRB.CONTINUOUS)
    
    ms.addConstr((x_1 - 2)*(x_1 - 2) - x_2 <= 0)
    ms.addConstr((x_1 - 2*y1) >= 0)
    ms.addConstr((x_1 - x_2) -3*(1-y1) >= 0)
    ms.addConstr(x_1 + y1 >= 1)
    ms.addConstr(x_1 - y2 >= 0)
    ms.addConstr((x_1 + x_2) >= 3*y1)
    ms.addConstr(x_1 <= 4)
    ms.addConstr(x_2 <= 4)
    obj = y1 + y2 + x_1*x_1 + x_2 * x_2
    ms.setObjective(obj, sense=GRB.MINIMIZE)
    ms.update()
    ms.Params.OutputFlag = 0
    ms.setParam('OutputFlag', False)
    ms.optimize()
    
    if ms.status == 2:
        return (round(ms.objVal, 4), round(x_1.x, 2), round(x_2.x, 2))

def Outer_Approximation(intial_guess = [1, 1], eps = 0.1, max_iter = 100):
    m = prepareMasterModel(); m.setParam('OutputFlag', False)
    UB = float("inf"); LB = - float("inf"); counter = 0
    while counter <= max_iter and UB - LB >= eps:           
        if counter == 0:
            y1, y2 = intial_guess[0], intial_guess[1]
        else:
            y1, y2 = round(m._y_1.x), round(m._y_2.x)
        obj, x1, x2 = subproblem(y1, y2)
        if obj < UB:
            UB = obj
        # Adding cut for the non-linear objective value
        m.addConstr(m._eta >= x1*x1 + x2*x2 + 2* x1 * (m._x_1 - x1) + 2 * x2 * (m._x_2 - x2))
        
        # Adding cut for the non-linear constraint
        m.addConstr((x1 -2)*(x1 - 2)  + 2*(x1 - 2) * (m._x_1 - x1) - m._x_2 <= 0)
        m.update()  
        m.optimize()
        if m.status == 2:
            if m.objVal > LB:
                LB = m.objVal
            counter += 1
        else:
            print("model cannot be solved")
            print("model status is %d"%m.status)
            sys.exit(0)      
        
        print("Iteration = %d, LB = %f, UB = %f, x1 = %f, x2 = %f, y1 = %d, y2 = %d"%(counter, LB, UB, x1, x2, y1, y2))


def OuterApproximationCallback(model, where):            
    if where == GRB.Callback.MIPSOL:
        model._cur_obj = model.cbGet(GRB.Callback.MIPSOL_OBJBST)
        model._cuts += 1
        y1 = model.cbGetSolution(model._y_1); y2 = model.cbGetSolution(model._y_2)
        obj, x1, x2 = subproblem(y1, y2)
        print("Cut = %d, LB = %f, UB = %f, x1 = %f, x2 = %f, y1 = %d, y2 = %d"%(model._cuts, model._cur_obj, obj, x1, x2, y1, y2))
        model.cbLazy(model._eta >= x1*x1 + x2*x2 + 2* x1 * (model._x_1 - x1) + 2 * x2 * (model._x_2 - x2))
        model.cbLazy((x1 -2)*(x1 - 2)  + 2*(x1 - 2) * (model._x_1 - x1) - model._x_2 <= 0)
        


############################################################################################################################
##############################################################################################################################

obj, y1, y2, x1, x2 = solveModelGurobi()
print("Gurobi solution:")
print("Objective = %f, x1 = %f, x2 = %f, y1 = %d, y2 = %d"%(obj, x1, x2, y1, y2))

print("Outer approximation:")
Outer_Approximation(intial_guess = [1, 1], eps = 0.1, max_iter = 100)


# Outer approximation callback implementation
print("Outer approximation using callback:")
m = prepareMasterModel(); m.setParam('OutputFlag', False)
m._y_1.start = 1; m._y_2.start = 1
m.optimize(OuterApproximationCallback)

