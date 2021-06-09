# -*- coding: utf-8 -*-
"""
Created on Wed Dec 30 22:26:17 2020

@author: Pramesh Kumar
"""
import numpy as np
import time

def generateFacilityLocationData(C, F):
    # Unbounded ray instance seed 159
    np.random.seed(15645)
    p =  np.random.randint(1000, size=(C, F))
    f = np.random.randint(1000, size=(F))
    for j in range(F):
        for i in range(C):
            f[j] += round(0.05*p[i,j])

    return C, F, p, f


# Step 1: Initialize variables
C = 1000
F = 10



# Step 2: Start clock
ts = time.time()
# Step 3: Generate instance
C, F, p, f = generateFacilityLocationData(C,F)

############################################################################################################################
##############################################################################################################################

def solveModelCplex():
    from docplex.mp.model import Model
    m2 = Model(name='benders', log_output=True)

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
    return m2.objVal, xVal, yVal    
