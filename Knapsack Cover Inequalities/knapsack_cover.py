# -*- coding: utf-8 -*-
"""
Solving the Knapsack problem using cover inequalities

Created on Tue Jun  8 07:44:14 2021

@author: Pramesh Kumar
"""

import numpy as np
from gurobipy import *
import time



def generateData():
    np.random.seed(600)
    a = np.random.randint(1000, size=(n))
    w = np.random.randint(1000, size=(n))
    b = np.sum(w)*0.60
    return a, w, b

def solveKnapsackUsingGurobi(n, a, w, b, verbose = 1):
    m = Model()
    x = {i : m.addVar(lb=0, vtype=GRB.BINARY) for i in range(n)}
    m.addConstr(sum([w[j] * x[j] for j in range(n)]) <= b)
    m.setObjective(sum([a[j] * x[j] for j in range(n)]), sense = GRB.MAXIMIZE)
    m.Params.OutputFlag = 0; m.optimize()
    if m.status == 2:  
        if verbose == 1:            
            print('*********************** Gurobi ***************************')
            print('Items selected: ', {k for k in x if x[k].x > 0.4})
            print('Total value: ', m.objVal)
            print('**********************************************************')
        return {k for k in x if x[k].x > 0.4}, m.objVal
    else:
        print('Infeasible model')
        return {}, -float("inf")
    
    
def preparePartialModel(n, a, w, b):
    m = Model()
    x = {i : m.addVar(lb=0.0, ub = 1.0, vtype=GRB.CONTINUOUS, name = str(i)) for i in range(n)}
    m.addConstr(sum([w[j] * x[j] for j in range(n)]) <= b)

    m.addConstr(sum([x[j] for j in range(n)]) <= n-1)
    m.setObjective(sum([a[j] * x[j] for j in range(n)]), sense = GRB.MAXIMIZE)
    m.Params.OutputFlag = 0; m.Params.IntegralityFocus = 1
    return m

def generateCover(n, a, w, b, x):
     m1 = Model()
     z = {i:m1.addVar(lb=0.0, vtype=GRB.BINARY) for i in range(n)}
     m1.addConstr(sum([w[j] * z[j] for j in range(n)]) >= b + 1)
     m1.setObjective(sum([(1-x[j]) * z[j] for j in range(n)]), sense = GRB.MINIMIZE)
     m1.Params.OutputFlag = 0; m1.optimize()
     return [j for j in range(n) if z[j].x > 0.4], m1.objVal
 
def liftCover(c, n, a, w, b):
    alpha = {k:0.0 for k in range(n) if k not in c}
    m2 = Model()
    x_dash = {i:m2.addVar(lb=0.0, vtype=GRB.BINARY) for i in c}
    toBechecked = list(alpha.keys()); checked = []
    if len(toBechecked) != 0:        
        i = toBechecked.pop(0); 
        x_dash[i] = m2.addVar(lb=0.0, vtype=GRB.BINARY)
        m2.update()
        tmpc = sum([x_dash[j] * w[j] for j in c]) + sum([alpha[j]*x_dash[j] for j in checked])
        tmpo = sum([x_dash[j]  for j in c]) + sum([alpha[j]*x_dash[j] for j in checked])
        m2.setObjective(tmpo, sense = GRB.MAXIMIZE)
        m2.addConstr(tmpc <= b - w[i])
        m2.Params.OutputFlag = 0; m2.optimize()       
        if len(c) - 1 - m2.objVal >= 1:
            alpha[i] = len(c) - 1 - m2.objVal
            checked.append(i)
            
    return alpha
    

         
    
def solveKnapsackUsingCoverInequalities(n, a, w, b, verbose = 1):
    m = preparePartialModel(n, a, w, b)
    m.optimize()
    if m.status == 2:
        x = {i:m.getVarByName(str(i)).x for i in range(n)}   
        c, zeta = generateCover(n, a, w, b, x)        
    else:    
        print('Infeasible model')       
    alpha = []
    while zeta < 1:
        m.addConstr(sum([m.getVarByName(str(i)) for i in c]) + sum([m.getVarByName(str(i)) * alpha[i] for i in alpha]) <= len(c) - 1)
        m.optimize()
        x = {i:m.getVarByName(str(i)).x for i in range(n)}   
        c, zeta = generateCover(n, a, w, b, x)  
        alpha = liftCover(c, n, a, w, b)
        zeta = round(zeta, 4)
        #print(c, zeta, alpha)
        
    if verbose == 1:
        print('******************* Cover Inequality *********************')
        print('Items selected: ', {i for i in range(n) if m.getVarByName(str(i)).x > 0.4} )
        print('Total value: ', m.objVal)
        print('**********************************************************')
    return {i for i in range(n) if m.getVarByName(str(i)).x > 0.4}   , m.objVal
        
    
     

################################################################################################
n = 10000 # No. of items
a, w, b = generateData()

start = time.time()
xg, objg = solveKnapsackUsingGurobi(n, a, w, b, verbose = 1)
print('Gurobi model took', time.time()-start); start = time.time()
xc, objc = solveKnapsackUsingCoverInequalities(n, a, w, b, verbose = 1)
print('CG model took', time.time()- start)
