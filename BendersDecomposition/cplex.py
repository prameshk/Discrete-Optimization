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
C = 5
F = 5


# Step 2: Start clock
ts = time.time()
# Step 3: Generate instance
C, F, p, f = generateFacilityLocationData(C,F)

############################################################################################################################
##############################################################################################################################
def solveModelCplex():
    from docplex.mp.model import Model
    m2 = Model(name='cplex', log_output=False)

    x = {j: m2.binary_var() for j in range(F)}
    y = {(i, j): m2.continuous_var(lb=0) for i in range(C) for j in range(F)}
    for i in range(C):
        m2.add_constraint(sum([y[i, j] for j in range(F)]) == 1)
    for j in range(F):
        for i in range(C):
            m2.add_constraint(y[i, j] <= bigM*x[j])
    obj = 0
    for j in range(F):
        obj = obj -f[j] * x[j]
        for i in range(C):
            obj += p[i, j] * y[i, j]
    m2.maximize(obj)
    msol = m2.solve() 
    assert msol is not None, "model can't solve"
    m2.report()
    xVal= [x[j].solution_value for j in range(F)]
    yVal =[y[i, j].solution_value for i in range(C) for j in range(F)]
    return m2.objective_value, xVal, yVal    

def solveModelCplexBenders():
    from docplex.mp.model import Model
    m2 = Model(name='benders_cplex', log_output=False)

    x = {j: m2.binary_var() for j in range(F)}
    y = {(i, j): m2.continuous_var(lb=0) for i in range(C) for j in range(F)}
    for i in range(C):
        m2.add_constraint(sum([y[i, j] for j in range(F)]) == 1)
    for j in range(F):
        for i in range(C):
            m2.add_constraint(y[i, j] <= bigM*x[j])
    obj = 0
    for j in range(F):
        obj = obj -f[j] * x[j]
        for i in range(C):
            obj += p[i, j] * y[i, j]
    m2.maximize(obj)
    m2.parameters.benders.strategy = 3
    msol = m2.solve() 
    assert msol is not None, "model can't solve"
    m2.report()
    xVal= [x[j].solution_value for j in range(F)]
    yVal =[y[i, j].solution_value for i in range(C) for j in range(F)]
    return m2.objective_value, xVal, yVal    


def subProblem(x, m):
    from docplex.mp.model import Model
    m2 = Model(name='subproblem', log_output=False)
    y = {(i, j): m2.continuous_var(lb=0) for i in range(C) for j in range(F)}
    constrMu = {}
    constrNu = {}
    for i in range(C):
        constrMu[i] = m2.add_constraint(sum([y[i, j] for j in range(F)]) == 1)
    for j in range(F):
        for i in range(C):
            constrNu[i, j] = m2.add_constraint(y[i, j] <= bigM*x[j])
            

    obj = sum([p[i, j] * y[i, j]  for j in range(F) for i in range(C)]) -sum([f[j] * x[j] for j in range(F)])

    m2.maximize(obj)
    m2.parameters.lpmethod = 2 #dual simplex
    msol = m2.solve() 
    
    if m2.solve_details.status_code == 1:
        tot = sum(m2.dual_values(constrMu.values()))
        for j in range(F):
            for i in range(C):
                tot += constrNu[i, j].dual_value*m.getVarByName(str(j))*bigM
    
        m.addConstr(m.getVarByName('eta') <= tot)
        for i in range(C):
            tot = sum(optCuts_mu.values())    
            for j in range(F):
                for i in range(C):
                    tot += optCuts_nu[i, j]*m._Model__vars_by_name[str(j)]             
            m.add_constraint(m._Model__vars_by_name('eta') <= tot)

    elif m2.solve_details.status_code == 3:
        if len(fesCuts_mu) != 0:            
            tot = sum(fesCuts_mu.values())
            for j in range(F):
                for i in range(C):
                    tot += fesCuts_nu[i, j] * m.getVarByName(str(j))* bigM
            m.addConstr (tot >= 0)
    
    else:
        print('Model status unknown')
        

    if m1.status == GRB.OPTIMAL:
        mu = {}
        nu = {}
        for i in range(C):
            mu[i] = constrMu[i].pi
        for j in range(F):
            for i in range(C):
                nu[i, j] = constrNu[i, j].pi
        return m1.objVal, mu, nu, [y[i, j].x for i in range(C) for j in range(F)], m1.status
    else:
        mu = {}
        nu = {}
        for i in range(C):
            mu[i] = constrMu[i].FarkasDual
        for j in range(F):
            for i in range(C):
                nu[i, j] = constrNu[i, j].FarkasDual
        return -float("inf"), mu, nu, [], m1.status
    
    

def setupMasterProblemModel():
    from docplex.mp.model import Model
    
    m = Model(name='benders', log_output=False)    
    eta =  m.continuous_var(name ='eta')
    x = {j: m.binary_var(name = str(j)) for j in range(F)}    
    m.maximize(eta -sum([f[j] * x[j] for j in range(F)]))

    return m

def solveMaster(m,  optCuts_mu, optCuts_nu, fesCuts_mu, fesCuts_nu):


    
    m.optimize()
    if m.status == GRB.OPTIMAL:
        return m.objVal, [m.getVarByName(str(k)).x for k in range(F)], m.getVarByName('eta'), m
    else:
        print("Sth went wrong in the master problem and it is ", m.status)





def solveUFLBenders(eps, x_initial, maxit, verbose=0):
    UB = float("inf"); LB = -float("inf"); optCuts_mu = []; optCuts_nu = []; fesCuts_mu = []
    fesCuts_nu = []
    tol = float("inf")
    it = 0
    x = x_initial
    m = setupMasterProblemModel()

    while eps < tol  and it < maxit :
        m, ob = subProblem(x, m)
        LB = max(LB, ob)
        if status == 2:
            #print(status, "optimality")
            obj, x, eta, m = solveMaster(m, mu, nu, [], [])
        else:
            #print(status, "feasibility")
            obj, x, eta, m = solveMaster(m,  [], [], mu, nu)
        
        UB = min(UB, obj)
        tol = UB - LB
        it += 1
        if verbose == 1:            
            print('----------------------iteration '  + str(it) +'-------------------' )
            print ('LB = ', LB, ', UB = ', UB, ', tol = ', tol)
            if len([k for k in range(F) if round(x[k]) != 0]) != 0:
                print('Opened Facilities: \t ', [k for k in range(F) if round(x[k]) != 0])
            else:
                print('No open facilities')
            
            '''
            print('Assignment....')
            for i in range(C):
                print('Customer \t | \t Facility')
                print(str(i)+ '\t | \t ' + str([k for k in for j in range(F) for i in range(C) if y[0] == i and y[k] == 1][0]))
            print(LB, UB, tol, x, it)
            '''
        else:
            continue
    return x, y, UB


##################################################################################################################################
bigM = 1
x_initial = np.zeros(F)
x_initial[1] = 1
x_initial[2] = 1
start = time.time()
'''
xb, yb, obb = solveUFLBenders(100, x_initial, 1000, 1)
print("Benders took...", round(time.time() - start, 2), "seconds")
'''
start = time.time()
obc, xc, yc = solveModelCplex()
print("CPLEX took...", round(time.time() - start, 2), "seconds")
start = time.time()
obcb, xcb, ycb = solveModelCplexBenders()
print("CPLEX Benders took...", round(time.time() - start, 2), "seconds")
#checkGurobiBendersSimilarity(xb, yb, xg, yg)
