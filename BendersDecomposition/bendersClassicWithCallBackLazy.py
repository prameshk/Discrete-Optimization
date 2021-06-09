'''
This script is written by Pramesh Kumar for testing various acceleration techniques for Benders Decomposition

Problem:

    
MAXIMIZE      \sum_{i \in C} \sum_{j \in F} p_{ij} y_{ij} - \sum_{j \in F} x_j
s.t.          \sum_{j \in F} y_{ij} = 1, \forall i \in C
              0 <= y_{ij} <= x_j, \forall i \in C, \forall j in F
              x binary, y binary
'''


import numpy as np
from gurobipy import *
import time






def generateFacilityLocationData(C, F):
    # Unbounded ray instance seed 159
    np.random.seed(3501)
    p =  np.random.randint(1000, size=(C, F))
    f = np.random.randint(1000, size=(F))
    for j in range(F):
        for i in range(C):
            f[j] += round(0.05*p[i,j])

    return C, F, p, f


# Step 1: Initialize variables
C = 100
F = 10



# Step 2: Start clock
ts = time.time()
# Step 3: Generate instance
C, F, p, f = generateFacilityLocationData(C,F)

############################################################################################################################
##############################################################################################################################

    


def solveModelGurobi():
    m2 = Model()
    x = {j: m2.addVar(lb=0, vtype=GRB.BINARY) for j in range(F)}
    y = {(i, j): m2.addVar(lb=0, vtype=GRB.BINARY) for i in range(C) for j in range(F)}
    for i in range(C):
        m2.addConstr(sum([y[i, j] for j in range(F)]) == 1)
    for j in range(F):
        for i in range(C):
            m2.addConstr(y[i, j] <= bigM*x[j])
    
    obj = sum([p[i, j] * y[i, j]  for j in range(F) for i in range(C)]) -sum([f[j] * x[j] for j in range(F)])

    m2.setObjective(obj, sense=GRB.MAXIMIZE)
    m2.update()
    m2.Params.OutputFlag = 0
    m2.setParam('OutputFlag', False)
    m2.optimize()
    xVal= [x[j].x for j in range(F)]
    yVal =[y[i, j].x for i in range(C) for j in range(F)]
    return m2.objVal, xVal, yVal


def subProblem(x):
    m1 = Model()
    y = {(i, j):  m1.addVar(lb=0, vtype=GRB.CONTINUOUS) for i in range(C) for j in range(F)}
    constrMu = {}
    constrNu = {}
    for i in range(C):
        constrMu[i] = m1.addConstr(sum([y[i, j] for j in range(F)]) == 1)
    for j in range(F):
        for i in range(C):
            constrNu[i, j] = m1.addConstr(y[i, j] <= bigM*x[j])

    obj = sum([p[i, j] * y[i, j]  for j in range(F) for i in range(C)]) -sum([f[j] * x[j] for j in range(F)])

    m1.setObjective(obj, sense=GRB.MAXIMIZE)

    m1.update()
    m1.Params.OutputFlag = 0
    m1.Params.InfUnbdInfo = 1
    m1.Params.DualReductions = 0

    m1.optimize()

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
    
    

            
            


def solveMaster(m,  optCuts_mu, optCuts_nu, fesCuts_mu, fesCuts_nu):

    '''
    Adding optimality cut
    '''
    
    if len(optCuts_nu) != 0:            
        tot = sum(optCuts_mu.values())    
        for j in range(F):
            for i in range(C):
                tot += optCuts_nu[i, j]*m.getVarByName(str(j))*bigM
    
        m.addConstr(m.getVarByName('eta') <= tot)
    '''
    Adding feasibility cut
    '''
    if len(fesCuts_mu) != 0: 
        tot = sum(fesCuts_mu.values())
        for j in range(F):
            for i in range(C):
                tot += fesCuts_nu[i, j] * m.getVarByName(str(j))* bigM
        m.addConstr (tot >= 0)
    

    
    m.optimize()
    if m.status == GRB.OPTIMAL:
        return m.objVal, [m.getVarByName(str(k)).x for k in range(F)], m.getVarByName('eta'), m
    else:
        print("Original model is infeasible", m.status)


def setupMasterProblemModel():
    m = Model()
    eta =  m.addVar(vtype=GRB.CONTINUOUS, ub = 1e5, name ='eta')
    x = {j: m.addVar(lb=0, vtype=GRB.BINARY, name = str(j)) for j in range(F)}
    
    m.setObjective(eta -sum([f[j] * x[j] for j in range(F)]), sense=GRB.MAXIMIZE)
    m.update()
    m.Params.OutputFlag = 0
    m.Params.lazyConstraints = 1
    return m


def callBackFunction(model, where):
    if where == GRB.Callback.MIPSOL: 
        xHat = [model.cbGetSolution(model.getVarByName(str(k))) for k in range(F)]
        print(xHat)
        UB = model.cbGet(GRB.Callback.MIPSOL_OBJ)
        LB, mu, nu, y, status = subProblem(xHat)
       
        if status == 3:
            tot = sum(mu.values())
            for j in range(F):
                for i in range(C):
                    tot += nu[i, j] * model.getVarByName(str(j))* bigM
            model.cbLazy(tot >= 0)
            
            
            
        if status == 2:            
            if round(LB) != round(UB):                    
                tot = sum(mu.values())  
                for j in range(F):
                    for i in range(C):
                        tot += nu[i, j]*model.getVarByName(str(j))*bigM  
                model.cbLazy(model.getVarByName('eta') <= tot)
            else:
                pass
                
                

            
       
        

def solveUFLBenders(eps, x_initial, maxit, verbose=0):
    UB = float("inf")
    LB = -float("inf")
    optCuts_mu = []
    optCuts_nu = []
    fesCuts_mu = []
    fesCuts_nu = []
    tol = float("inf")
    it = 0
    x = x_initial    
    m = setupMasterProblemModel()
    while eps < tol  and it < maxit :
        ob, mu, nu, y, status = subProblem(x)
        LB = max(LB, ob)
        if status == 2:
            obj, x, eta, m = solveMaster(m, mu, nu, [], [])
        else:
            obj, x, eta, m = solveMaster(m,  [], [], mu, nu)        
        #UB = min(UB, obj)
        UB = obj
        tol = UB - LB
        it += 1

        if verbose == 1:            
            print('---------------------- iteration '  + str(it) +' -------------------' )
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
    return x, y, LB

def runCallBackBenders():
    m = setupMasterProblemModel()
    m.optimize(callBackFunction)
    x = [m.getVarByName(str(j)).x for j in range(F)]
    ob, mu, nu, y, status = subProblem(x)
    if round(ob) == round(m.ObjVal):
        return x, y, m.ObjVal
    else:
        print('Callback procedure failed')
        print ('LB = ', ob, ', UB = ', m.ObjVal)
        return x, y, m.ObjVal

def checkGurobiBendersSimilarity(xb, yb, xg, yg):

    ind = 0
    for j in range(F):
        if xb[j] != xg[j]:
            ind += 1
    k = 0
    for j in range(F):
        for i in range(C):
            if yb[k] != yg[k]:

                ind += 1
            k += 1
    if ind == 0:
        print('Solution obtained from both methods are same')
    else:
        print('Solution obtained from both methods are different!!')


bigM = 100000000
x_initial = np.ones(F)
x_initial[1] = 1
x_initial[2] = 0
start = time.time()
xc, yc, obcb = runCallBackBenders()  
print("Classic Benders took with callbacks took...", round(time.time() - start, 2), "seconds")
start = time.time()
xb, yb, obb = solveUFLBenders(0, x_initial, 1000, 1)
print("Classic Benders took...", round(time.time() - start, 2), "seconds")



start = time.time()

obg, xg, yg = solveModelGurobi()
print("Gurobi took...", round(time.time() - start, 2), "seconds")
checkGurobiBendersSimilarity(xb, yb, xg, yg)



