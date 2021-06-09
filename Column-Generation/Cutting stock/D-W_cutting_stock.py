# -*- coding: utf-8 -*-
"""
Created on Mon May 24 06:47:31 2021

@author: Pramesh Kumar
"""

import numpy as np
from gurobipy import *
import time, math


def GenerateCuttingStockData():
    n=12;
    d=[20,13,17,18,20,45,13,50,20,35,20,30];
    w=[29,75,58,32,18,47,24,62,70,33,21,41];
    W=100;

    u=np.zeros(n) # no. of rolls of type i from one big roll
    for i in range(n):
        u[i]=math.floor(W/w[i]);
    

    nraw=0; # total big rolls needed
    for i in range(n):
        nraw+=math.ceil(d[i]/math.floor(W/w[i]))
    

    return n,d,w,u,W,nraw;


def SolveCuttingStockModel(n,d,w,u,W,nraw,types):

    val=-float("inf");
    # Step 0: start clock
    ts = time.time();

    # Step 1: Create the model
    m = Model()
    if (types=="IP"):
        x = {(i,j): m.addVar(lb = 0.0, vtype = GRB.INTEGER, name = str(i) + '_' + str(j)) for i in range(n) for j in range(nraw)}
        y = {j: m.addVar(lb = 0.0, vtype = GRB.BINARY, name = str(j)) for j in range(nraw)}
    elif (types=="LP"):
        x = {(i,j): m.addVar(lb = 0.0, vtype = GRB.CONTINUOUS, name = str(i) + '_' + str(j)) for i in range(n) for j in range(nraw)}
        y = {j: m.addVar(lb = 0.0, ub = 1.0, vtype = GRB.CONTINUOUS, name = str(j)) for j in range(nraw)}
    else:
        print("\nModel type ({:s}) is not recognized \n".format(types));
        val=-float("inf");
        return val
    
    obj = sum([y[k] for k in y])

    
    for i in range(n):
        m.addConstr(sum([x[i, j] for j in range(nraw)]) == d[i])
        
    for j in range(nraw):
        m.addConstr(sum([w[i] * x[i, j] for i in range(n)]) <= W*y[j])
    
    m.setObjective(obj, sense=GRB.MINIMIZE); m.update(); m.Params.OutputFlag = 1; m.Params.InfUnbdInfo = 1; m.Params.DualReductions = 0


    # Step 2: Solve the model
    m.optimize()
    

    # Step 3: Display solution if required
    if m.status == 2:        
        print("---------------------------------------\n")
        print("Natural Formulation:",types)
        print(" {:6.2f} \n".format(m.objVal))
        print("---------------------------------------\n")
    
        if (type=="IP"):
            for j in range(nraw):
                valy=y[j].x
                if (valy>1e-4):
                    print("* Pattern : [")
                    for i in range(n):
                        print("{:3d}, ".format(round(x[i,j].x)));
                    print("   .......   used   {:3d} times\n".format(round(y[j].x)))

        else:
            for j in range(nraw):
                valy=y[j].x
                if (valy>1e-4):
                    print("* Pattern : [")
                    for i in range(n):
                        print("{:3.2f}, ".format(x[i,j].x));                    
                    print("   .......   used   {:3f} times\n".format(round(y[j].x)))

    # Step 4: Calculate solution time
    te = time.time();
    Deltat = (te-ts);

    print("-----> Run time:      {:6.4f}\n\n".format(Deltat))
        
    # Step 5: output optimal value;
    return m.objVal


def SolveCuttingStockModelCG(n,d,w,u,W,nraw,type, verbsol):

    # Step 0: start clock
    ts = time.time();

    # Step 1: Initialize variables
    t=1;

    # Step 2: Display solution
    print("---------------------------------------\n")
    print("Column Generation Algorithm: \n")
    print("---------------------------------------\n")

    # Step 3: Generate a first collection of simple patterns
    A=np.zeros((n,n)) # pttern of type i rom pattern j
    for i in range(n):
        A[i,i]=u[i]
        
    

    # Step 4: Create and solve initial master problem
    m = Model()
    x = {i: m.addVar(lb = 0.0, vtype = GRB.CONTINUOUS, name = str(i)) for i in range(n)}
    m.update()
    demConstrs = {i: m.addConstr(sum([A[i,j]*x[j] for j in range(n)]) >= d[i]) for i in range(n)}
    obj = sum([x[j] for j in range(n)])
    m.setObjective(obj, sense=GRB.MINIMIZE); m.update(); m.Params.OutputFlag = verbsol; m.Params.InfUnbdInfo = 1; m.Params.DualReductions = 0
    m.optimize()
    if m.status == 2:
        val = m.objVal
        print("* Iteration {:3.1f}: {:7.3f}       ".format(t,val))

        
    # Step 5: Create and solve initial subproblem
    subm = Model()
    y = {i: subm.addVar(lb = 0.0, vtype = GRB.INTEGER, name = str(i)) for i in range(n)}
    for j in range(n):
        subm.addConstr(sum([w[i]*y[i] for i in range(n)]) <= W)
        
        
    subObj = sum([demConstrs[i].pi * y[i] for i in range(n)])
    subm.setObjective(subObj, sense=GRB.MAXIMIZE); m.update(); m.Params.OutputFlag = verbsol; m.Params.InfUnbdInfo = 1; m.Params.DualReductions = 0

    subm.optimize()
    if subm.status == 2:
        rc=subm.objVal;
        print("// |rc| = {:7.3f} //       Pattern = [".format(rc-1))
        for i in range(n):
            print("{:3.0f},".format(round(y[i].x)))
    
    

    


    # Step 6: Create names and matrices for new columns generated
    newcols = [];
    newpatterns =np.zeros((0,n), int)

    # Step 7: Perform column generation
    while (rc>1+1e-5):
        # Step 7.1: Increase number of patterns added
        t+=1;
        

        # Step 7.2: Add variable to master and record generated pattern
        # Column is used for the added varaibles that are used to modify the constraints
        c = Column()
        for i in demConstrs:
            c.addTerms(y[i].x, demConstrs[i])
        
        z = m.addVar(lb = 0.0, vtype = GRB.CONTINUOUS, name = str(len(x)), column = c)
        m.update()      
        
        newcols.append(z)
        newpatterns = np.append(newpatterns, np.array([[y[k].x for k in y]]), axis=0)
        
        
        # Step 7.3: Solve master model
        m.optimize()
        if m.status == 2:
            val=m.objVal
        print("* Iteration {:3.1f}: {:7.3f}       ".format(t,val))

        # Step 7.4: Set the objective and solve new subproblem
        subObj = sum([demConstrs[i].pi * y[i] for i in range(n)])
        subm.setObjective(subObj, sense=GRB.MAXIMIZE); m.update()
    
        subm.optimize()
        if subm.status == 2:
            rc=subm.objVal;
            print("// |rc| = {:7.3f} //       Pattern = [".format(rc-1))
            for i in range(n):
                print("{:3.0f},".format(round(y[i].x)))
            


    val=m.objVal;

    # Step 8: Write lp model to a file
    m.write("CGmodel.lp")

    # Step 9: Display solution obtained at the end of the CG algorithm
    print("---------------------------------------\n")
    print("Solution (Relaxation): {:7.3f} \n".format(val))
    print("---------------------------------------\n")

    for i in range(n):
        newpat = x[i].x;
        if (newpat>1e-4):
            print("* Pattern {:3.0f}: [".format(i));
            for j in range(n-1):
                print("{:3.0f}, ".format(A[i,j]))
            print("{:3.0f}]   .......   used {:7.3f} times \n".format(A[i,n-1],newpat))


    for i in range(len(newcols)):
        newpat= newcols[i].x;
        if (newpat>1e-4):
            print("* Pattern {:3.0f}: [".format(i+n));
            for j in range(n-1):
                print("{:3.0f}, ".format(round(newpatterns[i,j])))
            print("{:3.0f}]   .......   used {:7.3f} times \n".format(round(newpatterns[i,n-1]),newpat))

    print("---------------------------------------\n")


    # Step 11: Calculate solution time
    te = time.time();
    Deltat = (te-ts);
    if (verbsol>=1):
        print("-----> Run time:      {:6.4f}\n\n".format(Deltat))


    # Step 12: return values
    return(val)






############################################################################################################

# Step 1: Define dimensions and parameters

n,d,w,u,W,nraw = GenerateCuttingStockData();



# Step 2: Solve column generation
vallp3 = SolveCuttingStockModelCG(n,d,w,u,W,nraw,"IP", 0)

# Step 3: 
#vallp1 = SolveCuttingStockModel(n,d,w,u,W,nraw,"LP");
#valip1 = SolveCuttingStockModel(n,d,w,u,W,nraw,"IP");
