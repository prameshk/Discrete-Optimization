import numpy as np
from gurobipy import *
import time


###########################################################################################################
def generateFacilityLocationData(C, F):
    # Unbounded ray instance seed 159
    np.random.seed(1234)
    p =  np.random.randint(1000, size=(C, F))
    f = np.random.randint(1000, size=(F))
    for j in range(F):
        for i in range(C):
            f[j] += round(0.05*p[i,j])

    return(C, F, p, f)

def FacilityLocationHeuristic(C, F, p, f):
    # Step 1: Initalize output variables
    x=np.zeros(F); y=np.zeros((C,F));

    # Step 2: Find the best assignment of customer i to a Facility
    #         Implement assignment and open facility
    for i in range(C):        
        idx = np.where(p[i,:] == np.amin(p[i,:]))[0][0]
        y[i,idx]=1;
        x[idx]=1;

    # Step 3: Compute objective value of computed solution
    val=0;
    for j in range(F):
        val=val-f[j]*x[j];
        for i in range(C):
            val=val+p[i,j]*y[i,j];

    # Step 4: Ouput heuristic solution and value
    return(val, x, y);

def ComputeLagrangian(C,F,p,f,lamb):

    # Step 1: Initalize ouput variables
    sublambda=np.zeros(C); x=np.zeros(F);  y=np.zeros((C,F));

    # Step 2: Compute Lagrangian value
    zlambda=0;
    for j in range(F):
        t=0;
        for i in range(C):
            tt=p[i,j]-lamb[i];
            if (tt>=0):
                y[i,j]=1;
                t+=tt;
        ttt=t-f[j];
        if (ttt>=0):
            zlambda=zlambda+ttt;
            x[j]=1;
        else:
            for i in range(C):
                y[i,j]=0;

    for i in range(C):
        zlambda=zlambda+lamb[i];

    # Step 3: Compute subgradient
    for i in range(C):
        sublambda[i]=1;
        for j in range(F):
            sublambda[i]=sublambda[i]-y[i,j];


    # Step 4: Return output values
    return(zlambda,sublambda,x,y);

def LagrangianHeuristic(C,F,p,f,x,y):

# Step 1: Initialize output variables
    newx=np.zeros(F)
    newy=np.zeros((C,F))

# Step 2: If no facility is open, open one randomly
    nopen=0;
    for j in range(F):
        nopen+=x[j];
    
    if (nopen==0):
        np.random.seed(1234)
        idx = np.random.randint(10)
        x[idx] = 1
    

# Step 3: For each customer, keep a single assignment of the current solution
#         (if multiple assignments are in the lagrangian solution) or create
#         one if no assignment exist ... open facilities as required
    for i in range(C):
        minp=float("inf");minpp=float("inf");
        minidx=-1;minidxpp=-1;
        for j in range(F):
            if (y[i,j]>0.5):
                if (p[i,j]<minp):
                    minp=p[i,j];
                    minidx=j;

            if (x[j]>0.5):
                if (p[i,j]<minpp):
                    minpp=p[i,j];
                    minidxpp=j;
                    
        if (minidx>0):
            newy[i,minidx]=1;
            newx[minidx]=1;
        else:
            newy[i,minidxpp]=1;
            newx[minidxpp]=1;


# Step 4: Compute the value of the heuristic solution created
    obj=0;
    for j in range(F):
        obj-=f[j]*newx[j];
        for i in range(C):
            obj+=p[i,j]*newy[i,j];


# Step 5: Return output variables
    return obj,newx,newy;

def SubgradientStep(zlambda, C, lamb, sublambda, rho, K, g, gstar, age):

    # Step 1: initialize working variables
    alpha=0; done='n'; newlambda=np.zeros(C);

    # Step 2: If solution has not improved, increase its age
    if (zlambda>=g):
        age=age+1;
    else:
        age=0;
        g=zlambda;

    # Step 3: If best solution is old, reduce rho
    if (age>K):
        age=0;
        rho=rho/2;
        g=zlambda;
    
    print("        ... rho= {:5.2f},".format(rho))

    # Step 4: Compute square norm of subgradient
    s=0;
    for i in range(C):
        s+= sublambda[i]*sublambda[i];

    print(" s= {:5.2f},".format(s));

    # Step 5: If subgradient is zero, solution is optimal
    if (s<=1e-6):
        done='y';
        s=1e+20;

    # Step 6: Compute step size
    alpha=(zlambda-gstar)*rho/s;
    print(" alpha = {:5.2f} ...\n".format(alpha))

    # Step 7: Perform solution update
    for i in range(C):
        newlambda[i]=lamb[i]-alpha*sublambda[i];

    # Step 8: Return new solution and updated algo parameters
    return newlambda, rho, g, age, done;


###########################################################################################################
# Step 1: Initalize algorithms parameters
rho = 1; K = 5; age = 10; g = float("inf"); maxt = 10000; C = 10; F = 10; lamb = np.zeros(C);

# Step 2: Generate problem instance
C, F, p, f = generateFacilityLocationData(C,F)

# Step 3: Create a greedy heuristic solution
gstar, xx, yy = FacilityLocationHeuristic(C,F,p,f)

# Step 4: Perform maxt iterations of the subgradient algorithm
for t in range(maxt):
    # Step 4.2: Compute lagragian relaxation based on lambda
    zlambda,sublambda,x,y = ComputeLagrangian(C,F,p,f,lamb)
    print("# {:5d}: ".format(t))
    print("UB = {:8.3f} ".format(zlambda))

    # Step 4.3: Modify the lagrangian solution into a feasible solution
    obj,xx,yy = LagrangianHeuristic(C,F,p,f,x,y)
    if (obj>gstar):
        gstar=obj;

    print("/ LB = {:8.3f}".format(gstar))
    if ((zlambda-gstar)<=1e-4):
        print("\nCurrent solution is optimal \n");
        show(xx);print("\n");
        break;



    # Step 4.4: Perform a subgradient step
    lamb, rho, g, age, done = SubgradientStep(zlambda, C, lamb, sublambda, rho, K, g, gstar, age)
    if (done=='y'):
        print("\nCurrent solution is optimal \n");
        print(xx);print(yy);print(II),print("\n");
        break;
