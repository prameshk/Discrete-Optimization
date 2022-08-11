# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 17:47:25 2022

@author: Pramesh Kumar
"""
import random, math
from gurobipy import *
import networkx as nx


def subtour_constr(model, where):
    if where == GRB.Callback.MIPSOL:
        x0 = model.cbGetSolution(model._x)
        selected_edges = []
        for k in x0:
            if x0[k] > 0.4:
                if k not in selected_edges and (k[1], k[0]) not in selected_edges:
                    selected_edges.append(k)
        G = nx.Graph()
        G.add_edges_from(selected_edges)
        try:
            cycle = nx.find_cycle(G)
            if len(cycle) < len(nodes):
                cycleNodes = list(set([k[0] for k in cycle] + [k[1] for k in cycle]))
                tmp = 0
                for i in cycleNodes:
                    for j in nodes:
                        if i!= j:
                            if (i,j) not in cycle and (j,i) not in cycle:
                                tmp += model._x[i,j]
                model.cbLazy(tmp >= 2)
                
        except:
            print('no cycle found')
        
    

n = 100 # specify the number of nodes
nodes = range(n) # List of noes

random.seed(1)
points = [(random.randint(0, 100), random.randint(0, 100)) for i in range(n)]

# cost of moving from i to j
cost = {(i, j):
        math.sqrt(sum((points[i][k]-points[j][k])**2 for k in range(2)))
        for i in nodes for j in nodes}
    
m = Model()
m._x = {(i, j): m.addVar(vtype=GRB.BINARY) for i in nodes for j in nodes if i != j}

for (i,j) in m._x:
    m.addConstr(m._x[i, j] == m._x[j, i])
    
for i in nodes:
    tmp = sum([m._x[i, j] for j in nodes if i != j])
    m.addConstr(tmp == 2)
    
obj = sum([cost[i,j]*m._x[i,j] for (i,j) in cost if i!= j])

m.setObjective(obj, sense = GRB.MINIMIZE)
m.Params.lazyConstraints = 1

m.optimize(subtour_constr)


if m.status == 3:
    m.computeIIS()
    m.write("model.ilp")
elif m.status == 2:
    selected_edges = []
    for k in m._x:
        if m._x[k].x > 0.4:
            if k not in selected_edges and (k[1], k[0]) not in selected_edges:
                selected_edges.append(k)
    print(selected_edges)

