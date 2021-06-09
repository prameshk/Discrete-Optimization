# -*- coding: utf-8 -*-
"""
Created on Tue May 25 10:45:59 2021

@author: Pramesh Kumar
"""
import math, time, heapq
import numpy as np
from gurobipy import *

inputLocation = "Sioux Falls network/"

class Zone:
    def __init__(self, _tmpIn):
        self.zoneId = _tmpIn[0]
        self.lat = 0
        self.lon = 0
        self.destList = []


class Node:
    '''
    This class has attributes associated with any node
    '''
    def __init__(self, _tmpIn):
        self.Id = _tmpIn[0]
        self.lat = 0
        self.lon = 0
        self.outLinks = []
        self.inLinks = []
        self.label = float("inf")
        self.pred = ""
        self.inDegree = 0
        self.outDegree = 0
        self.order = 0 # Topological order
        self.wi = 0.0 # Weight of the node in Dial's algorithm
        self.xi = 0.0 # Toal flow crossing through this node in Dial's algorithm


class Link:
    '''
    This class has attributes associated with any link
    '''
    def __init__(self, _tmpIn):
        self.tailNode = _tmpIn[0]
        self.headNode = _tmpIn[1]
        self.capacity = float(_tmpIn[2]) # veh per hour
        self.length = float(_tmpIn[3]) # Length
        self.fft = float(_tmpIn[4]) # Free flow travel time (min)
        self.beta = float(_tmpIn[6])
        self.alpha = float(_tmpIn[5])
        self.speedLimit = float(_tmpIn[7])
        #self.toll = float(_tmpIn[9])
        #self.linkType = float(_tmpIn[10])
        self.flow = 0.0
        self.cost =  float(_tmpIn[4]) #float(_tmpIn[4])*(1 + float(_tmpIn[5])*math.pow((float(_tmpIn[7])/float(_tmpIn[2])), float(_tmpIn[6])))
        self.cost_derivative = 0.0 # For CFW 
        self.target = 0.0  # For CFW
        self.logLike = 0.0
        self.reasonable = True # This is for Dial's stochastic loading
        self.wij = 0.0 # Weight in the Dial's algorithm
        self.xij = 0.0 # Total flow on the link for Dial's algorithm
        self.resource = 0.0




def readNetwork():
    inFile = open(inputLocation + "network.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        linkSet[tmpIn[0], tmpIn[1]] = Link(tmpIn)
        if tmpIn[0] not in nodeSet:
            nodeSet[tmpIn[0]] = Node(tmpIn[0])
        if tmpIn[1] not in nodeSet:
            nodeSet[tmpIn[1]] = Node(tmpIn[1])
        if tmpIn[1] not in nodeSet[tmpIn[0]].outLinks:
            nodeSet[tmpIn[0]].outLinks.append(tmpIn[1])
        if tmpIn[0] not in nodeSet[tmpIn[1]].inLinks:
            nodeSet[tmpIn[1]].inLinks.append(tmpIn[0])

    inFile.close()
    print(len(nodeSet), "nodes")
    print(len(linkSet), "links")
    
def assignResources():
    np.random.seed(15645)
    for l in linkSet:
        linkSet[l].resource = np.random.randint(8)
    
###########################################################################################################################
def DijkstraHeap(origin, cost, dual):
    '''
    Calcualtes shortest path from an origin to all other destinations.
    The labels and preds are stored in node instances.
    '''
    for n in nodeSet:
        nodeSet[n].label = float("inf")
        nodeSet[n].pred = ""
    nodeSet[origin].label = 0.0
    nodeSet[origin].pred = "NA"
    SE = [(0, origin)]
    while SE:
        currentNode = heapq.heappop(SE)[1]
        currentLabel = nodeSet[currentNode].label
        for toNode in nodeSet[currentNode].outLinks:
            link = (currentNode, toNode)
            newNode = toNode
            newPred =  currentNode
            existingLabel = nodeSet[newNode].label
            if cost == 'fft':                
                newLabel = currentLabel + linkSet[link].fft
            elif cost == 'resource':
                newLabel = currentLabel + linkSet[link].resource
            elif cost == 'pricing':
                newLabel = currentLabel + (linkSet[link].fft - linkSet[link].resource*dual)
                
            if newLabel < existingLabel:
                heapq.heappush(SE, (newLabel, newNode))
                nodeSet[newNode].label = newLabel
                nodeSet[newNode].pred = newPred
                
def tracePreds(dest):
   '''
   This method traverses predecessor nodes in order to create a shortest path
   '''
   prevNode = nodeSet[dest].pred
   spLinks = []
   while nodeSet[dest].pred != "NA":
       spLinks.append((prevNode, dest))
       dest = prevNode
       prevNode = nodeSet[dest].pred
   return spLinks   

      
def gurobiModel(o, d):
    print("---------------------------------------\n")
    print("Solving using gurobi: \n")
    print("---------------------------------------\n")
    m = Model()
    x = {k: m.addVar(lb = 0.0, vtype = GRB.CONTINUOUS, name = str(k[0]) + "," + str(k[1])) for k in linkSet}
    x['o',o] = m.addVar(lb = 0.0, vtype = GRB.CONTINUOUS, name = "o," + str(o))
    x[d,'d'] = m.addVar(lb = 0.0, vtype = GRB.CONTINUOUS, name = "d," + str(d))
    m.update()
    m.addConstr(x['o',o] == 1); m.addConstr(x[d,'d'] == 1)
    
    for i in nodeSet:
        if i == o:
            m.addConstr(sum([x[i, j] for j in nodeSet[i].outLinks]) - sum([x[j, i] for j in nodeSet[i].inLinks]) - x['o',o] == 0)
        elif i == d:
            m.addConstr(sum([x[i, j] for j in nodeSet[i].outLinks]) + x[d,'d'] - sum([x[j, i] for j in nodeSet[i].inLinks]) == 0)
        else:                        
            m.addConstr(sum([x[i, j] for j in nodeSet[i].outLinks]) - sum([x[j, i] for j in nodeSet[i].inLinks]) == 0)

    m.addConstr(sum([x[k]*linkSet[k].resource for k in linkSet]) <= 4)    
    obj = sum([x[k]*linkSet[k].fft for k in linkSet])    
    m.setObjective(obj, sense=GRB.MINIMIZE); m.update(); m.Params.OutputFlag = 0; m.Params.InfUnbdInfo = 1; m.Params.DualReductions = 0
    m.optimize()
    print('Final path found')
    print({l for l in linkSet if x[l].x != 0})
    print('ObjVal', m.objVal)
def SPRC_CG(o, d):
    val=-float("inf");
    # Step 0: start clock
    ts = time.time();
    
    # Step 1: Initialize variables
    t=1;
    
    # Step 2: Display solution algorithm
    print("---------------------------------------\n")
    print("Column Generation Algorithm: \n")
    print("---------------------------------------\n")

    # Step 3: Generate a first collection of simple path using shortest path alg.
    DijkstraHeap(o, 'resource', 0)
    paths = {0:tracePreds(d)}
    
    # Step 4: 
    m = Model()
    lamb = {i: m.addVar(lb = 0.0, vtype = GRB.CONTINUOUS, name = str(i)) for i in paths}
    m.update()
    tempSum = 0; obj = 0
    for p in paths:
        tempSum += sum([linkSet[k].resource for k in paths[p]])*lamb[p] 
        obj += sum([linkSet[k].fft for k in paths[p]])*lamb[p] 
    resConstr = m.addConstr(tempSum <= 4)
    convConstr = m.addConstr(sum([lamb[p] for p in paths]) == 1)


    m.setObjective(obj, sense=GRB.MINIMIZE); m.update(); m.Params.OutputFlag = 0; m.Params.InfUnbdInfo = 1; m.Params.DualReductions = 0
    m.optimize()
    if m.status == 2:
        val = m.objVal
        print("* Iteration {:3.1f}: {:7.3f}       ".format(t,val))
        
    
    
    # Step 5: Solve initial subproblem using shortest path
    DijkstraHeap(o, 'pricing', resConstr.pi)
    reducedCost = nodeSet[d].label - convConstr.pi
    
    
    while reducedCost < 0-1e-5:
        
        # Step 7.1: Increase number of patterns added
        t+=1;
        
        # Step 7.1: Include new path with the reduced cost
        p =  tracePreds(d)
        paths[len(paths)] = p
        
        
        # Step 7.2: Add lambda variable of new path to master
        # Column is used for the added variables that are used to modify the constraints
        c = Column()
        
        c.addTerms(sum([linkSet[k].resource for k in p]), resConstr)
        c.addTerms(1.0, convConstr)
        
        lamb[len(lamb)] = m.addVar(lb = 0.0, vtype = GRB.CONTINUOUS, name = str(len(lamb)), column = c)
        m.update()
        
        # Step 7.3: Solve master model
        m.optimize()
        if m.status == 2:
            val=m.objVal
        print("* Iteration {:3.1f}: {:7.3f}       ".format(t,val))
    
        # Step 5: Solve initial subproblem using shortest path
        DijkstraHeap(o, 'pricing', resConstr.pi)
        reducedCost = nodeSet[d].label - convConstr.pi
        
    print("Paths added")
    for p in paths:
        print(paths[p],lamb[p].x)
    
    print("Obj val", m.objVal)
        
    
    
      
                
###########################################################################################################################
readStart = time.time()

tripSet = {}
linkSet = {}
nodeSet = {}
zoneSet = {}



readNetwork()
assignResources()

print("Reading the network data took", round(time.time() - readStart, 2), "secs")
###########################################################################################################################
o, d = ('5', '13')
gurobiModel(o, d)
SPRC_CG(o, d)
