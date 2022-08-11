# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 19:10:24 2022

@author: Pramesh Kumar
"""
import math, time, heapq
import numpy as np
from gurobipy import *


inputLocation = "Chicago Sketch network/"

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
        self.length = float(_tmpIn[3]) # Length in miles
        self.fft = float(_tmpIn[4]) # Free flow travel time (min)
        self.beta = float(_tmpIn[6])
        self.alpha = float(_tmpIn[5])
        self.speedLimit = float(_tmpIn[7])
        #self.toll = float(_tmpIn[9])
        #self.linkType = float(_tmpIn[10])
        self.flow = 0.0
        self.cost =  float(_tmpIn[4])*VOT #float(_tmpIn[4])*(1 + float(_tmpIn[5])*math.pow((float(_tmpIn[7])/float(_tmpIn[2])), float(_tmpIn[6])))
        self.construction_cost =  float(_tmpIn[3])*constr_cost_per_mi
        self.cost_derivative = 0.0 # For CFW 
        self.target = 0.0  # For CFW

class OD:
    '''
    This class has attributes associated to origin-destination pair
    '''
    def __init__(self, _tmpIn):
        self.origin = _tmpIn[0]
        self.dest = _tmpIn[1]
        self.dem = float(_tmpIn[2])




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
    
def readDemand():
    inFile = open(inputLocation + "demand.dat")
    tmpIn = inFile.readline().strip().split("\t")
    for x in inFile:
        tmpIn = x.strip().split("\t")
        if (tmpIn[0], tmpIn[1]) not in tripSet:
            tripSet[tmpIn[0], tmpIn[1]] = OD(tmpIn)
    

        
###########################################################################################################################

def mcfndp_gurobi():
    print("---------------------------------------\n")
    print("Solving using gurobi: \n")
    print("---------------------------------------\n")
    m = Model()
    x = {(k[0], k[1], o, d): m.addVar(lb = 0.0, vtype = GRB.CONTINUOUS, name = str(k[0]) + "," + str(k[1])) for k in linkSet for (o,d) in tripSet}
    
    m.update()
    for (o,d) in tripSet:
        for i in nodeSet:
            if i == o:
                m.addConstr(sum([x[i, j] for j in nodeSet[i].outLinks]) - sum([x[j, i] for j in nodeSet[i].inLinks])  == tripSet[o,])
            elif i == d:
                m.addConstr(sum([x[i, j] for j in nodeSet[i].outLinks]) - sum([x[j, i] for j in nodeSet[i].inLinks]) == 0)
            else:                        
                m.addConstr(sum([x[i, j] for j in nodeSet[i].outLinks]) - sum([x[j, i] for j in nodeSet[i].inLinks]) == 0)

    obj = sum([x[k]*linkSet[k].fft for k in linkSet])    
    m.setObjective(obj, sense=GRB.MINIMIZE); m.update(); m.Params.OutputFlag = 0; m.Params.InfUnbdInfo = 1; m.Params.DualReductions = 0
    m.optimize()
    print('Final path found')
    print({l for l in linkSet if x[l].x != 0})
    print('ObjVal', m.objVal)

###########################################################################################################################
readStart = time.time()
VOT = 23 # $/hr
constr_cost_per_mi = 5e6 # $/mi


tripSet = {}
linkSet = {}
nodeSet = {}
zoneSet = {}



readNetwork()
readDemand()


print("Reading the network data took", round(time.time() - readStart, 2), "secs")
###########################################################################################################################
