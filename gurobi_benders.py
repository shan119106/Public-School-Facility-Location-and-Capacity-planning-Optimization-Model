# -*- coding: utf-8 -*-
"""
Created on Mon May 16 12:55:38 2022

@author: Admin
"""
'''sol = master.getObjective().getValue() + obj2.getValue()'''
import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import math
#importing excel data
schl = pd.read_excel("data\school.xlsx")
vil = pd.read_excel("data\Village.xlsx")
Current_stat = pd.read_excel("data\datapoints.xlsx",sheet_name = "CurrentStat")
schl_size = len(schl)
vil_size = len(vil)
SCHOOL=[]
VILLAGE =[]
#school and village index
for i in range(schl_size):
    SCHOOL.append(schl["SCHOOL"][i])
for i in range(vil_size):
    VILLAGE.append(vil["VILLAGE"][i])
temp =[]
distance = 0
#updating the euclidien distance
for i in range(schl_size):
    for j in range(vil_size):
        distance = math.sqrt((vil["VILLAGE.X"][j]-schl["SCHOOL.X"][i])**2 +(vil["VILLAGE.Y"][j]-schl["SCHOOL.Y"][i])**2 )
        temp.append(math.ceil(distance))
    column_name ="DIST" +str(i+1)
    Current_stat[column_name] = temp
    temp =[]
    Current_stat.to_excel("data\datapoints.xlsx",sheet_name = "CurrentStat",index=False)
#variables required
length = 3
Max_dist = 80
alpha = 0.2
beta = 0.3
pi =100
fixed ={}
a ={}
UB = 1000000000
LB = 1
f ={}
M =1000000
Ie = [1,3,4] 
A ={}
B ={}
C={}
D ={}
E ={}
F ={}
G ={}
#fixedCost Calculation
for i in SCHOOL:
    fixed[i] = {}
    f[i] = schl["Fixed"][i-1]
    for j in range(1,length+1):
        fixed[i][j] = schl["Fixed" +str(j)][i-1]
#costtransport calculation
Costtransport ={}
for i in SCHOOL:
    Costtransport[i] ={}
    a[i] ={}
    for j in VILLAGE:
        Costtransport[i][j] = 0.1*Current_stat["DIST"+str(i)][j-1]
        if(Current_stat["DIST"+str(i)][j-1] < 8):
            a[i][j] = 1
        else:
            a[i][j]=0   

iteration =0
#LPMasterProblem
master = gp.Model("Masterproblem")
#LPMasterproblem Variable
Y ={}
for i in SCHOOL:
    Y[i] ={}
    for j in range(1,length+1):
        Y[i][j] = master.addVar(vtype = GRB.BINARY,name ="Y"+str(i)+str(j))
Yc ={}
for i in SCHOOL:
    Yc[i] =master.addVar(vtype = GRB.BINARY,name ="Yc"+str(i))
Ync ={}
for i in SCHOOL:
    Ync[i] = master.addVar(vtype = GRB.BINARY,name = "Ync"+str(i))
Phi = master.addVar(name = "AuxVar",vtype =GRB.CONTINUOUS,lb=0)   
lambda0=0.25
lambda1=0.25
lambda2=0.25
lambda3 =0.25     
master.setObjective(gp.quicksum(lambda0*fixed[i][j]*Y[i][j] for j in range(1,4) for i in SCHOOL)+ gp.quicksum(lambda0*f[i]*Yc[i] for i in SCHOOL)+Phi,GRB.MINIMIZE)
#Constraint1
for i in SCHOOL:
    master.addConstr(gp.quicksum(Y[i][j] for j in range(1,4)) + Yc[i] + Ync[i] == 1)
#constraint5
master.addConstrs(Yc[i] == 0 for i in Ie)
while((UB-LB)/UB >= 0.2):
    iteration  = iteration +1
    #cuts
    for itr in range(1,iteration):
         master.addConstr(Phi >= gp.quicksum(vil["POP"][j-1]*A[itr][j] for j in VILLAGE) + gp.quicksum(B[itr][i]*schl["Qi"][i-1]*Ync[i]+ gp.quicksum(B[itr][i]*schl["Qi"+ str(l)][i-1]*Y[i][l] for l in range(1,4)) for i in SCHOOL) + gp.quicksum(C[itr][i][j]*a[i][j]*vil["POP"][j-1] +D[itr][i][j]*Current_stat[i][j-1] for i in SCHOOL for j in VILLAGE) + gp.quicksum(E[itr][i]*alpha*schl["Qi"][i-1]*Ync[i] + gp.quicksum(E[itr][i]*schl["Qi"+ str(l)][i-1]*Y[i][l] for l in range(1,4)) for i in SCHOOL) + gp.quicksum(F[itr][j]*beta*vil["POP"][j-1] for j in VILLAGE))
    master.update()
    master.optimize()
    #objectValue
    obj = master.getObjective()
    LB = obj.getValue()
    #pspsubproblem
    psp = gp.Model("PrimalSubproblem")
    #pspvariables
    X ={}
    for i in SCHOOL:
        X[i] ={}
        for j in VILLAGE:
            X[i][j] =psp.addVar(lb =0,vtype =GRB.CONTINUOUS,name="X")
    U={}
    for i in SCHOOL:
        U[i] =psp.addVar(lb =0,vtype =GRB.CONTINUOUS,name ="U")
    V={}
    for j in VILLAGE:
        V[j] = psp.addVar(lb =0,vtype =GRB.CONTINUOUS,name ="V")
    rho = psp.addVar(lb =0,vtype =GRB.CONTINUOUS,name ="rho")
    plusS ={}
    for i in SCHOOL:
        plusS[i] ={}
        for j in VILLAGE:
            plusS[i][j] = psp.addVar(lb =0,vtype =GRB.CONTINUOUS,name ="plusS")
    MinusS ={}
    for i in SCHOOL:
        MinusS[i] ={}
        for j in VILLAGE:
            MinusS[i][j] = psp.addVar(lb =0,vtype =GRB.CONTINUOUS,name ="MinusS")
    #probconstraint2
    psp.addConstrs((gp.quicksum(X[i][j] for i in SCHOOL) + V[j] == vil["POP"][j-1]for j in VILLAGE),name = "probconstraint2")
    #probconstraint3
    psp.addConstrs((gp.quicksum(X[i][j] for j in VILLAGE) + U[i] == schl["Qi"][i-1]*Ync[i].x+ gp.quicksum(schl["Qi"+ str(l)][i-1]*Y[i][l].x for l in range(1,4)) for i in SCHOOL),name = "probconstraint3")
    #probconstraint4
    psp.addConstrs((X[i][j] <= a[i][j]*vil["POP"][j-1] for i in SCHOOL for j in VILLAGE),name = "probconstraint4")
    #pspconstraint6
    psp.addConstrs((X[i][j] - plusS[i][j] +MinusS[i][j] == Current_stat[i][j-1] for i in SCHOOL for j in VILLAGE),name = "probconstraint6")
    #probconstraint7
    psp.addConstrs((U[i] <= alpha*(schl["Qi"][i-1]*Ync[i].x + gp.quicksum(schl["Qi"+ str(l)][i-1]*Y[i][l].x for l in range(1,4))) for i in SCHOOL),name = "probconstraint7")
    #probconstraint8
    psp.addConstrs((V[j]<= beta*(vil["POP"][j-1]) for j in VILLAGE),name = "probconstraint8")
    #probconstraint9
    psp.addConstrs((gp.quicksum(Costtransport[i][j]*X[i][j] for i in SCHOOL)-vil["POP"][j-1]*rho <= 0 for j in VILLAGE),name = "probconstraint9")
    #Objective for psp
    psp.setObjective(gp.quicksum(lambda1*Costtransport[i][j]*X[i][j] for j in VILLAGE for i in SCHOOL) + gp.quicksum(lambda2*pi*MinusS[i][j] for j in VILLAGE for i in SCHOOL) + lambda3*rho,GRB.MINIMIZE)
    psp.Params.DualReductions =0
    psp.write("psp.lp")
    psp.optimize()
    print(psp.status)
    obj1 = psp.getObjective()
    UB = obj.getValue() + obj1.getValue()
    A[iteration] ={}
    for j in VILLAGE:
        A[iteration][j] = psp.getConstrByName("probconstraint2"+str([j])).pi
    for j in VILLAGE:
        print(A[iteration][j])
    B[iteration] = {}
    for i in SCHOOL:
        B[iteration][i] = psp.getConstrByName("probconstraint3"+str([i])).pi
    C[iteration]={}
    for i in SCHOOL:
        C[iteration][i] ={}
        for j in VILLAGE:
            C[iteration][i][j] = psp.getConstrByName("probconstraint4" + "["+str(i)+"," +str(j) +"]").pi
    D[iteration] ={}
    for i in SCHOOL:
        D[iteration][i] = {}
        for j in VILLAGE:
            D[iteration][i][j] =psp.getConstrByName("probconstraint6" +"["+str(i)+"," +str(j) +"]").pi
    E[iteration] ={}
    for i in SCHOOL:
        E[iteration][i] = psp.getConstrByName('probconstraint7' +str([i])).pi
    F[iteration] ={}
    for j in VILLAGE:
        F[iteration][j] = psp.getConstrByName("probconstraint8"+str([j])).pi
    G[iteration] ={}
    for j in VILLAGE:
        G[iteration][j] = psp.getConstrByName("probconstraint9"+str([j])).pi
