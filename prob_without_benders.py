
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 19:24:25 2022

@author: Admin
"""

import pulp as p
import pandas as pd
import math
#importing excel data
schl = pd.read_excel("data\school.xlsx")
vil = pd.read_excel("data\Village.xlsx")
Current_stat = pd.read_excel("data\Current_stat.xlsx")
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
a = 0
#updating the euclidien distance
for i in range(schl_size):
    for j in range(vil_size):
        a = math.sqrt((vil["VILLAGE.X"][j]-schl["SCHOOL.X"][i])**2 +(vil["VILLAGE.Y"][j]-schl["SCHOOL.Y"][i])**2 )
        temp.append(math.ceil(a))
    column_name ="DIST" +str(i+1)
    Current_stat[column_name] = temp
    temp =[]
    Current_stat.to_excel("data\Current_stat.xlsx",index=False)
#variables required
l = 3
Max_dist = 3
alpha = 0.2
beta = 0.3
pi =100
fixed ={}
a ={}
UB = 1000000000
LB = 1
f ={}
#fixedCost Calculation
for i in SCHOOL:
    fixed[i] = {}
    f[i] = schl["Fixed"][i-1]
    for j in range(l):
        fixed[i][j+1] = schl["Fixed" +str(i)][j]
#costtransport calculation
Costtransport ={}
for i in SCHOOL:
    Costtransport[i] ={}
    a[i] ={}
    for j in VILLAGE:
        Costtransport[i][j] = 0.1*Current_stat["DIST"+str(i)][j-1]
        if(Current_stat["DIST"+str(i)][j-1] < 3):
            a[i][j] = 1
        else:
            a[i][j]=0   
iteration = 0
prob = p.LpProblem("Whole_problem",p.LpMinimize)
Y = p.LpVariable.dicts("Y",[(i,j) for i in SCHOOL for j in range(1,4)],cat ="Binary")
Yc = p.LpVariable.dicts("Yc",SCHOOL,cat ="Binary")
Ync = p.LpVariable.dicts("Ync",SCHOOL,cat ="Binary")
X = p.LpVariable.dicts("Villagetoschool",[(i,j) for i in SCHOOL for j in VILLAGE],lowBound=0)
U =p.LpVariable.dicts("Unmet capacity of School",SCHOOL,lowBound =0)
V=p.LpVariable.dicts("Unmet demand of village",VILLAGE,lowBound =0)
rho = p.LpVariable("MxCost",lowBound =0)
plusS = p.LpVariable.dicts("DropinAllo",[(i,j) for i in SCHOOL for j in VILLAGE],lowBound =0)
MinusS = p.LpVariable.dicts("IncrementinAllo",[(i,j) for i in SCHOOL for j in VILLAGE],lowBound=0)
#LpObjectivefunction
lambda0=0.25
lambda1=0.25
lambda2=0.25
lambda3 =0.25
Ie = [1,2]
prob += p.lpSum(lambda0*fixed[i][j]*Y[(i,j)] for j in range(1,4) for i in SCHOOL)+ p.lpSum(lambda0*f[i]*Yc[i] for i in SCHOOL) + p.lpSum(lambda1*Costtransport[i][j]*X[(i,j)] for j in VILLAGE for i in SCHOOL) + p.lpSum(lambda2*pi*MinusS[(i,j)] for j in VILLAGE for i in SCHOOL) + lambda3*rho
#Constraint1
for i in SCHOOL:
    prob += p.lpSum(Y[(i,j)] for j in SCHOOL) + Yc[i] + Ync[i] == 1
#constraint5
for i in Ie:
    prob += Yc[i] == 0
#probconstraint2
for j in VILLAGE:
    prob += p.lpSum(X[(i,j)] for i in SCHOOL) + V[j] == vil["POP"][j-1]
#probconstraint3
for i in SCHOOL:
    prob += p.lpSum(X[(i,j)] for j in VILLAGE) + U[i] == schl["Qi"][i-1]*Ync[i]+ p.lpSum(schl["Qi"+ str(l)][i-1]*Y[(i,l)] for l in range(1,4))
#probconstraint4
for i in SCHOOL:
    for j in VILLAGE:
        prob += X[(i,j)] <= a[i][j]*vil["POP"][j-1]   
#probconstraint6
for i in SCHOOL:
    for j in VILLAGE:
        prob += X[(i,j)] - plusS[(i,j)] +MinusS[(i,j)] == Current_stat[i][j-1] 
#probconstraint7
for i in SCHOOL:
    prob += U[i] <= alpha*(schl["Qi"][i-1]*Ync[i] + p.lpSum(schl["Qi"+ str(l)][i-1]*Y[(i,l)] for l in range(1,4)))
#probconstraint8
for j in VILLAGE:
    prob += V[j]<= beta*vil["POP"][j-1]
#probconstraint9
for j in VILLAGE:
    prob += p.lpSum(Costtransport[i][j]*X[(i,j)] for i in SCHOOL)-vil["POP"][j-1]*rho <= 0
prob.solve()
print(p.value(prob.objective))
for i in SCHOOL:
    print(Yc[i].varValue)
    print(Ync[i].varValue)
    for j in range(1,l+1):
        print(Y[(i,j)].varValue)
    