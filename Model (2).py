#!/usr/bin/env python
# coding: utf-8

# In[3]:


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
#storing Subproblem values
aval ={}
bval ={}
cval ={}
dval ={}
Eval ={}
fval ={}
gval ={}
f ={}
#for storing every cut
cut ={}
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
while ((UB-LB)/UB) >= 0.42: 
    #MasterProblem Initialization
    master = p.LpProblem("Master_masterlem",p.LpMinimize)
    #Parameters
    lambda0=0.25
    lambda1=0.25
    lambda2=0.25
    lambda3 =0.25
    Ie = [1,3]            
    #LPMasterproblem Variable
    Y = p.LpVariable.dicts("Y",[(i,j) for i in SCHOOL for j in range(1,4)],lowBound = 0,upBound =1,cat ="Binary")
    Yc = p.LpVariable.dicts("Yc",SCHOOL,lowBound = 0,upBound =1,cat ="Binary")
    Ync = p.LpVariable.dicts("Ync",SCHOOL,lowBound = 0,upBound =1,cat ="Binary")
    Phi = p.LpVariable("AuxVar",lowBound=0)            
    #LpObjectivefunction
    master += p.lpSum(lambda0*fixed[i][j]*Y[(i,j)] for j in range(1,4) for i in SCHOOL)+ p.lpSum(lambda0*f[i]*Yc[i] for i in SCHOOL) + Phi
    #Constraint1
    for i in SCHOOL:
        master += p.lpSum(Y[(i,j)] for j in SCHOOL) + Yc[i] +Ync[i] == 1
    #constraint5
    for i in Ie:
        master += Yc[i] == 0
    iteration +=1
    #cuts
    for itr in range(1,iteration):
        print(itr)
        print(Eval[itr])
        master += Phi>= p.lpSum(vil["POP"][j-1]*aval[itr][j] for j in VILLAGE) + p.lpSum(bval[itr][i]*(schl["Qi"][i-1]*Ync[i]+ p.lpSum(schl["Qi"+ str(l)][i-1]*Y[(i,l)] for l in range(1,4))) for i in SCHOOL) + p.lpSum(cval[itr][i][j]*a[i][j]*vil["POP"][j-1] +dval[itr][i][j]*Current_stat[i][j-1] for i in SCHOOL for j in VILLAGE) + p.lpSum(Eval[itr][i]*alpha*schl["Qi"][i-1]*Ync[i] + Eval[itr][i]*alpha*schl["Qi"+ str(l)][i-1]*Y[(i,l)] for l in range(1,4) for i in SCHOOL) + p.lpSum(fval[itr][j]*beta*vil["POP"][j-1] for j in VILLAGE)
    #write it in output
    master.writeLP("Master")
    #solving master problem
    master.solve() 
    #setting Lower bound
    LB = max(LB,p.value(master.objective))
    #primal problem
    psp = p.LpProblem("Facility_Location",p.LpMinimize)
    #primal variables
    X = p.LpVariable.dicts("Villagetoschool",[(i,j) for i in SCHOOL for j in VILLAGE],lowBound=0)
    U =p.LpVariable.dicts("Unmet capacity of School",SCHOOL,lowBound =0)
    V=p.LpVariable.dicts("Unmet demand of village",VILLAGE,lowBound =0)
    rho = p.LpVariable("MxCost",lowBound =0)
    plusS = p.LpVariable.dicts("DropinAllo",[(i,j) for i in SCHOOL for j in VILLAGE],lowBound =0)
    MinusS = p.LpVariable.dicts("IncrementinAllo",[(i,j) for i in SCHOOL for j in VILLAGE],lowBound=0)
    #psp objective function
    psp += p.lpSum(lambda1*Costtransport[i][j]*X[(i,j)] for j in VILLAGE for i in SCHOOL) + p.lpSum(lambda2*pi*MinusS[(i,j)] for j in VILLAGE for i in SCHOOL) + lambda3*rho
    
    #pspconstraint2
    for j in VILLAGE:
        psp += p.lpSum(X[(i,j)] for i in SCHOOL) + V[j] == vil["POP"][j-1]
    #pspconstraint3
    for i in SCHOOL:
        psp += p.lpSum(X[(i,j)] for j in VILLAGE) + U[i] == schl["Qi"][i-1]*Ync[i].varValue+ p.lpSum(schl["Qi"+ str(i)][l-1]*Y[(i,l)].varValue for l in range(1,4))
    #pspconstraint4
    for i in SCHOOL:
        for j in VILLAGE:
            psp += X[(i,j)] <= a[i][j]*vil["POP"][j-1]   
    
    #pspconstraint6
    for i in SCHOOL:
        for j in VILLAGE:
            psp += X[(i,j)] - plusS[(i,j)] +MinusS[(i,j)] == Current_stat[i][j-1] 
    #pspconstraint7
    for i in SCHOOL:
        psp += U[i] <= alpha*(schl["Qi"][i-1]*Ync[i].varValue + p.lpSum(schl["Qi"+ str(l)][i-1]*Y[(i,l)].varValue for l in range(1,4)))
    #pspconstraint8
    for j in VILLAGE:
        psp += V[j]<= beta*vil["POP"][j-1]
    #pspconstraint9
    for j in VILLAGE:
        psp += p.lpSum(Costtransport[i][j]*X[(i,j)] for i in SCHOOL)-vil["POP"][j-1]*rho <= 0
    #pspsolve
    psp.solve()
    #dsp problem
    dsp = p.LpProblem("dsp",p.LpMaximize)
    #dsp variables
    A =p.LpVariable.dicts("A",VILLAGE)
    B = p.LpVariable.dicts("B",SCHOOL)
    C = p.LpVariable.dicts("C",[(i,j) for i in SCHOOL for j in VILLAGE],upBound =0)
    D =p.LpVariable.dicts("D",[(i,j) for i in SCHOOL for j in VILLAGE],lowBound=0)
    E = p.LpVariable.dicts("E",SCHOOL,upBound =0)
    F = p.LpVariable.dicts("F",VILLAGE,upBound =0)
    G = p.LpVariable.dicts("G",VILLAGE,upBound =0)
    #dsp objective
    dsp += p.lpSum(vil["POP"][j-1]*A[j] for j in VILLAGE) + p.lpSum(B[i]*schl["Qi"][i-1]*Ync[i].varValue+ p.lpSum(schl["Qi"+ str(i)][l-1]*Y[(i,l)].varValue for i in SCHOOL) for l in range(1,4)) + p.lpSum(C[(i,j)]*a[i][j]*vil["POP"][j-1] +D[(i,j)]*Current_stat[i][j-1] for i in SCHOOL for j in VILLAGE) + p.lpSum(E[i]*alpha*schl["Qi"][i-1]*Ync[i].varValue + p.lpSum(E[i]*alpha*schl["Qi"+ str(l)][i-1]*Y[(i,l)].varValue for l in range(1,4)) for i in SCHOOL) + p.lpSum(F[j]*beta*vil["POP"][j-1] for j in VILLAGE)
    
    #constraint1
    for i in SCHOOL:
        for j in VILLAGE:
            dsp += A[j] + B[i] +C[(i,j)]+ D[(i,j)] + G[j]*Costtransport[i][j] <= lambda1*Costtransport[i][j]
    #constraint2
    for i in SCHOOL:
        for j in VILLAGE:
            dsp += D[(i,j)] <= pi*lambda2
    #constraint3
    dsp += -p.lpSum(G[j]*vil["POP"][j-1]  for j in VILLAGE) <= lambda3
    #constraint4
    for j in VILLAGE:
        dsp += A[j] +F[j] <=0
    #constarint5
    for i in SCHOOL:
        dsp += B[i]+E[i] <= 0
    
    dsp.writeLP("Dual",writeSOS=1, mip=1, max_length=100)
    dsp.solve()
    sol = p.value(dsp.objective) + p.value(master.objective)
    UB = sol    
    aval[iteration] ={}
    bval[iteration] ={}
    cval[iteration] ={}
    dval[iteration] ={}
    Eval[iteration] ={}
    fval[iteration] ={} 
    gval[iteration] ={}
    for j in VILLAGE:
        aval[iteration][j]= A[j].varValue
        fval[iteration][j] =F[j].varValue
        gval[iteration][j] = G[j].varValue
    for i in SCHOOL:
        bval[iteration][i] =B[i].varValue
        Eval[iteration][i] = E[i].varValue
    for i in SCHOOL:
        cval[iteration][i] ={}
        dval[iteration][i] ={}
        for j in VILLAGE:
            cval[iteration][i][j] = C[(i,j)].varValue
            dval[iteration][i][j] =D[(i,j)].varValue
     #cuts
    for itr in range(1,iteration):
         master += Phi>= p.lpSum(vil["POP"][j-1]*aval[itr][j] for j in VILLAGE) + p.lpSum(bval[itr][i]*(schl["Qi"][i-1]*Ync[i]+ p.lpSum(schl["Qi"+ str(l)][i-1]*Y[(i,l)] for l in range(1,4))) for i in SCHOOL) + p.lpSum(cval[itr][i][j]*a[i][j]*vil["POP"][j-1] +dval[itr][i][j]*Current_stat[i][j-1] for i in SCHOOL for j in VILLAGE) + p.lpSum(Eval[itr][i]*alpha*schl["Qi"][i-1]*Ync[i] + Eval[itr][i]*alpha*schl["Qi"+ str(l)][i-1]*Y[(i,l)] for l in range(1,4) for i in SCHOOL) + p.lpSum(fval[itr][j]*beta*vil["POP"][j-1] for j in VILLAGE)
           
            
        
            
        




'''for i in SCHOOL:
    temp =[]
    for j in VILLAGE:
        temp.append(X[(i,j)].varValue)
    Current_stat["New"+str(i)] = temp
    Current_stat.to_excel("data\Current_stat.xlsx",index=False)'''


# In[ ]:





# In[ ]:




