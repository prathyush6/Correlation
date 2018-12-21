# @ author: Prathyush Sambaturu

from gurobipy import *
import sys
import time

start_time = time.time()

#create the model for correlation 
m = Model("Correlation")

V = [0, 1]
T = [0, 1, 2, 3, 4, 5, 6, 7]
N = [ [638, 635, 653, 842, 861, 1004, 1041, 1014], [1373, 1490, 1715, 1632, 1819, 1850, 1896, 2085] ]
L = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
alpha = 1900000
betasq = 5500000000000
#beta2sq = 30

#roots of the cluster
rA = 0
rB = 1

yA = {}
yB = {}
#add variables: yA[i] (resp. yB[i]) is 1 if node i is in cluster A (resp. B)
for i in V:
    yA[i] = m.addVar(vtype = GRB.BINARY, name = "yA("+str(i)+")")
    yB[i] = m.addVar(vtype = GRB.BINARY, name = "yB("+str(i)+")")

#root constraints
m.addConstr( yA[rA] == 1, name = "RC1")
m.addConstr( yA[rB] == 0, name = "RC2")
m.addConstr( yB[rA] == 0, name = "RC3")
m.addConstr( yB[rB] == 1, name = "RC4")
#yA = {}
#yB = {}
#yA[0] = 1
#yA[1] = 0
#yB[0] = 0
#yB[1] = 1 

XA = {}
XB = {}
#add variables XA and XB
for t in T:
    XA[t] = m.addVar(lb = 0.0, name = "XA["+str(t)+"]")
    XB[t] = m.addVar(lb = 0.0, name = "XB["+str(t)+"]")

for t in T:
    m.addConstr( XA[t] == quicksum(yA[i] * N[i][t] for i in V), name = "RC5["+str(t)+"]")
    m.addConstr( XB[t] == quicksum(yB[i] * N[i][t] for i in V), name = "RC6["+str(t)+"]")

Z1A = {}
Z1B = {}
#add variables Z1A and Z1B
for l in L:
    for t in T:
        Z1A[l, t] = m.addVar(vtype = GRB.BINARY, name = "Z1A["+str(l)+","+str(t)+"]")
        Z1B[l, t] = m.addVar(vtype = GRB.BINARY, name = "Z1B["+str(l)+","+str(t)+"]")

for t in T:
    m.addConstr( XA[t] -1 <= quicksum(Z1A[l, t]* pow(2, l) for l in L), name = "RC7["+str(t)+"]")
    m.addConstr( XA[t] >= quicksum(Z1A[l, t] * pow(2, l) for l in L), name = "RC8["+str(t)+"]")
    m.addConstr( XB[t] -1 <= quicksum(Z1B[l, t] * pow(2, l) for l in L), name = "RC9["+str(t)+"]")
    m.addConstr( XB[t] >= quicksum(Z1B[l, t] * pow(2, l) for l in L), name = "RC10["+str(t)+"]")

Z1 = {}
#add variables Z
for l in L:
    for ll in L:
        for t in T:
            Z1[l, ll, t] = m.addVar(vtype = GRB.BINARY, name = "Z1["+str(l)+","+str(ll)+","+str(t)+"]")
            m.addConstr( Z1[l, ll, t] <= Z1A[l, t], name = "RC11["+str(l)+","+str(ll)+","+str(t)+"]")
            m.addConstr( Z1[l, ll, t] <= Z1B[ll, t], name = "RC12["+str(l)+","+str(ll)+","+str(t)+"]")

#Setup SA and SB variables
SA = m.addVar(lb = 0, name = "SA")
SB = m.addVar(lb = 0, name = "SB")

m.addConstr( SA == quicksum(XA[t] for t in T), name = "RC13["+str(t)+"]")
m.addConstr( SB == quicksum(XB[t] for t in T), name = "RC14["+str(t)+"]")

Z2A = {}
Z2B = {}
#add variables Z2A and Z2B
for l in L:
    Z2A[l] = m.addVar(vtype = GRB.BINARY, name = "Z2A["+str(l)+"]")
    Z2B[l] = m.addVar(vtype = GRB.BINARY, name = "Z2B["+str(l)+"]")

m.addConstr( SA - 1 <= quicksum(Z2A[l]*pow(2, l) for l in L), name = "RC15")
m.addConstr( SA >= quicksum(Z2A[l]*pow(2, l) for l in L), name = "RC16")
m.addConstr( SB - 1 <= quicksum(Z2B[l]*pow(2, l) for l in L), name = "RC17")
m.addConstr( SB >= quicksum(Z2B[l]*pow(2, l) for l in L), name = "RC18")

Z2 = {}
#add variables Z2
for l in L:
    for ll in L:
        Z2[l, ll] = m.addVar(vtype = GRB.BINARY, name = "Z2["+str(l)+","+str(ll)+"]")
        m.addConstr( Z2[l, ll] <= Z2A[l], name = "RC19["+str(l)+","+str(ll)+"]")
        m.addConstr( Z2[l, ll] <= Z2B[ll], name = "RC20["+str(l)+","+str(ll)+"]")

#Numerator constraint
m.addConstr( len(T)* (quicksum(Z1[l, ll, t]*pow(2, l+ll) for l in L for ll in L for t in T)) - quicksum(Z2[l, ll]*pow(2, l+ll) for l in L for ll in L)  >= alpha, name = "RC21") 

Z3A = {}
Z4A = {}
#add variables Z3A and Z4A
for l in L:
    for ll in L:
        for t in T:
            Z3A[l, ll, t] = m.addVar(vtype = GRB.BINARY, name = "Z3A["+str(l)+","+str(ll)+","+str(t)+"]")
            m.addConstr( Z3A[l, ll, t] <= Z1A[l, t], name = "RC22["+str(l)+","+str(ll)+","+str(t)+"]")
            m.addConstr( Z3A[l, ll, t] <= Z1A[ll, t], name = "RC23["+str(l)+","+str(ll)+","+str(t)+"]")

for l in L:
    for ll in L:
        Z4A[l, ll] = m.addVar(vtype = GRB.BINARY, name = "Z4A["+str(l)+","+str(ll)+"]")
        m.addConstr( Z4A[l, ll] <= Z2A[l], name = "RC24["+str(l)+","+str(ll)+"]")
        m.addConstr( Z4A[l, ll] <= Z2A[ll], name = "RC25["+str(l)+","+str(ll)+"]")


#Denominator1 constraint
m.addConstr(len(T)* (quicksum(Z3A[l, ll, t]*pow(2, l+ll) for l in L for ll in L for t in T)) - quicksum(Z4A[l, ll]*pow(2, l+ll) for l in L for ll in L) <= beta1sq, name = "RC26")

Z3B = {}
Z4B = {}
#add variables Z3B and Z4B


for l in L:
    for ll in L:
        for t in T:
            Z3B[l, ll, t] = m.addVar(vtype = GRB.BINARY, name = "Z3B["+str(l)+","+str(ll)+","+str(t)+"]")
            m.addConstr( Z3B[l, ll, t] <= Z1B[l, t], name = "RC27["+str(l)+","+str(ll)+","+str(t)+"]")
            m.addConstr( Z3B[l, ll, t] <= Z1B[ll, t], name = "RC28["+str(l)+","+str(ll)+","+str(t)+"]")

for l in L:
    for ll in L:
        Z4B[l, ll] = m.addVar(vtype = GRB.BINARY, name = "Z4B["+str(l)+","+str(ll)+"]")
        m.addConstr( Z4B[l, ll] <= Z2B[l], name = "RC29["+str(l)+","+str(ll)+"]")
        m.addConstr( Z4B[l, ll] <= Z2B[ll], name = "RC30["+str(l)+","+str(ll)+"]")

#Denominator2 constraint
m.addConstr(len(T)* (quicksum(Z3B[l, ll, t]*pow(2, l+ll) for l in L for ll in L for t in T)) - quicksum(Z4B[l, ll]*pow(2, l+ll) for l in L for ll in L) <= beta2sq, name = "RC31")

#objective
m.setObjective( quicksum(Z1[l, ll, t] for l in L for ll in L for t in T) + quicksum(Z2[l, ll] for l in L for ll in L) + quicksum(Z3A[l, ll, t] for l in L for ll in L for t in T) + quicksum(Z4A[l, ll] for l in L for ll in L) + quicksum(Z3B[l, ll, t] for l in L for ll in L for t in T) + quicksum(Z4B[l, ll] for l in L for ll in L), GRB.MAXIMIZE)


m.update()
m.write('correlation.lp')


m.optimize()

elapsed_time = time.time() - start_time
print("Total time elapsed: "+str(elapsed_time))
