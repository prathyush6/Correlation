# @ author: Prathyush Sambaturu

from gurobipy import *
import sys
import time
import math
import matplotlib.pyplot as plt

start_time = time.time()

#create the model for correlation 
m = Model("Correlation")


T = []
#N = [ [638, 635, 653, 842, 861, 1004, 1041, 1014], [1373, 1490, 1715, 1632, 1819, 1850, 1896, 2085] ]
L = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
alpha = 1000000
betasq = 455000000000000000
#beta2sq = 30
fp = open(sys.argv[1], "r")
fq = open(sys.argv[2], "w")
fw = open(sys.argv[4], "w")
#dictionary for state and its number
statedict = {}
N = []
noofstates = 0
line = fp.readline()
line = fp.readline()
nmin = 100000000.0
nmax = 0.0
while line:
      col = line.strip().split(",")
      #print(line[14])
      st_name = col[1]
      week_no = col[3]
      population = float(col[14])
      if population < nmin:
         nmin = population
      if population > nmax:
         nmax = population
      if st_name not in statedict:
         statedict[st_name] = noofstates
         noofstates = noofstates+1
         N.append([])
         N[noofstates-1].append(population)
      else:
         l = statedict[st_name]
         N[l].append(population)
      line = fp.readline()

#normalize values to between 1 to 100
for i in range(0, len(N)):
    for t in range(0, len(N[i])):
        N[i][t] = float(9.0 * (N[i][t] - nmin)/(nmax -nmin)+ 1)
        #print(str(N[i][t]))    

#print(statedict)
#fq.write(str(N)+"\n")
fq.write("|V| = "+str(len(N))+"\n")
fq.write("|T| = "+str(len(N[0]))+"\n")
fp.close()


V = []
noV = len(N)
for i in range(0, noV):
    V.append(i)
fq.write(str(V))

noT = len(N[0])
for i in range(0, noT):
    T.append(i)
#roots of the cluster
#rA = 0
#rB = 1

yA = {}
yB = {}
#add variables: yA[i] (resp. yB[i]) is 1 if node i is in cluster A (resp. B)
for i in V:
    yA[i] = m.addVar(vtype = GRB.BINARY, name = "yA("+str(i)+")")
    yB[i] = m.addVar(vtype = GRB.BINARY, name = "yB("+str(i)+")")

#root constraints
#m.addConstr( yA[rA] == 1, name = "RC1")
#m.addConstr( yA[rB] == 0, name = "RC2")
#m.addConstr( yB[rA] == 0, name = "RC3")
#m.addConstr( yB[rB] == 1, name = "RC4")
for i in V:
    m.addConstr(yA[i] + yB[i] <= 1, name = "RC39["+str(i)+"]")

#yA = {}
#yB = {}
#yA[0]/ = 1
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
    m.addConstr( XA[t] >= quicksum(yA[i] * N[i][t] for i in V), name = "RC5["+str(t)+"]")
    m.addConstr( XA[t] - 1<= quicksum(yA[i] * N[i][t] for i in V), name = "RC51["+str(t)+"]")
    m.addConstr( XB[t] >= quicksum(yB[i] * N[i][t] for i in V), name = "RC6["+str(t)+"]")
    m.addConstr( XB[t] - 1<= quicksum(yB[i] * N[i][t] for i in V), name = "RC61["+str(t)+"]")

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
DA = m.addVar(lb = 0, name = "DA")
m.addConstr(DA == len(T)* (quicksum(Z3A[l, ll, t]*pow(2, l+ll) for l in L for ll in L for t in T)) - quicksum(Z4A[l, ll]*pow(2, l+ll) for l in L for ll in L), name = "RC26")

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
DB = m.addVar(lb = 0, name = "DB")
m.addConstr(DB == len(T)* (quicksum(Z3B[l, ll, t]*pow(2, l+ll) for l in L for ll in L for t in T)) - quicksum(Z4B[l, ll]*pow(2, l+ll) for l in L for ll in L), name = "RC31")

Z5A = {}
Z5B = {}
#add variables z5A and z5B
for l in L:
    Z5A[l] = m.addVar(vtype = GRB.BINARY, name = "Z5A["+str(l)+"]")
    Z5B[l] = m.addVar(vtype = GRB.BINARY, name = "Z5B["+str(l)+"]")

m.addConstr( DA - 1 <= quicksum(Z5A[l] for l in L), name = "RC32")
m.addConstr( DA >= quicksum(Z5B[l] for l in L), name = "RC33")
m.addConstr( DB - 1 <= quicksum(Z5B[l] for l in L), name = "RC34")
m.addConstr( DB >= quicksum(Z5B[l] for l in L), name = "RC35")

Z5 = {}
for l in L:
    for ll in L:
        Z5[l, ll] = m.addVar(vtype = GRB.BINARY, name = "Z5["+str(l)+","+str(ll)+"]")

for l in L:
    for ll in L:
        m.addConstr( Z5[l, ll] <= Z5A[l], name = "RC36["+str(l)+","+str(ll)+"]")
        m.addConstr( Z5[l, ll] <= Z5B[ll], name = "RC37["+str(l)+","+str(ll)+"]")

m.addConstr( quicksum(Z5[l, ll]*pow(2, l+ll) for l in L for ll in L) <= betasq, name = "RC38")           

m.addConstr( quicksum(yA[i] + yB[i] for i in V) <= 6, name = "RC40")
#objective
m.setObjective( quicksum(1.0 * Z1[l, ll, t] for l in L for ll in L for t in T) + quicksum(1.0 * Z2[l, ll] for l in L for ll in L) + quicksum(1.0 * Z3A[l, ll, t] for l in L for ll in L for t in T) + quicksum(1.0 * Z4A[l, ll] for l in L for ll in L) + quicksum(1.0 * Z3B[l, ll, t] for l in L for ll in L for t in T) + quicksum(1.0 * Z4B[l, ll] for l in L for ll in L) + quicksum(1.0 * Z5A[l] for l in L) + quicksum(Z5B[l] for l in L) + quicksum(1.0 * Z5[l, ll] for l in L for ll in L), GRB.MAXIMIZE) 

m.update()
m.write('correlation.lp')

gap = float(sys.argv[3])
fq.write("gap "+str(gap))
m.setParam(GRB.Param.MIPGap, gap)
m.optimize()

term1 = 0.0
for l in L:
    for ll in L:
        for t in T:
            term1 = term1 +  Z1[l, ll, t].X * pow(2, l+ll)
term1 = len(T)* term1
fq.write(str(term1)+"\n")

term2 = 0.0
for l in L:
    for ll in L:
        term2 = term2 + Z2[l, ll].X * pow(2, l+ll)
fq.write(str(term2)+"\n")

numerator = term1 - term2
fq.write("Numerator is "+str(numerator)+", while its guess was "+str(alpha)+"\n")

denominator = 0.0
for l in L:
    for ll in L:
        denominator = denominator + Z5[l, ll].X * pow(2, l+ll) 
Cluster1 = []
Cluster2 = []
for i in V:
    v1 = yA[i].X
    v2 = yB[i].X
    if v1 == 1:
       Cluster1.append(i)
    elif v2 == 1:
       Cluster2.append(i)

fq.write("Cluster 1: "+str(Cluster1)+"\n")
fq.write("Cluster 2: "+str(Cluster2)+"\n")


num1 = 0.0
num2 = 0.0
num3 = 0.0
den4 = 0.0
den5 = 0.0
for t in T:
    num1 = num1 + XA[t].X * XB[t].X
    num2 = num2 + XA[t].X
    num3 = num3 + XB[t].X
    den4 = den4 + XA[t].X * XA[t].X
    den5 = den5 + XB[t].X * XB[t].X
    
numerator1 = len(T) * num1 - num2 * num3
denominator1 = math.sqrt( len(T) * den4 - num2 * num2) * math.sqrt( len(T) * den5 - num3 * num3) 
correlation1  = numerator1/denominator1


fq.write("Denominator is "+str(denominator)+", while its guess was "+str(betasq)+"\n")
den_sqrt = math.sqrt(denominator)
fq.write(str(den_sqrt)+"\n")
#corr = numerator/den_sqrt
#fq.write("correlation by program "+str(corr)+"\n")
fq.write("Correlation1 "+str(correlation1)+"\n")

clusvals1 = []
clusvals2 = []
for t in T:
    clusvals1.append(XA[t].X)
    clusvals2.append(XB[t].X)
    fw.write(str(t)+","+str(XA[t].X)+","+str(XB[t].X)+"\n")


fq.write("Clusters time series\n")
fq.write(str(clusvals1))
fq.write("\n")
fq.write(str(clusvals2))

elapsed_time = time.time() - start_time
fq.write("Total time elapsed: "+str(elapsed_time)+"\n")

fq.close()
fw.close()

