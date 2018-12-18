import sys

fp = open(sys.argv[1], "r")

#dictionary for state and its number
statedict = {}
N = [] 
noofstates = 0
line = fp.readline()
line = fp.readline() 
while line:
      col = line.strip().split(",")
      #print(line[14])
      st_name = col[1]
      week_no = col[3]
      population = int(col[14])
      if st_name not in statedict:
         statedict[st_name] = noofstates
         noofstates = noofstates+1
         N.append([])
         N[noofstates-1].append(population)
      else: 
         l = statedict[st_name]
         N[l].append(population)
      line = fp.readline()

#print(statedict)
print(N)
print(len(N))
fp.close()

