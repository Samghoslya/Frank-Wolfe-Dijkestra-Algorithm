import numpy as np
import pandas as pd
from scipy import optimize
from scipy.optimize import minimize_scalar
import scipy.linalg as la
import numpy.linalg as la2
import scipy.integrate as integrate

def firstfunc(t0,xa,ca):
	ta = t0*(1+0.15*(xa/ca)**4) #calculating Cij values using Vij and Kij, constants 0.15 and 4 have been obtained from online sources
	return ta #xa represents Vij and ca represents Kij

def estimateZ(alpha,xa,ca,t0,ya):
	Z = 0
	for i in range(len(xa)):
		Z += integrate.quad(lambda x: firstfunc(t0[i],x,ca[i]),0,xa[i]+alpha*(ya[i]-xa[i]))[0]
	return Z

def linearsearch(xa,ca,t0,ya):
	alpha = minimize_scalar(estimateZ, args=(xa, ca, t0, ya), bounds = (0,1), method = 'Bounded') #Optimizes and finds the minimal value and returns alpha
	return alpha.x

#import the demand, free flow travel times, network and capacities
LinkNode = pd.read_csv("linknode.csv", header = None)
Q = pd.read_csv("demand.csv", header = None)
Q = Q.to_numpy()
y = np.copy(Q)
coeff = pd.read_csv("coefficient.csv", header = None)
coeff = coeff.to_numpy()
t0 = coeff[:,0]

   
n = 76 # number of total links
k = 24 # number of total nodes

	## creat link flow matrix for iterations (Y)
s = (n,k) # 76 links, 24 nodes/origins
	#Y = np.zeros(s) # each entry represents the flow on link a from origin i


t0 = coeff [:,0] # free flow travel time from 1st column of *Coeff*
ca = coeff[:,1] # capacity for each link

origq = np.sum(Q, axis = 1) # row sums of Q: total flow from origin i
destq = -Q
s2 = (k,k)
RHT = np.zeros(s2)# each row represents an origin, each column represents the flow on node k (with origin i )
RHT = -Q
np.fill_diagonal(RHT, origq)
	#print RHT

c0_0 = np.transpose(t0)
c_0 = np.tile(c0_0,k)

A0 = np.transpose(LinkNode) # Construct block matrix for A_eq
A1 = [A0]*k
A = la.block_diag(*A1)
	#print A

b0 = np.transpose(RHT)
b = np.ravel(b0, order = 'F')[:,np.newaxis] # construct long b_eq

ybounds = (0, None)
result = optimize.linprog(
c_0, A_eq = A, b_eq = b, bounds = (ybounds), options = {"disp":True, "maxiter":2000,"bland":True}
		)
	#print result
	#print len(result['x'])
result = np.reshape(result['x'],(k,n))
	#print result
xa = np.sum(result, axis = 0) # intialization xa
	#print xa
ta = firstfunc(t0,xa,ca)

step = 0
tanorm = 1000000
iteration = []
Z = []

while (tanorm>(n/10)):
	iteration.append(step)
	ta_old = ta
	c0 = np.transpose(ta)
	c = np.tile(c0,k)
	result = optimize.linprog(
	c, A_eq = A, b_eq = b, bounds = (ybounds), options = {"disp":True, "maxiter":2000,"bland":True}
		)
	resultx = np.reshape(result['x'],(k,n))
	ya = np.sum(resultx, axis = 0)
	alpha = linearsearch(xa,ca,t0,ya)
	print ("alpha is", alpha)
	xa = (1-alpha)*xa + alpha * ya
	print ("xa is ",xa)
	ta = firstfunc(t0,xa,ca)
	tanorm = la2.norm(ta-ta_old)
	z = np.dot(np.transpose(xa),ta)
	Z.append(z)
	print ("ta is ", ta)
	print ("norm of ta is ", tanorm)
	step +=1# allow each link has 0.1 diff. in ta on average
	
print("flow on each link is ", xa)
print("ttravel time on each link is", ta)

LinkNode = LinkNode.to_numpy()
taa = pd.DataFrame(ta)

def linktomatrix(LinkNode, ta):
    k=24
    l1 = []
    l2 = []
    for i in range (len(LinkNode)):
        for j in range(k):
    #            print("i=",i,"j=",j)
                if LinkNode[i][j] == 1:
                    
                    l1.append(j+1)
    #                print("l1", l1)
                    for n in range(k):
                        if LinkNode[i][n] == -1:
                            l2.append((n+1))
    #                        print("l2",l2)
                            

    k = pd.DataFrame({'l1':l1, 'l2':l2})

    od = pd.concat([k,ta],axis=1)
#    print("od is", od)
    #np.savetxt("od.csv", od, delimiter = ",")

    od.columns = ["Start", "end","cost"]
    od.head()

    def retro_dictify(df):
    	d = {}
    	for val in df.values:
    		key1 = str(int(val[0]))
    		if key1 not in list(d.keys()):
    			d[key1] = {}
    			key2 = str(int(val[1]))
    			d[key1][key2] = val[2]
    		else:
    			key2 = str(int(val[1]))
    			d[key1][key2] = val[2]
    	return d

    d = retro_dictify(od)
    #print(d)

    graph = d
    #print(d)
    def dijkestra(graph):
        w, h = 24, 24;
        Matrix = [[0 for x in range(w)] for y in range(h)]
        Matrix1 = [[0 for x in range(w)] for y in range(h)]
        for start in range(1,25):
            for goal in range(1,25):
                goal = str(goal)
                start = str(start)
                shortest_distance = {}
                predecessor = {}
                unseenNodes = graph.copy()
                infinity = 999999
                path = []
                for node in unseenNodes:
                    shortest_distance[node] = infinity
                shortest_distance[start] = 0
            
                while unseenNodes:
                    minNode = None
                    for node in unseenNodes:
                        if minNode is None:
                            minNode = node
                        elif shortest_distance[node] < shortest_distance[minNode]:
                            minNode = node
            
                    for childNode, weight in graph[minNode].items():
                        if weight + shortest_distance[minNode] < shortest_distance[childNode]:
                            shortest_distance[childNode] = weight + shortest_distance[minNode]
                            predecessor[childNode] = minNode
                    unseenNodes.pop(minNode)
    #            print(shortest_distance)
                
                currentNode = goal
                while currentNode != start:
                    try:
                        path.insert(0,int(currentNode))
                        currentNode = predecessor[currentNode]
                    except KeyError:
                        print('path not reachable')
                        break
                path.insert(0,start)
                if shortest_distance[goal] != infinity:
    #                print('shortest distance is ' + str(shortest_distance[goal]))
                    Matrix[int(start)-1][int(goal)-1] = shortest_distance[goal]
    #                print('and the path between '+ str(int(start)) + ' and ' + str(int(goal)) + ' is ' + str(path))
                    
                    currentNode = int(currentNode)
                    currentNode = int(goal) + 1
                    Matrix1[int(start)-1][int(goal)-1] = str(path)
        #print(Matrix)
        return Matrix, Matrix1
    time, routes = dijkestra(graph)
   
    return time, routes

t_time, Routes = linktomatrix(LinkNode, taa)

#Shortest route nodes to reach each destination from a origin
Routes = pd.DataFrame(Routes)
#Routes.to_csv("Routes.csv", index=False)

#Travel time for all the shortest route paths
Time = pd.DataFrame(t_time)
#Time.to_csv("Travel_time.csv", index=False)

#link specific flow
Link_flow = xa

#Link specific travel time
Link_time = ta
