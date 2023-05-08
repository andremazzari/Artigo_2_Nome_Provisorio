import gurobipy as gp
import numpy as np
import itertools
import os


'''Optimization
Description: This function uses gurobi to make a linear program to optimize the expression of a specific inequality of the Lnd polytope.
It is imposed tthe following constraints in the behaviour:
- Non-signaling between Alice and Bob
- Non-disturbing in Bob's results.
- The behaviour must satisfy all facets of the L_G (Local with general response functions) polytope.

Output:
It prints the value of the optimization and the behaviour obtained.
It also calls an function the verify of the obtained behaviour satisfy all the contraints, and prints if there is any error in it.
'''
def Optimization(Lndinequality):
    m = gp.Model()
    m.Params.LogToConsole = 0

    #VARIABLES
    #Add behaviour variables
    behaviour_vars = gp.tuplelist([])
    for A in [0, 1]:
        for B in [[0, 1], [1, 2]]:
            for a in [0, 1]:
                for b in itertools.product([0,1], repeat = 2):
                    behaviour_vars.append((a,b[0],b[1],A,B[0],B[1]))
    p = m.addVars(behaviour_vars, name = 'p')

    #CONSTRAINTS
    #Add behaviour normalization constraints
    BehaviourNormalization = []
    for A in [0, 1]:
        for B in [[0, 1], [1, 2]]:
            Constraint = 0
            for a in [0, 1]:
                for b in itertools.product([0,1], repeat = 2):
                    Constraint += p[a,b[0],b[1],A,B[0],B[1]]
            BehaviourNormalization.append(Constraint)
    m.addConstrs((BehaviourNormalization[m] == 1 for m in range(len(BehaviourNormalization))), name = "BehaviourNormalization")

    #Non-signaling constraints
    NonSignalingConstraints = []
    #For Alice
    for A in [0,1]:
        for B in [[0, 1]]:
            #Next context
            Bnext = [B[0] + 1, B[1] + 1]
            
            for a in [0,1]:
                Constraint  = 0
                for b in itertools.product([0,1], repeat = 2):
                    Constraint += p[a,b[0],b[1],A,B[0],B[1]]
                for b in itertools.product([0,1], repeat = 2):
                    Constraint -= p[a,b[0],b[1],A,Bnext[0],Bnext[1]]
                NonSignalingConstraints.append(Constraint)
    #For Bob
    for B in [[0, 1], [1, 2]]:
        for b in itertools.product([0,1], repeat = 2):
            Constraint = 0
            for a in [0,1]:
                Constraint += p[a,b[0],b[1],0,B[0],B[1]]
            for a in [0,1]:
                Constraint -= p[a,b[0],b[1],1,B[0],B[1]]
            NonSignalingConstraints.append(Constraint)
    m.addConstrs((NonSignalingConstraints[m] == 0 for m in range(len(NonSignalingConstraints))), name = "NonSignalingConstraints")

    #Non-disturbing constraints
    NonDisturbingContraints = []
    for A in [0,1]:
        for B in [1]: #Only B1 is in different contexts
            for b in [0,1]:
                Constraint = 0

                for a in [0,1]:
                    for b2 in [0,1]:
                        Constraint += p[a, b, b2, A, B, 2]

                for a in [0,1]:
                    for b0 in [0,1]:
                        Constraint -= p[a, b0, b, A, 0, B]
                
                NonDisturbingContraints.append(Constraint)
    m.addConstrs((NonDisturbingContraints[m] == 0 for m in range(len(NonDisturbingContraints))), name = "NonDisturbingConstraints")

    #Local disturbing contraints
    path = os.path.dirname(__file__)
    file = open(path + "/2V-L-Facets-Prob.txt", "r")
    inequalities = file.readlines()
    file.close()

    for inequality in inequalities:
        InequalityExpr = 0

        lhs, rhs = inequality.split(" <= ")
        lhs = lhs.strip()
        bound = int(rhs.strip())
        for term in lhs.split(" "):
            coefficient, info = term.split("p")
            
            if coefficient == '-':
                coefficient = -1
            elif coefficient == '+' or coefficient == '':
                coefficient = 1
            else:
                coefficient = int(coefficient)

            results, measurements = info.split("|")
            a = int(results[0])
            b0 = int(results[1])
            b1 = int(results[2])
            A = int(measurements[1])
            B0 = int(measurements[3])
            B1 = int(measurements[5])

            InequalityExpr += coefficient*p[a, b0, b1, A, B0, B1]
        m.addConstr(InequalityExpr <= bound)

    #Maximize Lnd inequality
    obj = 0

    lhs, rhs = Lndinequality.split(" <= ")
    lhs = lhs.strip()
    bound = int(rhs.strip())
    for term in lhs.split(" "):
        coefficient, info = term.split("p")
        
        if coefficient == '-':
            coefficient = -1
        elif coefficient == '+' or coefficient == '':
            coefficient = 1
        else:
            coefficient = int(coefficient)

        results, measurements = info.split("|")
        a = int(results[0])
        b0 = int(results[1])
        b1 = int(results[2])
        A = int(measurements[1])
        B0 = int(measurements[3])
        B1 = int(measurements[5])

        if B0 == 0 and B1 == 3:
            B0 = 3
            B1 = 0
            btemp = b0
            b0 = b1
            b1 = btemp

        obj += coefficient*p[a, b0, b1, A, B0, B1]
    m.setObjective(obj, gp.GRB.MAXIMIZE)
    m.optimize()

        
    print("Inequality bound: ", bound)
    print("Optimization value: ", m.ObjVal)

    #Get behaviour
    behaviour = {}
    for A in [0, 1]:
        for B in [[0, 1], [1, 2]]:
            for a in [0, 1]:
                for b in itertools.product([0,1], repeat = 2):
                    behaviour[a,b[0],b[1],A,B[0],B[1]] = p[a,b[0],b[1],A,B[0],B[1]].X
    print(behaviour)
    print("Error in behaviour: ",Verify_Behaviour(behaviour))

    return behaviour, bound, m.ObjVal


'''Verify_Behaviour
Description: verify if a behaviour satisfy the non-signaling and non-disturbing conditions, and all the facets of the LD polytope.
Input: behaviour in the form of a dictionary.
Output:
 - False if there is no error.
 - True if there is some error.
'''
def Verify_Behaviour(p):
    Error = False

    #Verify non-signaling
    #For Alice
    for A in [0,1]:
        for B in [[0, 1]]:
            #Next context
            Bnext = [B[0] + 1, B[1] + 1]
            
            for a in [0,1]:
                Expr1 = 0
                for b in itertools.product([0,1], repeat = 2):
                    Expr1 += p[a,b[0],b[1],A,B[0],B[1]]
                Expr2 = 0
                for b in itertools.product([0,1], repeat = 2):
                    Expr2 += p[a,b[0],b[1],A,Bnext[0],Bnext[1]]
                if abs(Expr1 - Expr2) > 1e-7:
                    print("Error in non-signaling constraint", Expr1, Expr2)
                    Error = True
    #For Bob
    for B in [[0, 1], [1, 2]]:
        for b in itertools.product([0,1], repeat = 2):
            Expr1 = 0
            for a in [0,1]:
                Expr1 += p[a,b[0],b[1],0,B[0],B[1]]
            Expr2 = 0
            for a in [0,1]:
                Expr2 += p[a,b[0],b[1],1,B[0],B[1]]
            if abs(Expr1 - Expr2) > 1e-7:
                    print("Error in non-signaling constraint", Expr1, Expr2)
                    Error = True


    #Verify non-disturbance
    for A in [0,1]:
        for B in [1]: #Only B1 is in different contexts
            for b in [0,1]:
                Expr1 = 0
                for b2 in [0,1]:
                    for a in [0,1]:
                        Expr1 += p[a, b, b2, A, B, 2]

                Expr2 = 0
                for b0 in [0,1]:
                    for a in [0,1]:
                        Expr2 += p[a, b0, b, A, 0, B]
                
                if abs(Expr1 - Expr2) > 1e-5:
                    print("Error in non-disturbing constraint", Expr1, Expr2)
                    Error = True
    
    #Local disturbing contraints
    path = os.path.dirname(__file__)
    file = open(path + "/2V-L-Facets-Prob.txt", "r")
    inequalities = file.readlines()
    file.close()

    i = 0
    for inequality in inequalities:
        InequalityExpr = 0

        lhs, rhs = inequality.split(" <= ")
        lhs = lhs.strip()
        bound = int(rhs.strip())
        for term in lhs.split(" "):
            coefficient, info = term.split("p")
            
            if coefficient == '-':
                coefficient = -1
            elif coefficient == '+' or coefficient == '':
                coefficient = 1
            else:
                coefficient = int(coefficient)

            results, measurements = info.split("|")
            a = int(results[0])
            b0 = int(results[1])
            b1 = int(results[2])
            A = int(measurements[1])
            B0 = int(measurements[3])
            B1 = int(measurements[5])

            InequalityExpr += coefficient*p[a, b0, b1, A, B0, B1]
        
        if (InequalityExpr > bound):
            print("Error in local disturbing constraint", InequalityExpr, bound)
            Error = True
    return Error

'''Lnd_Vertices_2V_Prob
Returns the extreme vertices of the Lnd polytope for the 2V scenario
'''
def Lnd_Vertices_2V_Prob():
    first = True

    #deterministic vertices
    for a_results in itertools.product([0,1], repeat = 2): #Deterministic results for Alice.
        for b_results in itertools.product([0,1], repeat = 3): #Deterministic results for Bob.
            vertice = []
            for A in [0,1]: #Alice's measurements
                for B in [[0,1], [1,2]]: #Bob's contexts
                    for a in [0,1]:
                        for b in itertools.product([0,1], repeat = 2):
                            if a_results[A] == a and b_results[B[0]] == b[0] and b_results[B[1]] == b[1]:
                                vertice.append(1)
                            else:
                                vertice.append(0)
            if first == True:
                Vertices = np.array([vertice], int)
                first = False
            else:
                Vertices = np.concatenate((Vertices, np.array([vertice], int)), axis = 0)
    return Vertices

'''Lnd_Vertices_2V_Prob
Returns the extreme vertices of the L_G polytope for the 2V scenario
'''
def L_Vertices_2V_Prob():
    first = True

    #Deterministic vertices (for each context)
    for a_results in itertools.product([0,1], repeat = 2): #Deterministic results for Alice's measurements.
        for b_results in itertools.product([[0,0], [0,1], [1,0], [1,1]], repeat = 2): #Deterministic results for Bob's contexts.
            vertice = []
            for A in [0,1]: #Alice's measurements
                for B in [0, 1]: #Bob's contexts
                    for a in [0,1]:
                        for b in itertools.product([0,1], repeat = 2):
                            if a_results[A] == a and b_results[B][0] == b[0] and b_results[B][1] == b[1]:
                                vertice.append(1)
                            else:
                                vertice.append(0)
            if first == True:
                Vertices = np.array([vertice], int)
                first = False
            else:
                Vertices = np.concatenate((Vertices, np.array([vertice], int)), axis = 0)
    return Vertices

'''Get_Lnd_Decomposition
Description: Use gurobi to try to find an Lnd decommposition for a behaviour in the 2-V scenario.
Output: 
    print status (status = 3 if infeasible)
    print vertices of the decomposition (if status != 3)
'''
def Get_Lnd_Decomposition(p):
    #Put behaviour in array form
    behaviour = []
    for A in [0, 1]:
        for B in [[0, 1], [1, 2]]:
            for a in [0, 1]:
                for b in itertools.product([0,1], repeat = 2):
                    behaviour.append(p[a,b[0],b[1],A,B[0],B[1]])
    behaviour = np.array(behaviour, float)

    m = gp.Model()
    m.Params.LogToConsole = 0

    #Get LD vertices
    ExtremeVertices = np.transpose(Lnd_Vertices_2V_Prob())

    #Coefficients variables
    l = m.addMVar(shape = ExtremeVertices.shape[1], name = "lambda")

    #Add normalization constraint
    m.addConstr(l.sum() == 1,  name = "Nomalization")

    #Add behaviour constraint
    m.addConstr(ExtremeVertices @ l == behaviour, name = "BehaviourConstraint")

    m.setObjective(1)
    m.optimize()

    print("Status: ", m.status)

    if m.status != 3:
        count = 0
        i = 0
        for entry in l.X:
            if entry > 0.0:
                print(entry)
                print(ExtremeVertices[:,i])
                count += 1
            i += 1
        print("Non zero entries: ", count) 


'''Get_L_Decomposition
Description: Use gurobi to try to find an L decommposition for a behaviour in the 2-V scenario.
Output: 
    print status (status = 3 if infeasible)
    print vertices of the decomposition
'''
def Get_L_Decomposition(p):
    #Put behaviour in array form
    behaviour = []
    for A in [0, 1]:
        for B in [[0, 1], [1, 2]]:
            for a in [0, 1]:
                for b in itertools.product([0,1], repeat = 2):
                    behaviour.append(p[a,b[0],b[1],A,B[0],B[1]])
    behaviour = np.array(behaviour, float)

    m = gp.Model()
    m.Params.LogToConsole = 0

    #Get LD vertices
    ExtremeVertices = np.transpose(L_Vertices_2V_Prob())

    #Coefficients variables
    l = m.addMVar(shape = ExtremeVertices.shape[1], name = "lambda")

    #Add normalization constraint
    m.addConstr(l.sum() == 1,  name = "Nomalization")

    #Add behaviour constraint
    m.addConstr(ExtremeVertices @ l == behaviour, name = "BehaviourConstraint")

    m.setObjective(1)
    m.optimize()

    print("Status: ", m.status)

    count = 0
    i = 0
    for entry in l.X:
        if entry > 0.0:
            print(entry)
            print(repr(ExtremeVertices[:,i]))
            count += 1
        i += 1
    print("Non zero entries: ", count)


def Vertice_To_Behaviour(vertice):
    p = {}
    i = 0
    for A in [0,1]:
        for B in [[0,1],[1,2]]:
            for a in [0,1]:
                for b in itertools.product([0,1], repeat = 2):
                    p[a,b[0],b[1],A,B[0],B[1]] = vertice[i]
                    i += 1
    return p

def Print_Behaviour_NonZero_Entries(p):
    for A in [0,1]:
        for B in [[0,1],[1,2]]:
            for a in [0,1]:
                for b in itertools.product([0,1], repeat = 2):
                    if p[a,b[0],b[1],A,B[0],B[1]] > 0:
                        print(a,b[0],b[1],A,B[0],B[1],": " ,p[a,b[0],b[1],A,B[0],B[1]])




#Behaviour and inequality present in the article

inequality = '-p110|A0B0B1 +p101|A0B1B2 +p010|A1B0B1 +p010|A1B1B2 +p011|A1B1B2 +p100|A1B1B2 +p110|A1B1B2 +p111|A1B1B2 <= 1'

behaviour, _, value = Optimization(inequality)
print(behaviour)
print("Inequality value: ", value)

print("Lnd decomposition:")
Get_Lnd_Decomposition(behaviour)
print("L decomposition:")
Get_L_Decomposition(behaviour)


#Results of the optimization
p = {(0, 0, 0, 0, 0, 1): 0.0, (0, 0, 1, 0, 0, 1): 0.5, (0, 1, 0, 0, 0, 1): 0.5, (0, 1, 1, 0, 0, 1): 0.0, (1, 0, 0, 0, 0, 1): 0.0, (1, 0, 1, 0, 0, 1): 0.0, (1, 1, 0, 0, 0, 1): 0.0, (1, 1, 1, 0, 0, 1): 0.0, (0, 0, 0, 0, 1, 2): 0.5, (0, 0, 1, 0, 1, 2): 0.0, (0, 1, 0, 0, 1, 2): 0.5, (0, 1, 1, 0, 1, 2): 0.0, (1, 0, 0, 0, 1, 2): 0.0, (1, 0, 1, 0, 1, 2): 0.0, (1, 1, 0, 0, 1, 2): 0.0, (1, 1, 1, 0, 1, 2): 0.0, (0, 0, 0, 1, 0, 1): 0.0, (0, 0, 1, 1, 0, 1): 0.0, (0, 1, 0, 1, 0, 1): 0.5, (0, 1, 1, 1, 0, 1): 0.0, (1, 0, 0, 1, 0, 1): 0.0, (1, 0, 1, 1, 0, 1): 0.5, (1, 1, 0, 1, 0, 1): 0.0, (1, 1, 1, 1, 0, 1): 0.0, (0, 0, 0, 1, 1, 2): 0.0, (0, 0, 1, 1, 1, 2): 0.0, (0, 1, 0, 1, 1, 2): 0.5, (0, 1, 1, 1, 1, 2): 0.0, (1, 0, 0, 1, 1, 2): 0.5, (1, 0, 1, 1, 1, 2): 0.0, (1, 1, 0, 1, 1, 2): 0.0, (1, 1, 1, 1, 1, 2): 0.0}

Vertices_Decomposition = []
Vertices_Decomposition.append(np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]))
Vertices_Decomposition.append(np.array([0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]))

for i in range(len(Vertices_Decomposition)):
    print("Vertice ",i + 1)
    Print_Behaviour_NonZero_Entries(Vertice_To_Behaviour(Vertices_Decomposition[i]))
    print("")
