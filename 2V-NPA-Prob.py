import numpy as np
import picos
import itertools
from time import time

t0 = time()


def powerset(iterable):
    #powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))


'''SortMeasurements
Description: receives an array with individual measurement strings, and implement a buble sort to sort the measurements.
Input:
    measurements: array of individual measurement strings.
Output:
    Sorted array with individual meassurement strings.
'''
def SortMeasurements(measurements):
    #bubble sort
    changed = True
    upperindex = len(measurements)
    tempupperindex = upperindex
    while changed == True:
        changed = False
        for j in range(upperindex - 1):
            TempMeasurement1 = measurements[j]
            TempMeasurement2 = measurements[j + 1]
            
            #Verify if there is a minus in the beggining of the correlator.
            if TempMeasurement1[0] == '-':
                TempMeasurement1 = TempMeasurement1[1:]
            if TempMeasurement2[0] == '-':
                TempMeasurement2 = TempMeasurement2[1:]
                
            if ord(TempMeasurement1[0]) > ord(TempMeasurement2[0]):
                #exchange measurements of different parties
                TempMeasurement = measurements[j + 1]
                measurements[j + 1] = measurements[j]
                measurements[j] = TempMeasurement
                tempupperindex = j + 1
                changed = True
            elif ord(TempMeasurement1[0]) == ord(TempMeasurement2[0]):
                TempMeasurement1 = TempMeasurement1[1:]
                TempMeasurement2 = TempMeasurement2[1:]
                if int(TempMeasurement1) > int(TempMeasurement2):
                    #exchange measurements in the same party
                    TempMeasurement = measurements[j + 1]
                    measurements[j + 1] = measurements[j]
                    measurements[j] = TempMeasurement
                    tempupperindex = j + 1
                    changed = True
        upperindex = tempupperindex
    
    return measurements

def Build_Dict_Probability_Indexes_2V():
    dict_tuple_to_index = {}
    dict_index_to_tuple = {}

    index = 0
    for A in [0, 1]:
        for B in [[0, 1], [1, 2]]:
            for a in [0, 1]:
                for b in itertools.product([0,1], repeat = 2):
                    dict_tuple_to_index[(a,b[0],b[1],A,B[0],B[1])] = index
                    dict_index_to_tuple[index] = (a,b[0],b[1],A,B[0],B[1])
                    index += 1
    
    return dict_tuple_to_index, dict_index_to_tuple


def NPA_2V_Q1(pA, pB, pAB, pBB):
    mA = 2 #Number of measurements for Alice
    ra = 2 #Number of results for Alice's measurements
    mB = 3 #Number of measurements for Bob 
    rb = 2 #Number of results of Bob's measurements
    N = 1 + (ra - 1)*mA + (rb - 1)*mB #Dimension of moment matrix (discards the last result of each measurement)
    
    P = picos.Problem()
    
    #Moment matrix
    M = picos.HermitianVariable('M', (N,N))

    l = picos.RealVariable('l')

    #constraints
    P.add_constraint(M[0,0] == 1)
    P.add_constraint(M[1,1] == pA[0,0])
    P.add_constraint(M[2,2] == pA[0,1])
    P.add_constraint(M[3,3] == pB[0,0])
    P.add_constraint(M[4,4] == pB[0,1])
    P.add_constraint(M[5,5] == pB[0,2])

    P.add_constraint(M[0,1] == pA[0,0])
    P.add_constraint(M[0,2] == pA[0,1])
    P.add_constraint(M[0,3] == pB[0,0])
    P.add_constraint(M[0,4] == pB[0,1])
    P.add_constraint(M[0,5] == pB[0,2])

    P.add_constraint(M[1,3] == pAB[0,0,0,0])
    P.add_constraint(M[1,4] == pAB[0,0,0,1])
    P.add_constraint(M[1,5] == pAB[0,0,0,2])

    P.add_constraint(M[2,3] == pAB[0,0,1,0])
    P.add_constraint(M[2,4] == pAB[0,0,1,1])
    P.add_constraint(M[2,5] == pAB[0,0,1,2])

    P.add_constraint(M[3,4] == pBB[0,0,0,1])

    P.add_constraint(M[4,5] == pBB[0,0,1,2])
    P.add_constraint(M - l*picos.I(N) >> 0)
    P.set_objective('max', l)

    S = P.solve()
    print(P.status)

    S.apply()
    print(l)
    #print(M)


'''NPA_2V
Builds NPA hierarchy for the 2V scenario
'''
def NPA_2V(p, level):
    mA = 2 #Number of measurements for Alice
    ra = 2 #Number of results for Alice's measurements
    mB = 3 #Number of measurements for Bob 
    rb = 2 #Number of results of Bob's measurements

    numA = (ra - 1)*mA #Number of projectors for Alice
    numB = (rb - 1)*mB #Number of projectors for Bob

    ContextObj = Contexts([['A0','B0','B1'], ['A0','B1','B2'], ['A1','B0','B1'], ['A1','B1','B2']])

    MatrixIndexes = MomentMatrix(mA, mB, numA, numB, level, ContextObj)
    
    P = picos.Problem()

    #Moment matrix
    M = picos.HermitianVariable('M', (MatrixIndexes.N,MatrixIndexes.N))

    l = picos.RealVariable('l')
    
    #constraints
    for i in range(MatrixIndexes.N):
        for j in range(i, MatrixIndexes.N):
            #print(i,j)
            constraint = MatrixIndexes.VerifyConstraint(i, j)
            if constraint != None:
                if constraint == '0':
                    P.add_constraint(M[i,j] == 0)
                elif constraint == '1':
                    P.add_constraint(M[i,j] == 1)
                else:
                    results, measurements = constraint.split('|')
                    prob_type = ''.join([measurements[k] for k in range(0,len(measurements),2)])
                    prob_index = tuple(int(results[k]) for k in range(len(results))) + tuple(int(measurements[k]) for k in range(1,len(measurements),2))
                    P.add_constraint(M[i,j] == p[prob_type][prob_index])
    
    
    P.add_constraint(M - l*picos.I(MatrixIndexes.N) >> 0)
    P.set_objective('max', l)

    S = P.solve()
    print(P.status)

    S.apply()
    print(l)


def prepare_behaviours_PICOS_constraint(p, prob_type, prob_index, dict_tuple_to_index):
    behaviour = {}

    if prob_type == 'A':
        a = prob_index[0]
        A = prob_index[1]

        pA = 0.0
        for b in itertools.product([0,1], repeat = 2):
            pA += p[dict_tuple_to_index[a,b[0],b[1],A,0,1]]
        return pA

    elif prob_type == 'B':
        b = prob_index[0]
        B = prob_index[1]

        pB = 0.0
        if B == 0 or B == 1:
            for a in [0,1]:
                for btemp in [0,1]:
                    pB += p[dict_tuple_to_index[a,b,btemp,0,B,B+1]]
        elif B == 2:
            for a in [0,1]:
                for btemp in [0,1]:
                    pB += p[dict_tuple_to_index[a,btemp,b,0,1,2]]
        return pB

    elif prob_type == 'AB':
        a = prob_index[0]
        b = prob_index[1]
        A = prob_index[2]
        B = prob_index[3]

        pAB = 0.0
        if B == 0 or B == 1:
            for btemp in [0,1]:
                pAB += p[dict_tuple_to_index[a,b,btemp,A,B,B+1]]
        elif B == 2:
            for btemp in [0,1]:
                pAB += p[dict_tuple_to_index[a,btemp,b,A,1,2]]
        return pAB

    elif prob_type == 'BB':
        b = [prob_index[0], prob_index[1]]
        B = [prob_index[2], prob_index[3]]

        pBB = 0.0
        for a in range(2):
            pBB += p[dict_tuple_to_index[a,b[0],b[1],0,B[0],B[1]]]
        return pBB

    elif prob_type == 'ABB':
        return p[dict_tuple_to_index[prob_index]]

def NPA_2V_Inequality(level, inequality):
    mA = 2 #Number of measurements for Alice
    ra = 2 #Number of results for Alice's measurements
    mB = 3 #Number of measurements for Bob 
    rb = 2 #Number of results of Bob's measurements

    numA = (ra - 1)*mA #Number of projectors for Alice
    numB = (rb - 1)*mB #Number of projectors for Bob

    ContextObj = Contexts([['A0','B0','B1'], ['A0','B1','B2'], ['A1','B0','B1'], ['A1','B1','B2']])

    MatrixIndexes = MomentMatrix(mA, mB, numA, numB, level, ContextObj)
    
    P = picos.Problem()

    #Moment matrix variable
    M = picos.HermitianVariable('M', (MatrixIndexes.N,MatrixIndexes.N))

    #Probability variables
    dimension = ra*rb*rb*mA*2 #Dimension of behaviour vector; last factor of 2 represents the contexts of Bob.
    prob = picos.RealVariable('p', dimension)

    #PROBABILITY CONSTRAINTS
    dict_tuple_to_index, dict_index_to_tuple = Build_Dict_Probability_Indexes_2V()

    #non-negativity constraint
    P.add_constraint(prob >= 0)
    
    #Normalization constraints
    BehaviourNormalization = []
    for A in [0, 1]:
        for B in [[0, 1], [1, 2]]:
            Constraint = 0
            for a in [0,1]:
                for b in itertools.product([0,1], repeat = 2):
                    Constraint += prob[dict_tuple_to_index[a,b[0],b[1],A,B[0],B[1]]]
            BehaviourNormalization.append(Constraint)
    P.add_list_of_constraints([Constraint == 1 for Constraint in BehaviourNormalization])

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
                    Constraint += prob[dict_tuple_to_index[a,b[0],b[1],A,B[0],B[1]]]
                for b in itertools.product([0,1], repeat = 2):
                    Constraint -= prob[dict_tuple_to_index[a,b[0],b[1],A,Bnext[0],Bnext[1]]]
                NonSignalingConstraints.append(Constraint)
    #For Bob
    for B in [[0, 1], [1, 2]]:
        for b in itertools.product([0,1], repeat = 2):
            Constraint = 0
            for a in [0,1]:
                Constraint += prob[dict_tuple_to_index[a,b[0],b[1],0,B[0],B[1]]]
            for a in [0,1]:
                Constraint -= prob[dict_tuple_to_index[a,b[0],b[1],1,B[0],B[1]]]
            NonSignalingConstraints.append(Constraint)
    P.add_list_of_constraints([Constraint == 0 for Constraint in NonSignalingConstraints])

    #Non-disturbing constraints
    NonDisturbingContraints = []
    for A in [0,1]:
        for B in [1]: #Only B1 is in different contexts
            for a in [0,1]:
                for b in [0,1]:
                    Constraint = 0

                    for b2 in [0,1]:
                        Constraint += prob[dict_tuple_to_index[a, b, b2, A, B, 2]]

                    for b0 in [0,1]:
                        Constraint -= prob[dict_tuple_to_index[a, b0, b, A, 0, B]]
                    
                    NonDisturbingContraints.append(Constraint)
    P.add_list_of_constraints([Constraint == 0 for Constraint in NonDisturbingContraints])

    
    #MOMENT MATRIX CONSTRAINTS
    for i in range(MatrixIndexes.N):
        for j in range(i, MatrixIndexes.N):
            #print(i,j)
            constraint = MatrixIndexes.VerifyConstraint(i, j)
            if constraint != None:
                if constraint == '0':
                    P.add_constraint(M[i,j] == 0)
                elif constraint == '1':
                    P.add_constraint(M[i,j] == 1)
                else:
                    results, measurements = constraint.split('|')
                    prob_type = ''.join([measurements[k] for k in range(0,len(measurements),2)])
                    prob_index = tuple(int(results[k]) for k in range(len(results))) + tuple(int(measurements[k]) for k in range(1,len(measurements),2))
                    P.add_constraint(M[i,j] == prepare_behaviours_PICOS_constraint(prob, prob_type, prob_index, dict_tuple_to_index))
    
    
    P.add_constraint(M >> 0)


    #OBJECTIVE
    obj = 0

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

        if B0 == 0 and B1 == 3:
            B0 = 3
            B1 = 0
            btemp = b0
            b0 = b1
            b1 = btemp

        obj += coefficient*prob[dict_tuple_to_index[a, b0, b1, A, B0, B1]]
    
    P.set_objective('max', obj)
    
    S = P.solve()
    print("Status: ",P.status)

    S.apply()
    print("Max inequality value: ",obj.value)
    


class MomentMatrix():
    def __init__(self, ma, mb, numA, numB, level, contexts):
        self.MaxLevel = level #Maximum NPA level considered

        self.num = {'A':numA, 'B':numB} #number of projector for Alice and Bob
        self.m = {'A':ma, 'B':mb} #number of measurements of Alice and Bob
        self.LastParty = 'B' #FUTURE UPDATE: To generalize to more parties, this should be an argument of the class
        self.contexts = contexts

        self.sequencegroups = []
        self.CalculateSequenceGroups()

        self.sizeQ = {} #size of each NPA set considered
        self.CalculateSetSizes()
        self.N = sum(self.sizeQ.values()) #Dimension of moment matrix

        self.lowerindex = {} #dictonary to save the lower index of each sequence group
        self.upperindex = {} #dictonary to save the upper index of each sequence group
        self.CalculateIndexLimits()

        self.sequences = {}
        self.CalculateSequences()
        self.SimplifySequences()

    def CalculateSetSizes(self):
        #initially indexes for each size (excluding for Q0)
        for size in range(self.MaxLevel + 1):
            self.sizeQ[size] = 0

        for group in self.sequencegroups:
            if group == '1':
                self.sizeQ[0] = 1
            else:
                quantity = 1
                for party in group:
                    quantity *= self.num[party]
                if len(group) == 4:
                    print(group, quantity)
                self.sizeQ[len(group)] += quantity


    def CalculateSequenceGroups(self):
        for length in range(self.MaxLevel + 1):
            if length == 0:
                self.sequencegroups.append('1')
                #self.lowerindex['1'] = 0
                #previousgroup = '1'
                #previousgroup_quantity = 1
            else:
                #inital group
                group = list(np.full(length, 'A'))
                self.sequencegroups.append(''.join(group))
                #self.lowerindex[''.join(group)] = self.lowerindex[previousgroup] + previousgroup_quantity
                #self.upperindex[previousgroup] = self.lowerindex[''.join(group)]

                #previousgroup = ''.join(group)
                #previousgroup_quantity = self.CalculateGroupQuantity(''.join(group))
                
                sum_index = length - 1
                while sum_index >= 0:
                    if ord(group[sum_index]) < ord(self.LastParty): #compare ascii codes
                        group[sum_index] = chr(ord(group[sum_index]) + 1)
                        for next_index in range(sum_index + 1, length):
                            group[next_index] = group[sum_index]
                        self.sequencegroups.append(''.join(group))
                        #self.lowerindex[''.join(group)] = self.lowerindex[previousgroup] + previousgroup_quantity
                        #self.upperindex[previousgroup] = self.lowerindex[''.join(group)]

                        #previousgroup = ''.join(group)
                        #previousgroup_quantity = self.CalculateGroupQuantity(''.join(group))

                        sum_index = length - 1
                    else:
                        sum_index -= 1
        #self.upperindex[previousgroup] = self.N

    def CalculateIndexLimits(self):
        for group in self.sequencegroups:
            if group == '1':
                self.lowerindex['1'] = 0
                previousgroup = '1'
                previousgroup_quantity = 1
            else:
                self.lowerindex[group] = self.lowerindex[previousgroup] + previousgroup_quantity
                self.upperindex[previousgroup] = self.lowerindex[group]

                previousgroup = group
                previousgroup_quantity = self.CalculateGroupQuantity(group)
        self.upperindex[previousgroup] = self.N



    def CalculateGroupQuantity(self, group):
        quantity = 1
        for i in range(len(group)):
            if group[i] == 'A':
                quantity *= self.num['A']
            elif group[i] == 'B':
                quantity *= self.num['B']
        return quantity


    def GetLevel(self, index):
        for level in range(self.MaxLevel + 1):
            if index < sum(sizeQ[:(level + 1)]):
                return level 

    def GetSequenceGroup(self, index):
        for group in self.lowerindex.keys():
            if index >= self.lowerindex[group] and index < self.upperindex[group]:
                return group

    def BuildContext(self, sequence):
        context = []
        for projector in sequence:
            context.append(projector[2] + str(projector[1]))
        return SortMeasurements(context)

    def CalculateSequences(self):
        index = 0
        while index < self.N: #iterate in each row of the matrix
            SequenceGroup = self.GetSequenceGroup(index)
            GroupSize = len(SequenceGroup)

            #initial sequence
            sequence = []
            for party in SequenceGroup:
                sequence.append([0, 0, party]) #first entry: result; second entry: measurement; third entry: party;

            self.sequences[index] = [tuple(projector) for projector in sequence]
            index += 1

            while index < self.upperindex[SequenceGroup]: #iterate in each sequence of the group
                elementindex = GroupSize - 1
                while elementindex >= 0: #iterate in each element of the sequence, until it is found where to sum.
                    if sequence[elementindex][0] < (2 - 1 - 1): #Check if it is not the last result considered (Colins-Gisin representation, last result is discarded)
                        #With dichotomic measurements, this condition is never satisfied (UPDATE IN THE FUTUTRE)
                        sequence[elementindex][0] += 1

                        for next_index in range(elementindex + 1, GroupSize): #reset next projectors
                            sequence[next_index] = [0, 0, sequence[next_index][2]]
                        
                        break
                    elif sequence[elementindex][1] < (self.m[sequence[elementindex][2]] - 1): #Check if it is not the last measurement
                        sequence[elementindex][1] += 1
                        sequence[elementindex][0] = 0 #first result

                        for next_index in range(elementindex + 1, GroupSize): #reset next projectors
                            sequence[next_index] = [0, 0, sequence[next_index][2]]

                        break
                    elementindex -= 1

                self.sequences[index] = [tuple(projector) for projector in sequence]
                index += 1

    def SimplifySequences(self):
        #FUTURE UPDATE: set null sequences
        #For now: search for identity and for equal projectors
        for index in range(len(self.sequences)):
            size = len(self.sequences[index])
            #search for identity
            i = 0
            while i < size:
                if self.sequences[index][i][2] == '1' and size > 1:
                    self.sequences[index].pop(i)
                    size = len(self.sequences[index])
                else:
                    i += 1

            #search for equal projectors
            i = 0
            while i < size:
                j = i + 1
                while j < size:
                    if self.sequences[index][i][0] == self.sequences[index][j][0] and self.sequences[index][i][1] == self.sequences[index][j][1] and self.sequences[index][i][2] == self.sequences[index][j][2]: #Projectors i and j are the same
                        if self.VerifyPath(self.sequences[index], i, j): #Verify if it is possible to 'walk' from one projector to the other (i.e. if they commute with the intermediate projectors)
                            #remove projector j
                            self.sequences[index].pop(j)
                            size = len(self.sequences[index])
                        else:
                            break
                    else:
                        j += 1
                i += 1

    
    def GetSequence(self, index):
        return self.sequences[index]

    def VerifyConstraint(self, index_row, index_column):
        sequence_row = list(reversed(self.GetSequence(index_row))) #Hermitian conjugate
        sequence_column = self.GetSequence(index_column)
        sequence = sequence_row + sequence_column
        size = len(sequence)
        
        if size > 1:
            #verify if identity operator is present
            i = 0
            while i < size:
                if sequence[i][2] == '1' and size > 1:
                    sequence.pop(i)
                    size = len(sequence)
                else:
                    i += 1
        
        if size == 1: #Should not be an 'else' of the above 'if'
            #UPDATE THIS PART LATER
            if sequence[0][2] == '1': #Normalization factor
                return '1'
            else: #One measurement probability
                return str(sequence[0][0]) + '|' + sequence[0][2] + str(sequence[0][1])
        
        #Verify if the sequence is null or if there are equal projectors
        i = 0
        while i < size:
            j = i + 1
            while j < size:
                if sequence[i][0] == sequence[j][0] and sequence[i][1] == sequence[j][1] and sequence[i][2] == sequence[j][2]: #Projectors i and j are the same
                    if self.VerifyPath(sequence, i, j): #Verify if it is possible to 'walk' from one projector to the other (i.e. if they commute with the intermediate projectors)
                        #remove projector j
                        sequence.pop(j)
                        size = len(sequence)
                    else:
                        break
                elif sequence[i][0] != sequence[j][0] and sequence[i][1] == sequence[j][1] and sequence[i][2] == sequence[j][2]: #Projectors are orthogonal
                    if self.VerifyPath(sequence, i, j): #Verify if it is possible to 'walk' from one projector to the other (i.e. if they commute with the intermediate projectors)
                        #the sequence is null
                        return '0'
                    else:
                        break
                else:
                    j += 1
            i += 1
        
        #Verify if the remaining sequence is a measurement context
        sorted_sequence = self.SortSequence(sequence)
        if self.contexts.IsContext(self.BuildContext(sorted_sequence)):
            results = ''
            measurements = ''
            for projector in sorted_sequence:
                results += str(projector[0])
                measurements += projector[2] + str(projector[1])
            return results + '|' + measurements
        else:
            return None


    def VerifyPath(self, sequence, i, j):
        index = i + 1
        Measurement_i = sequence[i][2] + str(sequence[i][1])
        while index < j:
            Measurement_temp = sequence[index][2] + str(sequence[index][1])
            if self.contexts.IsContext([Measurement_i, Measurement_temp]) == False:
                return False
            index += 1
        return True

    def SortSequence(self, sequence):
        sequence = sequence.copy()
        #bubble sort
        exchanged = True
        while exchanged:
            exchanged = False
            last_index = len(sequence) - 1
            for i in range(last_index):
                if ord(sequence[i][2]) > ord(sequence[i + 1][2]): #check alphabetic order using ascii code
                    sequence[i], sequence[i + 1] = sequence[i + 1], sequence[i]
                    exchanged = True
                elif ord(sequence[i][2]) == ord(sequence[i + 1][2]): #measurements of the same party
                    if sequence[i][1] > sequence[i + 1][1]: #sort measurements
                        sequence[i], sequence[i + 1] = sequence[i + 1], sequence[i]
                        exchanged = True
                    elif sequence[i][1] == sequence[i + 1][1]: #same emasurements
                        if sequence[i][0] > sequence[i + 1][0]: #sort results
                            sequence[i], sequence[i + 1] = sequence[i + 1], sequence[i]
                            exchanged = True
            last_index -= 1
        return sequence
    

class Contexts():
    def __init__(self, MaximalContexts):
        #All maximal contexts should be present, including between Alice and Bob.
        self.MaximalContexts = MaximalContexts

        #Sort maximal contexts
        for i in range(len(self.MaximalContexts)):
            self.MaximalContexts[i] = SortMeasurements(self.MaximalContexts[i])

        self.contexts = {}
        #Prepare contexts lists (one list for each context size)
        self.max_length = 0
        for context in self.MaximalContexts:
            if len(context) > self.max_length:
                self.max_length = len(context)
        for length in range(1, self.max_length + 1):
            self.contexts[length] = []

        self.BuildContexts()

    
    def BuildContexts(self):
        for maximalcontext in self.MaximalContexts:
            for context in powerset(maximalcontext):
                if len(context) == 0:
                    continue
                context = list(context)
                if context not in self.contexts[len(context)]:
                    self.contexts[len(context)].append(context)

    def IsContext(self, context):
        if len(context) > self.max_length:
            return False
        
        context = SortMeasurements(context)
        if context in self.contexts[len(context)]:
            return True
        else:
            return False
    

def prepare_behaviours(p):
    pA = {}
    for A in range(2):
        for a in range(2):
            pA[a,A] = 0.0
            for b in itertools.product([0,1], repeat = 2):
                pA[a,A] += p[a,b[0],b[1],A,0,1]
    
    pB = {}
    for B in range(2):
        for b in range(2):
            pB[b,B] = 0.0
            for a in [0,1]:
                for btemp in [0,1]:
                    pB[b,B] += p[a,b,btemp,0,B,B+1]
    for b in range(2):
        pB[b,2] = 0.0
        for a in [0,1]:
            for btemp in [0,1]:
                pB[b,2] += p[a,btemp,b,0,1,2]

    pAB = {}
    for A in range(2):
        for B in range(2):
            for a in range(2):
                for b in range(2):
                    pAB[a,b,A,B] = 0.0
                    for btemp in [0,1]:
                        pAB[a,b,A,B] += p[a,b,btemp,A,B,B+1]
    for A in range(2):
        for a in range(2):
            for b in range(2):
                pAB[a,b,A,2] = 0.0
                for btemp in [0,1]:
                    pAB[a,b,A,2] += p[a,btemp,b,A,1,2]

    pBB = {}
    for B in [[0,1], [1,2]]:
        for b in itertools.product([0,1], repeat = 2):
            pBB[b[0],b[1],B[0],B[1]] = 0.0
            for a in range(2):
                pBB[b[0],b[1],B[0],B[1]] += p[a,b[0],b[1],A,B[0],B[1]]

    return pA, pB, pAB, pBB

def prepare_behaviours_V2(p):
    behaviour = {}
    pA = {}
    for A in range(2):
        for a in range(2):
            pA[a,A] = 0.0
            for b in itertools.product([0,1], repeat = 2):
                pA[a,A] += p[a,b[0],b[1],A,0,1]
    behaviour['A'] = pA
    
    pB = {}
    for B in range(2):
        for b in range(2):
            pB[b,B] = 0.0
            for a in [0,1]:
                for btemp in [0,1]:
                    pB[b,B] += p[a,b,btemp,0,B,B+1]
    for b in range(2):
        pB[b,2] = 0.0
        for a in [0,1]:
            for btemp in [0,1]:
                pB[b,2] += p[a,btemp,b,0,1,2]
    behaviour['B'] = pB

    pAB = {}
    for A in range(2):
        for B in range(2):
            for a in range(2):
                for b in range(2):
                    pAB[a,b,A,B] = 0.0
                    for btemp in [0,1]:
                        pAB[a,b,A,B] += p[a,b,btemp,A,B,B+1]
    for A in range(2):
        for a in range(2):
            for b in range(2):
                pAB[a,b,A,2] = 0.0
                for btemp in [0,1]:
                    pAB[a,b,A,2] += p[a,btemp,b,A,1,2]
    behaviour['AB'] = pAB

    pBB = {}
    for B in [[0,1], [1,2]]:
        for b in itertools.product([0,1], repeat = 2):
            pBB[b[0],b[1],B[0],B[1]] = 0.0
            for a in range(2):
                pBB[b[0],b[1],B[0],B[1]] += p[a,b[0],b[1],A,B[0],B[1]]
    behaviour['BB'] = pBB

    behaviour['ABB'] = p

    return behaviour

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



def Test_Vertices_MostGeneralND(level):
    #Local, non-disturbing and non-deterministic vertice
    p = {(0, 0, 0, 0, 0, 1): 0.0, (0, 0, 1, 0, 0, 1): 0.5, (0, 1, 0, 0, 0, 1): 0.5, (0, 1, 1, 0, 0, 1): 0.0, (1, 0, 0, 0, 0, 1): 0.0, (1, 0, 1, 0, 0, 1): 0.0, (1, 1, 0, 0, 0, 1): 0.0, (1, 1, 1, 0, 0, 1): 0.0, (0, 0, 0, 0, 1, 2): 0.5, (0, 0, 1, 0, 1, 2): 0.0, (0, 1, 0, 0, 1, 2): 0.5, (0, 1, 1, 0, 1, 2): 0.0, (1, 0, 0, 0, 1, 2): 0.0, (1, 0, 1, 0, 1, 2): 0.0, (1, 1, 0, 0, 1, 2): 0.0, (1, 1, 1, 0, 1, 2): 0.0, (0, 0, 0, 1, 0, 1): 0.0, (0, 0, 1, 1, 0, 1): 0.0, (0, 1, 0, 1, 0, 1): 0.5, (0, 1, 1, 1, 0, 1): 0.0, (1, 0, 0, 1, 0, 1): 0.0, (1, 0, 1, 1, 0, 1): 0.5, (1, 1, 0, 1, 0, 1): 0.0, (1, 1, 1, 1, 0, 1): 0.0, (0, 0, 0, 1, 1, 2): 0.0, (0, 0, 1, 1, 1, 2): 0.0, (0, 1, 0, 1, 1, 2): 0.5, (0, 1, 1, 1, 1, 2): 0.0, (1, 0, 0, 1, 1, 2): 0.5, (1, 0, 1, 1, 1, 2): 0.0, (1, 1, 0, 1, 1, 2): 0.0, (1, 1, 1, 1, 1, 2): 0.0}
    NPA_2V(prepare_behaviours_V2(p), level)

    #Local and deterministic vertice
    vertice1 = '0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0'
    vertice1 = np.array(vertice1.split(" ")).astype(float)
    p1 = Vertice_To_Behaviour(vertice1)
    NPA_2V(prepare_behaviours_V2(p1), level)

    #Nonlocal NSND vertice
    vertice2 = '0.5 0.0 0.0 0.0 0.0 0.0 0.0 0.5 0.0 0.0 0.5 0.0 0.5 0.0 0.0 0.0 0.5 0.0 0.0 0.0 0.0 0.0 0.0 0.5 0.5 0.0 0.0 0.0 0.0 0.0 0.5 0.0'
    vertice2 = np.array(vertice2.split(" ")).astype(float)
    p2 = Vertice_To_Behaviour(vertice2)
    NPA_2V(prepare_behaviours_V2(p2), level)

def Test_Moment_Matrix():
    mA = 2 #Number of measurements for Alice
    ra = 2 #Number of results for Alice's measurements
    mB = 3 #Number of measurements for Bob 
    rb = 2 #Number of results of Bob's measurements

    numA = (ra - 1)*mA #Number of projectors for Alice
    numB = (rb - 1)*mB #Number of projectors for Bob

    sizeQ = []
    sizeQ.append(1) #size of set Q0
    sizeQ.append((ra - 1)*mA + (rb - 1)*mB) #size of set Q1
    sizeQ.append(((ra - 1)*mA)**2 + ((rb - 1)*mB)**2 + (ra - 1)*mA*(rb - 1)*mB) #size of set Q2

    level = 4

    ContextObj = Contexts([['A0','B0','B1'], ['A0','B1','B2'], ['A1','B0','B1'], ['A1','B1','B2']])

    MatrixIndexes = MomentMatrix(mA, mB, numA, numB, level, ContextObj)
    print(MatrixIndexes.N)
    print(MatrixIndexes.sizeQ)
    print(MatrixIndexes.sequencegroups)
    print(MatrixIndexes.lowerindex)
    print(MatrixIndexes.upperindex)

    #print(MatrixIndexes.sequences)




inequality = '-p110|A0B0B1 +p101|A0B1B2 +p010|A1B0B1 +p010|A1B1B2 +p011|A1B1B2 +p100|A1B1B2 +p110|A1B1B2 +p111|A1B1B2 <= 1'
NPA_2V_Inequality(2, inequality)

#Local, non-disturbing and non-deterministic vertice
p = {(0, 0, 0, 0, 0, 1): 0.0, (0, 0, 1, 0, 0, 1): 0.5, (0, 1, 0, 0, 0, 1): 0.5, (0, 1, 1, 0, 0, 1): 0.0, (1, 0, 0, 0, 0, 1): 0.0, (1, 0, 1, 0, 0, 1): 0.0, (1, 1, 0, 0, 0, 1): 0.0, (1, 1, 1, 0, 0, 1): 0.0, (0, 0, 0, 0, 1, 2): 0.5, (0, 0, 1, 0, 1, 2): 0.0, (0, 1, 0, 0, 1, 2): 0.5, (0, 1, 1, 0, 1, 2): 0.0, (1, 0, 0, 0, 1, 2): 0.0, (1, 0, 1, 0, 1, 2): 0.0, (1, 1, 0, 0, 1, 2): 0.0, (1, 1, 1, 0, 1, 2): 0.0, (0, 0, 0, 1, 0, 1): 0.0, (0, 0, 1, 1, 0, 1): 0.0, (0, 1, 0, 1, 0, 1): 0.5, (0, 1, 1, 1, 0, 1): 0.0, (1, 0, 0, 1, 0, 1): 0.0, (1, 0, 1, 1, 0, 1): 0.5, (1, 1, 0, 1, 0, 1): 0.0, (1, 1, 1, 1, 0, 1): 0.0, (0, 0, 0, 1, 1, 2): 0.0, (0, 0, 1, 1, 1, 2): 0.0, (0, 1, 0, 1, 1, 2): 0.5, (0, 1, 1, 1, 1, 2): 0.0, (1, 0, 0, 1, 1, 2): 0.5, (1, 0, 1, 1, 1, 2): 0.0, (1, 1, 0, 1, 1, 2): 0.0, (1, 1, 1, 1, 1, 2): 0.0}
NPA_2V(prepare_behaviours_V2(p), 2)

#Result smaller than zero means that the behaviour does not belong to the respective set of the NPA hierarchy