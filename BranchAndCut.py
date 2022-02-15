#!/usr/bin/env python3.7

import gurobipy as gp
from gurobipy import GRB

import math
import numpy

import networkx as nx
from networkx.algorithms.approximation import clique

from BinPackingData import *
from BinPacking import *
from Model import *
from OrthogonalPacking import *
from PlacementPoints import *

""" Datasets at https://github.com/Oscar-Oliveira/OR-Datasets/tree/master/Cutting-and-Packing/2D/Datasets """
import time

class BinPackingCallback:
    def __init__(self):
        pass

    @staticmethod
    def AddCut(items, itemVariables, bins, model):
        for b in bins:
            expr = gp.LinExpr()
            for i in items:
                expr += itemVariables[b][i]

            model.cbLazy(expr <= len(items) - 1)

    @staticmethod
    def SolveKnapsack2D(liftedSet, fixedItemId, fixItemToBin, bin, objectiveCoefficients, model, timeLimit = 1.0):
        itemsToFix = [fixedItemId]
        items = [model._Items[fixedItemId]]
        for i in liftedSet:
            items.append(model._Items[i])
            if fixItemToBin[i]:
                itemsToFix.append(i)

        t1 = time.time()

        knapsack2D = Knapsack2D()
        status, upperBound, objective = knapsack2D.Solve(items, itemsToFix, objectiveCoefficients, bin, timeLimit)

        t2 = time.time()

        model._KnapsackTimeLifting += (t2 - t1)

        if status == 'INFEASIBLE' or status == 'INVALID':
            raise ValueError('Knapsack2D infeasible or invalid.')
        elif status == 'UNKOWN':
            # time limit reached
            model._LiftingTimeLimit += 1
            return len(liftedSet) - 1
        elif status == 'OPTIMAL':
            model._LiftingOptimal += 1
        elif status == 'FEASIBLE':
            model._LiftingFeasible += 1
            #print(f"Lifting feasible: {[item.Id for i, item in enumerate(items)]}")

        return upperBound

    @staticmethod
    def AddLiftedCoverInequality(infeasibleItemSubset, additionalItems, liftingCoefficients, itemVariables, bins, model):
        for b in bins:
            expr = gp.LinExpr()
            for i in infeasibleItemSubset:
                expr += itemVariables[b][i]

            for i in additionalItems:
                expr += liftingCoefficients[i] * itemVariables[b][i]

            model.cbLazy(expr <= len(infeasibleItemSubset) - 1)
        
        #print(f'Add cut in bins {max(bins)}: {infeasibleItemSubset} + {[liftingCoefficients[i] for i in additionalItems]} * {additionalItems} <= {len(infeasibleItemSubset) - 1}')

    """ See section 7.3 in Cote, Iori (2021). """
    @staticmethod
    def AddLiftedCut(infeasibleItemSubset, itemVariables, binId, model):
        objectiveCoefficients = [1 if i in infeasibleItemSubset else 0 for i in range(len(itemVariables))]
        liftingCoefficients = list(objectiveCoefficients)

        liftedSet = list(infeasibleItemSubset)
        additionalItems = []
        for i in range(len(itemVariables)):
            if i in infeasibleItemSubset or i < binId or model._Preprocess.FixItemToBin[i]:
                continue

            skipItem = False
            for j in infeasibleItemSubset:
                #isIncompatible = frozenset(i, j) in model._Preprocess.IncompatibleItems
                if model._Preprocess.FixItemToBin[j] and frozenset((i, j)) in model._Preprocess.IncompatibleItems:
                    skipItem = True
                    break
            
            if skipItem:
                continue

            # check against _FeasibleSets or _InfeasibleSets

            upperBound = BinPackingCallback.SolveKnapsack2D(liftedSet, i, model._Preprocess.FixItemToBin, model._Bins[binId], objectiveCoefficients, model, 1.0)
            
            liftingCoefficient = max(0, len(infeasibleItemSubset) - 1 - upperBound)
            if liftingCoefficient > 0:
                # add to model._InfeasibleSets.add(frozenset(liftedSet)) ?
                liftedSet.append(i)
                additionalItems.append(i)

                model._InfeasibleSets.add(frozenset(liftedSet))

                liftingCoefficients[i] = liftingCoefficient
                objectiveCoefficients[i] = liftingCoefficient
            elif liftingCoefficient < 0.0:
                print(f'Negative lifting coefficient = {liftingCoefficient}')
                break
            else:
                infeasibleSet = list(liftedSet)
                infeasibleSet.append(i)
                model._InfeasibleSets.add(frozenset(infeasibleSet))

        # Not true if any coefficient < 1.0. Which should not occur when the 2D Knapsack problem is solved with CP.
        #model._InfeasibleSets.add(frozenset(liftedSet))
        minItemId = min(liftedSet)

        model._LiftedCuts += 1 if len(liftedSet) > len(infeasibleItemSubset) else 0

        upperBoundBin = model.cbGet(GRB.Callback.MIPSOL_OBJ)

        #filteredBins = [b for b, bin in enumerate(model._Bins) if b <= minItemId and b < upperBoundBin]
        #BinPackingCallback.AddLiftedCoverInequality(infeasibleItemSubset, additionalItems, liftingCoefficients, itemVariables, filteredBins, model)
        
        if model._Preprocess.FixItemToBin[binId]:
            BinPackingCallback.AddLiftedCoverInequality(infeasibleItemSubset, additionalItems, liftingCoefficients, itemVariables, [binId], model)
        else:
            filteredBins = [b for b, bin in enumerate(model._Bins) if b <= minItemId and b < upperBoundBin]# and not model._Preprocess.FixItemToBin[b]]
            BinPackingCallback.AddLiftedCoverInequality(infeasibleItemSubset, additionalItems, liftingCoefficients, itemVariables, filteredBins, model)

    @staticmethod
    def AddStrengthenedCut(items, binId, model):
        sortedItems = sorted(items, reverse = True)
        itemIndices = [item.Id for item in sortedItems]

        numberOfItems = len(items)

        isFeasible = False
        while not isFeasible:
            removedItem = sortedItems.pop()
            removedIndex = itemIndices.pop()

            if len(itemIndices) <= 1:
                # if 1, then incompatibility cuts should have been added in preprocessing
                # if 0, item should have been removed in preprocessing as large item
                raise ValueError("Error in cut strengthening.")

            if frozenset(itemIndices) in model._FeasibleSets:
                # extend feasible subset by previously removed item to make it infeasible again 
                itemIndices.append(removedIndex)
                
                model._InfeasibleSets.add(frozenset(itemIndices))

                model._StrengthenedCuts += 1 if len(itemIndices) < numberOfItems else 0

                if model._EnableLifting:
                    BinPackingCallback.AddLiftedCut(itemIndices, model._VarsX, binId, model)
                else:
                    minItemId = min(itemIndices)
                    BinPackingCallback.AddCut(itemIndices, model._VarsX, [b for b in range(len(model._Bins)) if b <= minItemId], model)

                return

            if frozenset(itemIndices) in model._InfeasibleSets or frozenset(itemIndices) in model._Preprocess.IncompatibleItems:
                #raise ValueError('No previously discovered infeasible subsets should be found during cut strengthening.')
                model._InfeasibleDoubleCount += 1
                #print(f"Double count in bin {binId} (fixed = {model._Preprocess.FixItemToBin[binId]}): {itemIndices} vs. {[item.Id for i, item in enumerate(items)]}")
                continue
                #BinPackingCallback.AddCut(itemIndices, model._VarsX, model._Bins, model)
                #return

            t1 = time.time()

            orthogonalPacking = OrthogonalPackingSolver(sortedItems, model._Bins[0], model._PlacementPointStrategy)
            isFeasible = orthogonalPacking.Solve(model._InstanceId)

            t2 = time.time()

            model._BinPackingTimeStrengthening += (t2 - t1)

            if isFeasible:
                model._FeasibleSets.add(frozenset(itemIndices))

                itemIndices.append(removedIndex)
                
                model._InfeasibleSets.add(frozenset(itemIndices))

                model._StrengthenedCuts += 1 if len(itemIndices) < numberOfItems else 0

                if model._EnableLifting:
                    BinPackingCallback.AddLiftedCut(itemIndices, model._VarsX, binId, model)
                else:
                    minItemId = min(itemIndices)
                    BinPackingCallback.AddCut(itemIndices, model._VarsX, [b for b in range(len(model._Bins)) if b <= minItemId], model)

                return
                
            model._InfeasibleSets.add(frozenset(itemIndices))

    @staticmethod
    def FindIntegerAssignments(model):
        itemIndicesArray = []
        itemsInBinArray = []
        for b, binVariables in enumerate(model._Bins):
            area = 0.0
            itemsInBin = []
            itemIndices = []
            for i, x in enumerate(model._VarsX[b]):
                xVal = model.cbGetSolution(x)
                if xVal < 0.5:
                    continue

                item = model._Items[i]

                area += item.Weight
                itemsInBin.append(item)
                itemIndices.append(i)

            itemIndices.sort()

            if area > model._Bins[b].WeightLimit:
                raise ValueError('Capacity constraints violated')

            itemIndicesArray.append(itemIndices)
            itemsInBinArray.append(itemsInBin)

        return itemIndicesArray, itemsInBinArray
            
    @staticmethod
    def AddCuts(model, itemIndicesArray, itemsInBinArray):
        infeasibleSubproblems = 0
        for b, itemIndices in enumerate(itemIndicesArray):

            if infeasibleSubproblems >= model._InfeasibleSuproblemCutThreshold:
                return

            itemsInBin = itemsInBinArray[b]

            if len(itemIndices) == 0:
                continue

            if frozenset(itemIndices) in model._FeasibleSets:
                continue

            if frozenset(itemIndices) in model._InfeasibleSets:# or frozenset(itemIndices) in model._Preprocess.IncompatibleItems:
                #raise ValueError('Infeasible items occured twice')
                #print(f"Double count in bin {b} (fixed = {model._Preprocess.FixItemToBin[b]}): {itemIndices}")
                model._InfeasibleDoubleCount += 1
                infeasibleSubproblems += 1 # only count new cuts
                BinPackingCallback.AddCut(itemIndices, model._VarsX, [b], model)
                continue

            t1 = time.time()

            orthogonalPacking = OrthogonalPackingSolver(itemsInBin, model._Bins[b], model._PlacementPointStrategy)
            isFeasible = orthogonalPacking.Solve(model._InstanceId)

            t2 = time.time()

            model._BinPackingTimeSeparation += (t2 - t1)

            if isFeasible:
                model._FeasibleSets.add(frozenset(itemIndices))
            else:
                model._InfeasibleSets.add(frozenset(itemIndices))
                model._CutCount += 1
                
                # Only add cut if sensible w.r.t. compatible items
                # Only add cuts for bins < UB
                if model._EnableCutStrengthening:
                    BinPackingCallback.AddStrengthenedCut(itemsInBin, b, model)
                else:
                    minItemId = min(itemIndices)
                    BinPackingCallback.AddCut(itemIndices, model._VarsX, [id for id, bin in enumerate(model._Bins) if id <= minItemId], model)

                infeasibleSubproblems += 1

    @staticmethod
    def callback(model, where):
        if where == GRB.Callback.MIPSOL:
            itemIndicesArray, itemsInBinArray = BinPackingCallback.FindIntegerAssignments(model)
            BinPackingCallback.AddCuts(model, itemIndicesArray, itemsInBinArray)

class BinPackingMip:
    def __init__(self, enable2D = True):
        self.Model = gp.Model("BinPacking")

        self.Model.Params.OutputFlag = 0
        if enable2D:
            self.Model.Params.OutputFlag = 0
            
        self.Model.Params.lazyConstraints = 1
        self.Model.Params.MIPFocus = 2
        self.Model.Params.TimeLimit = 1800 if enable2D else 180

        self.Items = []
        self.Bins = []

        self.ItemVariables = []
        self.BinVariables = []

        self.Callback = BinPackingCallback()

        self.Enable2D = enable2D

        self.LowerBoundBin = 1

        self.Preprocess = None
        
        InstanceId = -1

    def AddItems(self, items):
        self.Items = items

    def AddBins(self, binDxArray, binDyArray):
        for i in range(len(binDxArray)):
            dx = binDxArray[i]
            dy = binDyArray[i]

            item = Bin(dx, dy)
            self.Bins.append(item)

    def CreateVariables(self):
        for b, bin in enumerate(self.Bins):
            self.ItemVariables.append([])
            for i, item in enumerate(self.Items):
                itemVar = self.Model.addVar(vtype=GRB.BINARY, name=f'x{b}{i}')

                self.ItemVariables[b].append(itemVar)
                #print(item.Weight)
                
        for b, bin in enumerate(self.Bins):
            binVar = self.Model.addVar(obj=1.0, vtype=GRB.BINARY, name=f'y{b}')

            self.BinVariables.append(binVar)
            #print(bin.WeightLimit)

        w = [item.Dx for item in self.Items]
        h = [item.Dy for item in self.Items]

        totalArea = numpy.dot(w,h)
        lowerBoundBinArea = math.ceil(float(totalArea) / float(self.Bins[0].Dx * self.Bins[0].Dy)) # homogeneous bins!

        lb1D = lowerBoundBinArea
        if self.Enable2D:
            binPacking1D = BinPackingMip(False)
            binPacking1D.Preprocess = self.Preprocess

            binPacking1D.AddItems(list(self.Items))
            binPacking1D.AddBins([self.Bins[0].Dx] * len(self.Bins), [self.Bins[0].Dy] * len(self.Bins))

            binPacking1D.CreateVariables()
            binPacking1D.CreateConstraints()

            binPacking1D.Solve()

            model1D = binPacking1D.Model  
            lb1D = math.ceil(model1D.objBound)
            # TODO add bound lifted indicator to solution statistic

        self.LowerBoundBin = max(lowerBoundBinArea, lb1D)

        for b in range(self.LowerBoundBin):
            self.BinVariables[b].LB = 1.0

    def CreateConstraints(self):
        for i, item in enumerate(self.Items):
            binExpr = gp.LinExpr()
            for b, bin in enumerate(self.Bins):
                binExpr += self.ItemVariables[b][i]
            self.Model.addConstr(binExpr == 1)

        for b, bin in enumerate(self.Bins):
            """
            for x in self.ItemVariables[b]:
                itemsInBin += x

            #itemsInBin.AddTerm(self.ItemVariables[b])
            self.Model.addConstr(itemsInBin == 1)
            """

            itemsInBin = gp.LinExpr()
            for i, item in enumerate(self.Items):
                #self.Model.addConstr(x * item.Weight <= bin.WeightLimit * self.BinVariables[b])

                itemsInBin += self.ItemVariables[b][i] * item.Weight

            self.Model.addConstr(itemsInBin <= bin.WeightLimit * self.BinVariables[b])

            # symmetry breaking, only valid if bins are homogeneous
            if b < len(self.Bins) - 1:
                self.Model.addConstr(self.BinVariables[b + 1] <= self.BinVariables[b])

            # symmetry breaking
            for i, item in enumerate(self.Items):
                if b > i:
                    self.Model.addConstr(self.ItemVariables[b][i] == 0)

    def SetCallbackData(self):
        self.Model._Items = self.Items
        self.Model._Bins = self.Bins
        self.Model._VarsX = self.ItemVariables
        self.Model._Preprocess = self.Preprocess

        self.Model._FeasibleSets = set()
        self.Model._InfeasibleSets = set()

        self.Model._CutCount = 0
        self.Model._StrengthenedCuts = 0
        self.Model._LiftedCuts = 0
        self.Model._InfeasibleDoubleCount = 0
        self.Model._LiftingOptimal = 0
        self.Model._LiftingFeasible = 0
        self.Model._LiftingTimeLimit = 0
        self.Model._KnapsackTimeLifting = 0.0
        self.Model._BinPackingTimeStrengthening = 0.0
        self.Model._BinPackingTimeSeparation = 0.0
        self.Model._EnableLifting = False
        self.Model._EnableCutStrengthening = False
        self.Model._InfeasibleSuproblemCutThreshold = 1
        
        self.Model._EnablePreprocessLifting = False
        self.Model._PlacementPointStrategy = PlacementPointStrategy.NormalPatterns

        self.Model._InstanceId = self.InstanceId

    def DeterminePositions(self, itemIndicesArray, itemsInBinArray):
        rectanglesArray = []
        for b, itemIndices in enumerate(itemIndicesArray):
            numberOfItems = len(itemIndices)
            if numberOfItems == 0:
                continue

            itemsInBin = itemsInBinArray[b]

            bin = self.Model._Bins[b]

            orthogonalPacking = OrthogonalPackingSolver(itemsInBin, bin, self.Model._PlacementPointStrategy) #
            isFeasible = orthogonalPacking.Solve(self.InstanceId)

            if not isFeasible:
                raise ValueError('Packing must not be infeasible during solution extraction.')
            
            xArray = [(b * bin.Dx) + orthogonalPacking.PositionsX[i] for i in range(numberOfItems)]
            yArray = [orthogonalPacking.PositionsY[i] for i in range(numberOfItems)]

            rectangles = ExtractDataForPlot(xArray, yArray, itemsInBin, bin.Dx, bin.Dy)

            rectanglesArray.extend(rectangles)

        return rectanglesArray

    @staticmethod
    def FindIntegerAssignments(model):
        itemIndicesArray = []
        itemsInBinArray = []
        for b, binVariables in enumerate(model._Bins):
            area = 0.0
            itemsInBin = []
            itemIndices = []
            for i, x in enumerate(model._VarsX[b]):
                xVal = x.x
                if xVal < 0.5:
                    continue

                item = model._Items[i]

                area += item.Weight
                itemsInBin.append(item)
                itemIndices.append(i)

            itemIndices.sort()

            if area > model._Bins[b].WeightLimit:
                raise ValueError('Capacity constraints violated')

            itemIndicesArray.append(itemIndices)
            itemsInBinArray.append(itemsInBin)

        return itemIndicesArray, itemsInBinArray

    def ExtractSolution(self):
        itemIndicesArray, itemsInBinArray = BinPackingMip.FindIntegerAssignments(self.Model)

        rectangles = self.DeterminePositions(itemIndicesArray, itemsInBinArray)

        return rectangles

    def SetStartSolution(self, itemBinAssignments, lb, ub):
        for b, bin in enumerate(self.Bins):
            if b < math.ceil(lb):
                self.BinVariables[b].LB = 1.0

            if b < int(ub):
                self.BinVariables[b].Start = 1.0

            if b >= int(ub):
                self.BinVariables[b].UB = 0.0
                for i, item in enumerate(self.Items):
                    self.ItemVariables[b][i].UB = 0.0

            for i, item in enumerate(self.Items):
                self.ItemVariables[b][i].Start = 0.0

        assignments = {}
        for b in range(int(ub)):
            assignments[b] = []

        for i, b in enumerate(itemBinAssignments):
            #self.ItemVariables[b][i].Start = 1.0
            assignments[b].append(i)

        sortedAssigments = sorted(assignments.values())

        for b, assigment in enumerate(sortedAssigments):
            for i in assigment:
                self.ItemVariables[b][i].Start = 1.0

    def FixIncompatibleItems(self, fixItemToBin):
        numberOfItems = len(self.Items)
        # This is a similar logic as in section 6.2 in Cote, Haouari, Iori (2019). 
        for i in range(0, numberOfItems):
            if fixItemToBin[i]:
                for b, bin in enumerate(self.Bins):
                    if b == i:
                        self.ItemVariables[b][i].LB = 1.0
                    else:
                        self.ItemVariables[b][i].UB = 0.0
            else:
                return

    def AddLiftedCoverInequality(self, infeasibleItemSubset, additionalItems, liftingCoefficients, itemVariables, bins, model):
        for b in bins:
            expr = gp.LinExpr()
            for i in infeasibleItemSubset:
                expr += itemVariables[b][i]

            for i in additionalItems:
                expr += liftingCoefficients[i] * itemVariables[b][i]

            model.addConstr(expr <= len(infeasibleItemSubset) - 1)

    def AddLiftedCut(self, infeasibleItemSubset, itemVariables, model):

        minIdInitial = min(infeasibleItemSubset)
        for b, bin in enumerate(model._Bins):
            if b > minIdInitial:
                continue
            
            objectiveCoefficients = [1 if i in infeasibleItemSubset else 0 for i in range(len(itemVariables))]
            liftingCoefficients = list(objectiveCoefficients)

            liftedSet = list(infeasibleItemSubset)
            additionalItems = []
            for i in range(len(itemVariables)):
                if i in infeasibleItemSubset or i < b or (model._Preprocess.FixItemToBin[i] and b != i):# or i < minIdInitial: or model._Preprocess.FixItemToBin[i] 
                    continue

                skipItem = False
                for j in liftedSet:
                    #isIncompatible = frozenset(i, j) in model._Preprocess.IncompatibleItems
                    if (model._Preprocess.FixItemToBin[j] and frozenset((i, j)) in model._Preprocess.IncompatibleItems):
                        skipItem = True
                        break
                
                if skipItem:
                    continue

                upperBound = BinPackingCallback.SolveKnapsack2D(liftedSet, i, model._Preprocess.FixItemToBin, model._Bins[0], objectiveCoefficients, model, 1.0)
                
                liftingCoefficient = max(0, len(infeasibleItemSubset) - 1 - upperBound)
                if liftingCoefficient > 0:
                    liftingCoefficients[i] = liftingCoefficient
                    liftedSet.append(i)
                    additionalItems.append(i)

                    model._InfeasibleSets.add(frozenset(liftedSet))

                    objectiveCoefficients[i] = liftingCoefficient

            #model._InfeasibleSets.add(frozenset(liftedSet))
            minItemId = min(liftedSet)

            model._LiftedCuts += 1 if len(liftedSet) > len(infeasibleItemSubset) else 0

            upperBoundBin = len(model._Items) #model.cbGet(GRB.Callback.MIPSOL_OBJ)

            #filteredBins = [b for b, bin in enumerate(model._Bins) if b <= minIdInitial and b < upperBoundBin]
            filteredBins = [b]
            self.AddLiftedCoverInequality(infeasibleItemSubset, additionalItems, liftingCoefficients, itemVariables, filteredBins, model)

        return liftedSet, liftingCoefficients

    def AddIncompatibilityCuts(self, incompatibleItems, fixItemToBin):
        #liftedSets = set()
        for incompatibleSet in incompatibleItems:
            allFixed = True
            for k in incompatibleSet:
                if not fixItemToBin[k]:
                    allFixed = False
                    break

            #if fixItemToBin[i] and fixItemToBin[j]:
            #    continue
            if allFixed:
                continue
            
            minItemId = min(incompatibleSet)

            liftedSet, liftingCoefficients = self.AddLiftedCut(incompatibleSet, self.ItemVariables, self.Model)
            #liftedSets.add(frozenset(liftedSet))
            
            """
            for b, bin in enumerate(self.Bins):
                if b > minItemId:
                    break
                
                # are the same cuts added multiple times if contained multiple times in incompatibleItems?
                # TODO.Logic: lifting 
                #self.Model.addConstr(self.ItemVariables[b][i] + self.ItemVariables[b][j] <= 1)
                expr = gp.LinExpr()
                for k in incompatibleSet:
                    expr += liftingCoefficients[k] * self.ItemVariables[b][k]

                self.Model.addConstr(expr <= len(incompatibleSet) - 1)
            """

    def ApplyPreprocessChanges(self, preprocess):
        fixItemToBin = preprocess.FixItemToBin
        self.FixIncompatibleItems(fixItemToBin)

        self.SetCallbackData()

        incompatibleItems = preprocess.IncompatibleItems
        self.AddIncompatibilityCuts(incompatibleItems, fixItemToBin)
        
    def Solve(self):
        #self.SetCallbackData()

        if self.Enable2D:
            self.Model.optimize(BinPackingCallback.callback)
        else:
            self.Model.optimize()

        #self.Model.printAttr('x')

        #self.Model.write('BPP.lp')

        statusCode = 1 if self.Model.Status == GRB.OPTIMAL else 0
        
        rectangles = None
        if self.Enable2D:
            rectangles = self.ExtractSolution()

            """
            print(f'Strengthened cuts: {self.Model._StrengthenedCuts}/{self.Model._CutCount}')
            print(f'Lifted cuts: {self.Model._LiftedCuts}/{self.Model._CutCount}')
            print(f'Lifting optimal: {self.Model._LiftingOptimal}')
            print(f'Lifting feasible: {self.Model._LiftingFeasible}')
            print(f'Lifting optimal: {self.Model._LiftingTimeLimit}')
            print(f'Bin packing time separation: {self.Model._BinPackingTimeSeparation}')
            print(f'Bin packing time strengthening: {self.Model._BinPackingTimeStrengthening}')
            print(f'Knapsack time lifting: {self.Model._KnapsackTimeLifting}')
            """

        return rectangles

class BinPackingBranchAndCutSolver:
    def __init__(self, instanceId, placementPointStrategy = PlacementPointStrategy.UnitDiscretization):
        self.BinPacking = BinPackingMip()
        self.BinPacking.InstanceId = instanceId

        self.IsOptimal = False
        self.LB = -1
        self.UB = -1
        self.Runtime = -1
        self.SolverType = ""
        self.CutCount = -1

        self.PlacementPointStrategy = placementPointStrategy
        self.preprocess = None

    def RetrieveSolutionStatistics(self):
        if self.IsOptimal:
            return 

        model = self.BinPacking.Model

        self.IsOptimal = 1 if model.Status == GRB.OPTIMAL else 0
        self.LB = model.objBound
        self.UB = model.objVal
        self.Runtime = model.Runtime
        self.SolverType = "B&C"

    def DetermineStartSolution(self, items, H, W, lowerBoundBin):
        solverCP = BinPackingSolverCP(items, H, W, lowerBoundBin, len(items), self.PlacementPointStrategy, 180, False, self.preprocess)
        rectangles = solverCP.Solve('OneBigBin')

        if solverCP.LB == solverCP.UB:
            return True, solverCP.LB, "CP", rectangles

        self.BinPacking.SetStartSolution(solverCP.ItemBinAssignments, solverCP.LB, solverCP.UB)

        return False, solverCP.LB, "B&C", rectangles

    def Run(self, items, H, W):   
        bin = Bin(W, H)

        self.preprocess = Preprocess(items, bin, self.PlacementPointStrategy)
        self.preprocess.Run()
        newItems = self.preprocess.ProcessedItems
        numberOfItems = len(newItems)

        self.BinPacking.Preprocess = self.preprocess
        self.BinPacking.AddItems(newItems)
        self.BinPacking.AddBins([W] * numberOfItems, [H] * numberOfItems)

        self.BinPacking.CreateVariables()
        self.BinPacking.CreateConstraints()

        self.BinPacking.ApplyPreprocessChanges(self.preprocess)
        
        #self.BinPacking.Model.write('BPP.lp')

        isOptimalCP, lbCP, solverType, rectangles = self.DetermineStartSolution(newItems, H, W, self.BinPacking.LowerBoundBin)

        if isOptimalCP:
            self.IsOptimal = 1
            self.LB = lbCP + len(self.preprocess.RemovedItems)
            self.UB = lbCP + len(self.preprocess.RemovedItems)
            self.SolverType = solverType
            self.Runtime = -1

            return rectangles

        rectangles = self.BinPacking.Solve()

        return rectangles

def main():
    solutions = {}
    # Single bin 2D-BPP CP model takes ages to prove feasibility/infeasibility on instances: 173, 262, 292, 297, 298, 317
    # Meet in the middle produces erroneous 1D KP instances for instances 11
    # Postsolve error: 51, 45, 125/126, 31, 26 or 27, 311
    # Double count: 21, 22, 26, 27, 31, 32
    # Negative lifting coefficient: 120, 123, 
    #hardInstances = [226, 232, 242, 242, 244, 245, 247, 248, 249, 261, 292, 313, 314, 332, 173, 187, 188, 191, 197, 142, 143, 145, 146, 149]
    #mediumInstance = [149, 174]

    for instance in range(87, 500):
    #for instance in hardInstances:
        currentInstanceId = instance
        items, H, W = ReadBenchmarkData(instance)
        
        solver = BinPackingBranchAndCutSolver(instance)
        rectangles = solver.Run(items, H, W)

        solver.RetrieveSolutionStatistics()

        # TODO: introduce solution statistics struct
        bestBoundMIP = solver.LB
        upperBoundMIP = solver.UB
        solverType = solver.SolverType
        isOptimalMIP = solver.IsOptimal
        
        PlotSolution(upperBoundMIP * W, H, rectangles)

        if isOptimalMIP:
            print(f'Instance {instance}: Optimal solution = {int(bestBoundMIP)} found by {solverType} (#items = {len(items)})')
        else:
            print(f'Instance {instance}: No optimal solution found, [lb, ub] = [{bestBoundMIP}, {upperBoundMIP}] (#items = {len(items)})')
    
        solutions[instance] = {'LB': bestBoundMIP, 'UB': upperBoundMIP, 'Solver': solverType}

    solutionsJson = json.dumps(solutions, indent = 4)
    with open("Solutions.json", "w") as outfile:
        outfile.write(solutionsJson)

if __name__ == "__main__":
    main()