#!/usr/bin/env python3.7

import gurobipy as gp
from gurobipy import GRB

import math
import numpy
import json

from ortools.sat.python import cp_model
from ortools.sat.python.cp_model import Domain

import networkx as nx
from networkx.algorithms.approximation import clique

from BinPackingData import *
from ErwinCP import *

""" Datasets at https://github.com/Oscar-Oliveira/OR-Datasets/tree/master/Cutting-and-Packing/2D/Datasets """
import time

from enum import IntEnum

class PlacementPointStrategy(IntEnum):
    UnitDiscretization = 0
    NormalPatterns = 1
    MeetInTheMiddlePatterns = 2
    MinimalMeetInTheMiddlePatterns = 3

class Item:
    def __init__(self, id, dx, dy):
        self.Id = id
        self.Dx = dx
        self.Dy = dy
        self.Weight = dx * dy
        
    def __eq__(self, other):
        return self.Dx == other.Dx and self.Dy == other.Dy

    def __lt__(self, other):
        return ((self.Weight, self.Dx) < (other.Weight, other.Dx))
        

class Bin:
    def __init__(self, dx, dy):
        self.Dx = dx
        self.Dy = dy
        self.WeightLimit = dx * dy

class Knapsack2D:
    def __init__(self):
        self.Model = cp_model.CpModel()
        self.Solver = None

    def Solve(self, items, itemsToFix, objectiveCoefficients, bin, timeLimit = 1.0):
        n = len(items)
        H = bin.Dy
        W = bin.Dx

        model = cp_model.CpModel()

        #preprocessing
        binDomains = []
        fixItemToBin = [False] * len(items)
        for i, item in enumerate(items):
            binDomains.append([0])
            if item.Id not in itemsToFix:
                binDomains[i].extend([i + 1])
            else:
                fixItemToBin[i] = True

        # variables

        b = []
        xb1 = []
        xb2 = []
        y1 = []
        y2 = []
        for i, item in enumerate(items):

            f = model.NewBoolVar(f'b{i}')
            b.append(f)

            yStart = model.NewIntVar(0, H - item.Dy,f'y1.{i}')
            yEnd = model.NewIntVar(item.Dy, H,f'y2.{i}')

            y1.append(yStart)
            y2.append(yEnd)

            if fixItemToBin[i]:
                d = model.NewIntVar(0, bin.Dx - item.Dx, f'xb1.{i}')
                e = model.NewIntVar(item.Dx, bin.Dx, f'xb2.{i}')

                xb1.append(d)
                xb2.append(e)
            else:
                # TODO: apply bin domains to these variables
                binStart = (i + 1) * bin.Dx
                d = model.NewIntVarFromDomain(Domain.FromIntervals([[0, bin.Dx - item.Dx], [binStart, binStart]]), f'xb1.{i}')
                e = model.NewIntVarFromDomain(Domain.FromIntervals([[item.Dx, bin.Dx], [binStart + item.Dx, binStart + item.Dx]]), f'xb2.{i}')
                
                xb1.append(d)
                xb2.append(e)


        # interval variables
        xival = [model.NewIntervalVar(xb1[i], items[i].Dx, xb2[i],f'xival{i}') for i in range(n)]
        yival = [model.NewIntervalVar(y1[i], items[i].Dy, y2[i],f'yival{i}') for i in range(n)]

        # constraints
        model.AddNoOverlap2D(xival, yival)

        for i, item in enumerate(items):
            #model.AddImplication(xb2[i] <= bin.Dx, b[i])
            #model.AddImplication(xb2[i] >= bin.Dx + item.Dx, b[i].Not())
            
            #model.AddLessOrEqual(xb2[i], bin.Dx).OnlyEnforceIf(b[i])
            #model.AddGreaterOrEqual(xb2[i], bin.Dx + item.Dx).OnlyEnforceIf(b[i].Not())
            
            model.Add(xb2[i] <= bin.Dx).OnlyEnforceIf(b[i])
            model.Add(xb2[i] >= bin.Dx + item.Dx).OnlyEnforceIf(b[i].Not())
            
            #model.Add(b[i] == 1).OnlyEnforceIf(xb2[i] <= bin.Dx)
            #model.Add(b[i] == 0).OnlyEnforceIf(xb2[i] >= bin.Dx + item.Dx)
            if item.Id in itemsToFix:
                model.Add(b[i] == 1)

        #packedItems = model.NewIntVar(0, len(items), 'itemCount')
        #model.Add(packedItems == sum(b[i] * objectiveCoefficients[i] for i in range(len(items))))

        # objective
        #model.Maximize(cp_model.LinearExpr.Sum(b))  
        #model.Maximize(cp_model.LinearExpr.Sum(sum(var for var in b)))  
        #model.Maximize(cp_model.LinearExpr.BooleanSum(b))    
        #model.Maximize(packedItems)
        model.Maximize(sum(b[i] * int(objectiveCoefficients[item.Id]) for i, item in enumerate(items)))

        # solve model
        solver = cp_model.CpSolver()
        solver.parameters.log_search_progress = False
        solver.parameters.max_time_in_seconds = timeLimit
        solver.parameters.num_search_workers = 1
        
        #with open(f"Model_{0}.txt","a") as f:
        #    f.write(str(model.Proto()))

        rc = solver.Solve(model)

        return solver.StatusName(), solver.BestObjectiveBound(), solver.ObjectiveValue()

# https://www.xiang.dev/cp-sat/
class BinPacking2D:
    def __init__(self):
        self.Model = cp_model.CpModel()
        self.Solver = None

        self.StartX = []
        self.EndX = []
        self.IntervalX = []

        self.StartY = []
        self.EndY = []
        self.IntervalY = []
    
    def AddItems(self, items):
        self.Items = items

    def AddBin(self, bin):
        self.Bin = bin

    @staticmethod
    def DetermineNormalPatternsX(items, binDx, offsetX = 0):
        if binDx <= 0:
            return [0]

        X = [0] * (binDx + 1)
        X[0] = 1

        for i, item in enumerate(items):
            for p in range(binDx - item.Dx, -1, -1):
                if X[p] == 1:
                    X[p + item.Dx] = 1

        normalPatternsX = []
        for p in range (binDx, -1, -1):
            if X[p] == 1:
                normalPatternsX.append(offsetX + p)

        return normalPatternsX

    @staticmethod
    def DetermineNormalPatternsY(items, binDy):
        if binDy <= 0:
            return [0]

        Y = [0] * (binDy + 1)
        Y[0] = 1

        for i, item in enumerate(items):
            for p in range(binDy - item.Dy, -1, -1):
                if Y[p] == 1:
                    Y[p + item.Dy] = 1

        normalPatternsY = []
        for p in range (binDy, -1, -1):
            if Y[p] == 1:
                normalPatternsY.append(p)

        return normalPatternsY

    @staticmethod
    def DetermineNormalPatterns(items, binDx, binDy, offsetX = 0):
        normalPatternsX = BinPacking2D.DetermineNormalPatternsX(items, binDx)
        normalPatternsY = BinPacking2D.DetermineNormalPatternsY(items, binDy)

        return normalPatternsX, normalPatternsY

    @staticmethod
    def DetermineMeetInTheMiddlePatternsX(items, itemI, binDx, t, offsetX = 0):
        """
        itemI = items[selectedItemIndex]
        filteredItems = []
        for i, item in enumerate(items):
            if i == selectedItemIndex:
                continue
            filteredItems.append(item)
        """

        meetInTheMiddlePoints = BinPacking2D.DetermineNormalPatternsX(items, min(t - 1, binDx - itemI.Dx), offsetX) # placemenetPointsLeft
        placemenetPointsRightPrime = BinPacking2D.DetermineNormalPatternsX(items, binDx - itemI.Dx - t, offsetX)

        for p in placemenetPointsRightPrime:
            meetInTheMiddlePoints.append(binDx - itemI.Dx - p)

        return meetInTheMiddlePoints

    @staticmethod
    def DetermineMeetInTheMiddlePatternsY(items, itemI, binDy, t):
        """
        itemI = items[selectedItemIndex]
        filteredItems = []
        for i, item in enumerate(items):
            if i == selectedItemIndex:
                continue
            filteredItems.append(item)
        """

        meetInTheMiddlePoints = BinPacking2D.DetermineNormalPatternsY(items, min(t - 1, binDy - itemI.Dy)) # placemenetPointsLeft
        placemenetPointsRightPrime = BinPacking2D.DetermineNormalPatternsY(items, binDy - itemI.Dy - t)

        for p in placemenetPointsRightPrime:
            meetInTheMiddlePoints.append(binDy - itemI.Dy - p)

        return meetInTheMiddlePoints

    @staticmethod
    def DetermineMeetInTheMiddlePatterns(items, itemI, binDx, binDy, offsetX = 0):
        meetInTheMiddlePointsX = set()
        meetInTheMiddlePointsY = set()
        for t in range(1, binDx + 1):
            meetInTheMiddlePoints = BinPacking2D.DetermineMeetInTheMiddlePatternsX(items, itemI, binDx, t, offsetX)
            #meetInTheMiddlePointsX.extend(meetInTheMiddlePoints)
            meetInTheMiddlePointsX.update(meetInTheMiddlePoints)

        for t in range(1, binDy + 1):
            meetInTheMiddlePoints = BinPacking2D.DetermineMeetInTheMiddlePatternsY(items, itemI, binDy, t)
            #meetInTheMiddlePointsY.extend(meetInTheMiddlePoints)
            meetInTheMiddlePointsY.update(meetInTheMiddlePoints)

        return list(meetInTheMiddlePointsX), list(meetInTheMiddlePointsY)
        
    @staticmethod
    def DetermineMinimalMeetInTheMiddlePatternsX(items, itemI, binDx, offsetX = 0):
        meetInTheMiddlePointsLeftX = [0] * (binDx + 1)
        meetInTheMiddlePointsRightX = [0] * (binDx + 1)

        regularNormalPatterns = BinPacking2D.DetermineNormalPatternsX(items, binDx - itemI.Dx, offsetX)
        #normalPatternsX = BinPacking2D.DetermineNormalPatternsX(items, binDx - itemI.Dx, offsetX)

        #regularNormalPatterns = set(normalPatternsX)
        for p in regularNormalPatterns:
            meetInTheMiddlePointsLeftX[p] = 1 # meetInTheMiddlePointsLeftX[p] + 1
            meetInTheMiddlePointsRightX[binDx - itemI.Dx - p] = 1

        for p in range(1, binDx + 1):
            meetInTheMiddlePointsLeftX[p] = meetInTheMiddlePointsLeftX[p] + meetInTheMiddlePointsLeftX[p - 1]
            meetInTheMiddlePointsRightX[binDx - p] = meetInTheMiddlePointsRightX[binDx - p] + meetInTheMiddlePointsRightX[binDx - (p - 1)]

        tMin = 1
        minX = meetInTheMiddlePointsLeftX[0] + meetInTheMiddlePointsRightX[1]

        for p in range(2, binDx + 1):
            if meetInTheMiddlePointsLeftX[p - 1] + meetInTheMiddlePointsRightX[p] < minX:
                minX = meetInTheMiddlePointsLeftX[p - 1] + meetInTheMiddlePointsRightX[p]
                tMin = p

        meetInTheMiddlePointsX = []
        for p in regularNormalPatterns:
            if p < tMin:
                meetInTheMiddlePointsX.append(p)

            if binDx - itemI.Dx - p >= tMin:
                meetInTheMiddlePointsX.append(binDx - itemI.Dx - p)

        return meetInTheMiddlePointsX
        
    @staticmethod
    def DetermineMinimalMeetInTheMiddlePatternsY(items, itemI, binDy):
        meetInTheMiddlePointsLeftY = [0] * (binDy + 1)
        meetInTheMiddlePointsRightY = [0] * (binDy + 1)

        normalPatternsY = BinPacking2D.DetermineNormalPatternsY(items, binDy - itemI.Dy)

        regularNormalPatterns = set(normalPatternsY)
        for p in regularNormalPatterns:
            meetInTheMiddlePointsLeftY[p] = 1 # meetInTheMiddlePointsLeftX[p] + 1
            meetInTheMiddlePointsRightY[binDy - itemI.Dy - p] = 1

        for p in range(1, binDy + 1):
            meetInTheMiddlePointsLeftY[p] = meetInTheMiddlePointsLeftY[p] + meetInTheMiddlePointsLeftY[p - 1]
            meetInTheMiddlePointsRightY[binDy - p] = meetInTheMiddlePointsRightY[binDy - p] + meetInTheMiddlePointsRightY[binDy - (p - 1)]

        tMin = 1
        minY = meetInTheMiddlePointsLeftY[0] + meetInTheMiddlePointsRightY[1]

        for p in range(2, binDy + 1):
            if meetInTheMiddlePointsLeftY[p - 1] + meetInTheMiddlePointsRightY[p] < minY:
                minY = meetInTheMiddlePointsLeftY[p - 1] + meetInTheMiddlePointsRightY[p]
                tMin = p

        meetInTheMiddlePointsY = []
        for p in regularNormalPatterns:
            if p < tMin:
                meetInTheMiddlePointsY.append(p)

            if binDy - itemI.Dy - p >= tMin:
                meetInTheMiddlePointsY.append(binDy - itemI.Dy - p)

        return meetInTheMiddlePointsY
        
    @staticmethod
    def DetermineMinimalMeetInTheMiddlePatterns(items, itemI, binDx, binDy, offsetX = 0):
        meetInTheMiddlePointsX = BinPacking2D.DetermineMinimalMeetInTheMiddlePatternsX(items, itemI, binDx, offsetX)
        meetInTheMiddlePointsY = BinPacking2D.DetermineMinimalMeetInTheMiddlePatternsY(items, itemI, binDy)

        return meetInTheMiddlePointsX, meetInTheMiddlePointsY

    def CreateVariables(self, placementPointStrategy):
        binDx = self.Bin.Dx
        binDy = self.Bin.Dy
        # Consider modifying domains according to Cote, Iori (2018): The meet-in-the-middle principle for cutting and packing problems.
        for i, item in enumerate(self.Items):
        
            filteredItems = [itemJ for j, itemJ in enumerate(self.Items) if i != j]

            if placementPointStrategy == PlacementPointStrategy.UnitDiscretization:
                placementPointsX = range(0, binDx + 1 - item.Dx)
                placementPointsY = range(0, binDy + 1 - item.Dy)
            elif placementPointStrategy == PlacementPointStrategy.NormalPatterns:
                placementPointsX, placementPointsY = BinPacking2D.DetermineNormalPatterns(filteredItems, binDx - item.Dx, binDy - item.Dy)
            elif placementPointStrategy == PlacementPointStrategy.MeetInTheMiddlePatterns:
                placementPointsX, placementPointsY = BinPacking2D.DetermineMeetInTheMiddlePatterns(filteredItems, item, binDx, binDy)
            elif placementPointStrategy == PlacementPointStrategy.MinimalMeetInTheMiddlePatterns:
                placementPointsX, placementPointsY = BinPacking2D.DetermineMinimalMeetInTheMiddlePatterns(filteredItems, item, binDx, binDy)
            else:
                raise ValueError('UnkownPlacementPointStrategy')

            placementPointsStartX = [p for p in placementPointsX if p + item.Dx <= binDx] # unnecessary, is satisfied by construction
            placementPointsEndX = [p + item.Dx for p in placementPointsStartX]
            x1 = self.Model.NewIntVarFromDomain(Domain.FromValues(placementPointsStartX), f'x1.{i}')
            x2 = self.Model.NewIntVarFromDomain(Domain.FromValues(placementPointsEndX), f'x2.{i}')

            self.StartX.append(x1)
            self.EndX.append(x2)

            placementPointsStartY = [p for p in placementPointsY if p + item.Dy <= binDy]  # is satisfied by construction
            placementPointsEndY = [p + item.Dy for p in placementPointsStartY]
            y1 = self.Model.NewIntVarFromDomain(Domain.FromValues(placementPointsStartY), f'y1.{i}')
            y2 = self.Model.NewIntVarFromDomain(Domain.FromValues(placementPointsEndY), f'y2.{i}')

            self.StartY.append(y1)
            self.EndY.append(y2)

            intervalX = self.Model.NewIntervalVar(x1, item.Dx, x2, f'xI{i}')
            intervalY = self.Model.NewIntervalVar(y1, item.Dy, y2, f'yI{i}')

            self.IntervalX.append(intervalX)
            self.IntervalY.append(intervalY)

    def CreateConstraints(self):
        """for i, itemI in enumerate(self.Items):
            for j, itemJ in enumerate(self.Items):
                if i == j:
                    continue

                self.Model.AddNoOverlap2D([self.IntervalX[i], self.IntervalX[j]],[self.IntervalY[i], self.IntervalY[j]])
                """

        self.Model.AddNoOverlap2D(self.IntervalX, self.IntervalY)

    def Solve(self):
        if len(self.Items) == 1:
            if self.Items[0].Dx <= self.Bin.Dx and self.Items[0].Dy <= self.Bin.Dy:
                return True
            else:
                return False

        self.Solver = cp_model.CpSolver()
        self.Solver.parameters.log_search_progress = False 
        self.Solver.parameters.num_search_workers = 8
        #self.Solver.parameters.max_time_in_seconds = 10
        #solver.parameters.cp_model_presolve = False

        """
        numberOfItems = len(self.Items)
        with open(f"Model_{numberOfItems}.txt","a") as f:
            f.write(str(self.Model.Proto()))
        """

        status = self.Solver.Solve(self.Model)

        if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
            return True
        elif status == cp_model.INFEASIBLE:
            return False
        elif status == cp_model.MODEL_INVALID:
            raise ValueError("Model invalid")
        elif status == cp_model.UNKNOWN:
            raise ValueError("Unkown status type")
        else:
            raise ValueError("Invalid status type")

        print(status)
        print(self.Solver.StatusName())

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
            if i in infeasibleItemSubset or i < binId or model._FixItemToBin[i]:
                continue

            skipItem = False
            for j in infeasibleItemSubset:
                #isIncompatible = frozenset(i, j) in model._IncompatibleItems
                if model._FixItemToBin[j] and frozenset((i, j)) in model._IncompatibleItems:
                    skipItem = True
                    break
            
            if skipItem:
                continue

            # check against _FeasibleSets or _InfeasibleSets

            upperBound = BinPackingCallback.SolveKnapsack2D(liftedSet, i, model._FixItemToBin, model._Bins[binId], objectiveCoefficients, model, 1.0)
            
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
        
        if model._FixItemToBin[binId]:
            BinPackingCallback.AddLiftedCoverInequality(infeasibleItemSubset, additionalItems, liftingCoefficients, itemVariables, [binId], model)
        else:
            filteredBins = [b for b, bin in enumerate(model._Bins) if b <= minItemId and b < upperBoundBin]# and not model._FixItemToBin[b]]
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

            if frozenset(itemIndices) in model._InfeasibleSets or frozenset(itemIndices) in model._IncompatibleItems:
                #raise ValueError('No previously discovered infeasible subsets should be found during cut strengthening.')
                model._InfeasibleDoubleCount += 1
                #print(f"Double count in bin {binId} (fixed = {model._FixItemToBin[binId]}): {itemIndices} vs. {[item.Id for i, item in enumerate(items)]}")
                continue
                #BinPackingCallback.AddCut(itemIndices, model._VarsX, model._Bins, model)
                #return

            t1 = time.time()

            binPacking2D = BinPacking2D()
            binPacking2D.AddBin(model._Bins[0]) # homoegeneous bins
            binPacking2D.AddItems(sortedItems)

            binPacking2D.CreateVariables(model._PlacementPointStrategy)
            binPacking2D.CreateConstraints()

            isFeasible = binPacking2D.Solve()

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

            if frozenset(itemIndices) in model._InfeasibleSets:# or frozenset(itemIndices) in model._IncompatibleItems:
                #raise ValueError('Infeasible items occured twice')
                #print(f"Double count in bin {b} (fixed = {model._FixItemToBin[b]}): {itemIndices}")
                model._InfeasibleDoubleCount += 1
                infeasibleSubproblems += 1 # only count new cuts
                BinPackingCallback.AddCut(itemIndices, model._VarsX, [b], model)
                continue

            t1 = time.time()

            binPacking2D = BinPacking2D()
            binPacking2D.AddBin(model._Bins[b]) # homoegeneous bins
            binPacking2D.AddItems(itemsInBin)

            binPacking2D.CreateVariables(model._PlacementPointStrategy)
            binPacking2D.CreateConstraints()

            isFeasible = binPacking2D.Solve()

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

        self.LowerBoundBin = 0
        self.FixItemToBin = []
        self.IncompatibleItems = set()

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

            # symmetry breaking, only valid if bins homogeneous
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
        self.Model._IncompatibleItems = self.IncompatibleItems
        self.Model._FixItemToBin = self.FixItemToBin

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
        self.Model._EnableCutStrengthening = True
        self.Model._InfeasibleSuproblemCutThreshold = 1
        
        self.Model._EnablePreprocessLifting = False
        self.Model._PlacementPointStrategy = PlacementPointStrategy.MinimalMeetInTheMiddlePatterns

    def DeterminePositions(self, itemIndicesArray, itemsInBinArray):
        rectanglesArray = []
        for b, itemIndices in enumerate(itemIndicesArray):
            numberOfItems = len(itemIndices)
            if numberOfItems == 0:
                continue

            itemsInBin = itemsInBinArray[b]

            bin = self.Model._Bins[b]

            binPacking2D = BinPacking2D()
            binPacking2D.AddBin(bin) # homoegeneous bins
            binPacking2D.AddItems(itemsInBin)

            binPacking2D.CreateVariables(self.Model._PlacementPointStrategy)
            binPacking2D.CreateConstraints()

            isFeasible = binPacking2D.Solve()

            if not isFeasible:
                raise ValueError('Packing must not be infeasible during solution extraction.')

            #binPacking2D.StartX
            
            xArray = []
            yArray = []
            if numberOfItems == 1:
                xArray = [(b * bin.Dx) + 0.0]
                yArray = [0.0]
            else:
                xArray = [(b * bin.Dx) + binPacking2D.Solver.Value(binPacking2D.StartX[i]) for i in range(numberOfItems)]
                yArray = [binPacking2D.Solver.Value(binPacking2D.StartY[i]) for i in range(numberOfItems)]

            w = []
            h = []
            for item in itemsInBin:
                w.append(item.Dx)
                h.append(item.Dy)

            rectangles = ExtractDataForPlot(xArray, yArray, w, h, bin.Dx, bin.Dy)

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

    def FixIncompatibleItems(self, incompatibleItems):
        numberOfItems = len(self.Items)
        fixItemToBin = [False] * len(self.Items)

        # This is a similar logic as in section 6.2 in Cote, Haouari, Iori (2019). 
        self.ItemVariables[0][0].LB = 1.0
        fixItemToBin[0] = True

        for i in range(1, numberOfItems):
            itemI = self.Items[i]

            isIncompatible = True
            for j in range(0, i):
                if frozenset((i, j)) not in incompatibleItems:
                    isIncompatible = False
                    return fixItemToBin
            
            if isIncompatible:
                fixItemToBin[i] = True
                for b, bin in enumerate(self.Bins):
                    if b == i:
                        self.ItemVariables[b][i].LB = 1.0
                    else:
                        self.ItemVariables[b][i].UB = 0.0

        return fixItemToBin

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
                if i in infeasibleItemSubset or i < b or (model._FixItemToBin[i] and b != i):# or i < minIdInitial: or model._FixItemToBin[i] 
                    continue

                skipItem = False
                for j in liftedSet:
                    #isIncompatible = frozenset(i, j) in model._IncompatibleItems
                    if (model._FixItemToBin[j] and frozenset((i, j)) in model._IncompatibleItems):
                        skipItem = True
                        break
                
                if skipItem:
                    continue

                upperBound = BinPackingCallback.SolveKnapsack2D(liftedSet, i, model._FixItemToBin, model._Bins[0], objectiveCoefficients, model, 1.0)
                
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

    def ApplyPreprocessChanges(self, incompatibleItems):
        fixItemToBin = self.FixIncompatibleItems(incompatibleItems)

        self.FixItemToBin = fixItemToBin
        self.IncompatibleItems = incompatibleItems

        self.SetCallbackData()

        self.AddIncompatibilityCuts(incompatibleItems, self.FixItemToBin)
        
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

        #xArray = [solver.Value(xb1[i]) for i in range(n)]
        #yArray = [solver.Value(y1[i]) for i in range(n)]

        #rectangles = ExtractDataForPlot(xArray, yArray, w, h, W, H)

        return rectangles

class BinPackingBranchAndCutSolver:
    def __init__(self):
        self.BinPacking = BinPackingMip()

        self.RemovedItems = []
        self.IncompatibleItems = set()

        self.IsOptimal = False
        self.LB = -1
        self.UB = -1
        self.Runtime = -1
        self.SolverType = ""
        self.CutCount = -1

    def RetrieveSolutionStatistics(self):
        if self.IsOptimal:
            return 

        model = self.BinPacking.Model

        self.IsOptimal = 1 if model.Status == GRB.OPTIMAL else 0
        self.LB = model.objBound
        self.UB = model.objVal
        self.Runtime = model.Runtime
        self.SolverType = "B&C"

    def DetermineStartSolution(self, items, H, W, lowerBoundBin, m):
        solverCP = BinPackingSolverCP()
        
        h = []
        w = []
        for i, item in enumerate(items):
            h.append(item.Dy)
            w.append(item.Dx)

        rectangles = solverCP.BinPackingErwin(items, h, w, H, W, lowerBoundBin, m, 30, False, self.IncompatibleItems)

        if solverCP.LB == solverCP.UB:
            return True, solverCP.LB, "CP", rectangles

        self.BinPacking.SetStartSolution(solverCP.ItemBinAssignments, solverCP.LB, solverCP.UB)

        return False, solverCP.LB, "B&C", rectangles

    def RemoveLargeItems(self, items, H, W, m):
        filteredItemIndices = []
        for i, item in enumerate(items):
            dy = item.Dy
            dx = item.Dx

            if dy == H and dx == W:
                #print(f'Item {i} has the same dimensions as the bin and will be removed.')
                self.RemovedItems.append(Item(i, dx, dy))
                continue
            
            isFullyIncompatible = True
            for j, itemJ in enumerate(items):
                if i == j:
                    continue

                if item.Dx + itemJ.Dx > W and item.Dy + itemJ.Dy > H:
                    continue

                isFullyIncompatible = False
                break

            if isFullyIncompatible:
                #print(f'Item {i} is fully incompatible and will be removed.')
                self.RemovedItems.append(Item(i, dx, dy))
                continue

            filteredItemIndices.append(i)

        self.BinPacking.Model.ObjCon = len(self.RemovedItems)

        newItems = []
        for index, i in enumerate(filteredItemIndices):
            newItems.append(Item(index, items[i].Dx, items[i].Dy))

        return newItems, len(filteredItemIndices)

    def DetermineConflicts(self, items, H, W):
        for i, itemI in enumerate(items):
            for j in range(i + 1, len(items)):
                itemJ = items[j]

                if itemI.Dx + itemJ.Dx > W and itemI.Dy + itemJ.Dy > H:
                    self.IncompatibleItems.add(frozenset((i, j)))
                    continue

    def Preprocess(self, items, H, W, m):
        items = sorted(items, reverse=True) # TODO: build conflict graph and compute maximal clique
        items, m = self.RemoveLargeItems(items, H, W, m)

        self.DetermineConflicts(items, H, W)

        conflictGraph = nx.Graph()
        conflictGraph.add_nodes_from([item.Id for i, item in enumerate(items)])
        for i in range(len(items)):
            itemI = items[i]
            for j in range(len(items)):
                itemJ = items[j]
                if frozenset((itemI.Id, itemJ.Id)) in self.IncompatibleItems:
                    conflictGraph.add_edge(itemI.Id, itemJ.Id)

        maxClique = clique.max_clique(conflictGraph)
        sortedMaxClique = sorted(maxClique)

        newItems = []
        for i, oldIndex in enumerate(sortedMaxClique):
            newItems.append(Item(i, items[oldIndex].Dx, items[oldIndex].Dy))
        
        for i, item in enumerate(items):
            if item.Id in sortedMaxClique:
                continue
            newItems.append(Item(len(newItems), items[item.Id].Dx, items[item.Id].Dy))
            
        items = newItems

        self.IncompatibleItems.clear()
        self.DetermineConflicts(items, H, W)

        return items, m

    def Run(self, h, w, H, W, m):   
        items = []
        for i in range(len(h)):
            items.append(Item(i, w[i], h[i]))

        items, m = self.Preprocess(items, H, W, m)

        self.BinPacking.AddItems(items)
        self.BinPacking.AddBins([W] * m, [H] * m)

        self.BinPacking.CreateVariables()
        self.BinPacking.CreateConstraints()

        self.BinPacking.ApplyPreprocessChanges(self.IncompatibleItems)
        
        #self.BinPacking.Model.write('BPP.lp')

        isOptimalCP, lbCP, solverType, rectangles = self.DetermineStartSolution(items, H, W, self.BinPacking.LowerBoundBin, m)

        if False:
            self.IsOptimal = 1
            self.LB = lbCP + len(self.RemovedItems)
            self.UB = lbCP + len(self.RemovedItems)
            self.SolverType = solverType
            self.Runtime = -1

            return rectangles

        rectangles = self.BinPacking.Solve()

        return rectangles

#h, w, H, W, m = ReadBenchmarkData(10)
#Run(h, w, H, W, m)

def main():
    #h, w, H, W, m = ReadExampleData()
    solutions = {}
    # Single bin 2D-BPP CP model takes ages to prove feasibility/infeasibility on instances: 173, 262, 292, 297, 298, 317
    # Meet in the middle produces erroneous 1D KP instances for instances 11
    # Postsolve error: 51, 45, 125/126, 31, 26 or 27, 311
    # Double count: 21, 22, 26, 27, 31, 32
    # Negative lifting coefficient: 120, 123, 
    #hardInstances = [226, 232, 242, 242, 244, 245, 247, 248, 249, 261, 292, 313, 314, 332, 173, 187, 188, 191, 197, 142, 143, 145, 146, 149]
    #mediumInstance = [149, 174]
    for instance in range(1, 174):
    #for instance in hardInstances:
        h, w, H, W, m = ReadBenchmarkData(instance)
        
        solver = BinPackingBranchAndCutSolver()
        rectangles = solver.Run(h, w, H, W, m)

        solver.RetrieveSolutionStatistics()

        # TODO: introduce solution statistics struct
        bestBoundMIP = solver.LB
        upperBoundMIP = solver.UB
        solverType = solver.SolverType
        isOptimalMIP = solver.IsOptimal
        
        #PlotSolution(upperBoundMIP * W, H, rectangles)

        if isOptimalMIP:
            print(f'Instance {instance}: Optimal solution = {int(bestBoundMIP)} found by {solverType} (#items = {len(h)})')
        else:
            #raise ValueError(f'Instance {instance}: No optimal solution found, [lb, ub] = [{bestBoundMIP}, {upperBoundMIP}]')
            print(f'Instance {instance}: No optimal solution found, [lb, ub] = [{bestBoundMIP}, {upperBoundMIP}] (#items = {len(h)})')
    
        solutions[instance] = {'LB': bestBoundMIP, 'UB': upperBoundMIP, 'Solver': solverType}

    solutionsJson = json.dumps(solutions, indent = 4)
    
    # Writing to sample.json
    with open("Solutions.json", "w") as outfile:
        outfile.write(solutionsJson)



if __name__ == "__main__":
    main()