#!/usr/bin/env python3.7

import gurobipy as gp
from gurobipy import GRB

import math
import numpy
import json

from ortools.sat.python import cp_model
from ortools.sat.python.cp_model import Domain

from BinPackingData import *
from ErwinCP import *

""" Datasets at https://github.com/Oscar-Oliveira/OR-Datasets/tree/master/Cutting-and-Packing/2D/Datasets """

class Item:
    def __init__(self, dx, dy):
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

# https://www.xiang.dev/cp-sat/
class BinPacking2D:
    def __init__(self):
        self.Items = []
        self.Bin = None

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
    def DetermineNormalPatterns(items, bin, offsetX = 0):
        X = [0] * bin.Dx
        Y = [0] * bin.Dy
        X[0] = 1
        Y[0] = 1

        minDx = bin.Dx
        minDy = bin.Dy

        for i, item in enumerate(items):
            minDx = min(minDx, item.Dx)
            minDy = min(minDy, item.Dy)
            for p in range(bin.Dx - item.Dx, -1, -1):
                if X[p] == 1:
                    if p + item.Dx < bin.Dx:
                        X[p + item.Dx] = 1

            for p in range(bin.Dy - item.Dy, -1, -1):
                if Y[p] == 1:
                    if p + item.Dy < bin.Dy:
                        Y[p + item.Dy] = 1

        normalPatternsX = []
        for p in range (bin.Dx - minDx, -1, -1):
            if X[p] == 1:
                normalPatternsX.append(offsetX + p)

        normalPatternsY = []
        for p in range (bin.Dy - minDy, -1, -1):
            if Y[p] == 1:
                normalPatternsY.append(p)

        return normalPatternsX, normalPatternsY

    def CreateVariables(self):
        normalPatternsX, normalPatternsY = BinPacking2D.DetermineNormalPatterns(self.Items, self.Bin)

        # Consider modifying domains according to Cote, Iori (2018): The meet-in-the-middle principle for cutting and packing problems.
        for i, item in enumerate(self.Items):

            normalPatternsStartX = [p for p in normalPatternsX if p + item.Dx <= self.Bin.Dx]
            normalPatternsEndX = [p + item.Dx for p in normalPatternsStartX]
            x1 = self.Model.NewIntVarFromDomain(Domain.FromValues(normalPatternsStartX), f'x1.{i}')
            x2 = self.Model.NewIntVarFromDomain(Domain.FromValues(normalPatternsEndX), f'x2.{i}')

            """
            x1 = self.Model.NewIntVar(0, self.Bin.Dx - item.Dx, f'x1.{i}')
            x2 = self.Model.NewIntVar(item.Dx, self.Bin.Dx, f'x2.{i}')
            """

            self.StartX.append(x1)
            self.EndX.append(x2)

            normalPatternsStartY = [p for p in normalPatternsY if p + item.Dy <= self.Bin.Dy]
            normalPatternsEndY = [p + item.Dy for p in normalPatternsStartY]
            y1 = self.Model.NewIntVarFromDomain(Domain.FromValues(normalPatternsStartY), f'y1.{i}')
            y2 = self.Model.NewIntVarFromDomain(Domain.FromValues(normalPatternsEndY), f'y2.{i}')

            """
            y1 = self.Model.NewIntVar(0, self.Bin.Dy - item.Dy, f'y1.{i}')
            y2 = self.Model.NewIntVar(item.Dy, self.Bin.Dy, f'y2.{i}')
            """

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
        self.Solver.parameters.num_search_workers = 1
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
        for b, bin in enumerate(bins):
            expr = gp.LinExpr()
            for i in items:
                expr += itemVariables[b][i]

            model.cbLazy(expr <= len(items) - 1)

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
        for b, itemIndices in enumerate(itemIndicesArray):
            itemsInBin = itemsInBinArray[b]

            if len(itemIndices) == 0:
                continue

            if frozenset(itemIndices) in model._FeasibleSets:
                continue

            if frozenset(itemIndices) in model._InfeasibleSets:
                #raise ValueError('Infeasible items occured twice')
                #print("Double Count")
                model._InfeasibleDoubleCount += 1
                BinPackingCallback.AddCut(itemIndices, model._VarsX, model._Bins, model)
                continue

            binPacking2D = BinPacking2D()
            binPacking2D.AddBin(model._Bins[b]) # homoegeneous bins
            binPacking2D.AddItems(itemsInBin)

            binPacking2D.CreateVariables()
            binPacking2D.CreateConstraints()

            isFeasible = binPacking2D.Solve()

            if isFeasible:
                model._FeasibleSets.add(frozenset(itemIndices))
            else:
                model._InfeasibleSets.add(frozenset(itemIndices))
                model._CutCount += 1
                BinPackingCallback.AddCut(itemIndices, model._VarsX, model._Bins, model)

    @staticmethod
    def callback(model, where):
        if where == GRB.Callback.MIPSOL:
            itemIndicesArray, itemsInBinArray = BinPackingCallback.FindIntegerAssignments(model)
            BinPackingCallback.AddCuts(model, itemIndicesArray, itemsInBinArray)

class BinPackingMip:
    def __init__(self, enable2D = True):
        self.Model = gp.Model("BinPacking")

        self.Model.Params.OutputFlag = 1
        self.Model.Params.lazyConstraints = 1
        self.Model.Params.MIPFocus = 2
        self.Model.Params.TimeLimit = 480 if enable2D else 30

        self.Items = []
        self.Bins = []

        self.ItemVariables = []
        self.BinVariables = []

        self.Callback = BinPackingCallback()

        self.Enable2D = enable2D

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

        lowerBoundBin = max(lowerBoundBinArea, lb1D)

        for b in range(lowerBoundBin):
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
        self.Model._FeasibleSets = set()
        self.Model._InfeasibleSets = set()
        self.Model._CutCount = 0
        self.Model._InfeasibleDoubleCount = 0

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

            binPacking2D.CreateVariables()
            binPacking2D.CreateConstraints()

            isFeasible = binPacking2D.Solve()

            if not isFeasible:
                raise ValueError('Packing must not be infeasible during solution extraction.')

            binPacking2D.StartX
            
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

        # This is a similar logic as in section 6.2 in Cote, Haouari, Iori (2019). 
        self.ItemVariables[0][0].LB = 1.0

        for i in range(1, numberOfItems):
            itemI = self.Items[i]

            isIncompatible = True
            for j in range(0, i):
                if frozenset((i, j)) not in incompatibleItems:
                    isIncompatible = False
                    return i
            
            if isIncompatible:
                for b, bin in enumerate(self.Bins):
                    if b == i:
                        self.ItemVariables[b][i].LB = 1.0
                    else:
                        self.ItemVariables[b][i].UB = 0.0

        return numberOfItems - 1

    def AddIncompatibilityCuts(self, incompatibleItems, lastFixedItemIndex):
        for i, j in incompatibleItems:
            if i < lastFixedItemIndex and j < lastFixedItemIndex:
                continue

            for b, bin in enumerate(self.Bins):
                # are the same cuts added multiple times if contained multiple times in incompatibleItems?
                # TODO.Logic: lifting 
                self.Model.addConstr(self.ItemVariables[b][i] + self.ItemVariables[b][j] <= 1)

    def ApplyPreprocessChanges(self, incompatibleItems):
        lastFixedItemIndex = self.FixIncompatibleItems(incompatibleItems)
        self.AddIncompatibilityCuts(incompatibleItems, lastFixedItemIndex)
        

    def Solve(self):
        self.SetCallbackData()

        if self.Enable2D:
            self.Model.optimize(BinPackingCallback.callback)
        else:
            self.Model.optimize()

        #self.Model.printAttr('x')

        #self.Model.write('BPP.lp')

        statusCode = 1 if self.Model.Status == GRB.OPTIMAL else 0

        if self.Model._InfeasibleDoubleCount > 0:
            #print(f'Infeasible double count: {self.Model._InfeasibleDoubleCount}')
            pass
        
        rectangles = None
        if self.Enable2D:
            rectangles = self.ExtractSolution()

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

    def DetermineStartSolution(self, items, H, W, m):
        solverCP = BinPackingSolverCP()
        
        h = []
        w = []
        for i, item in enumerate(items):
            h.append(item.Dy)
            w.append(item.Dx)

        rectangles = solverCP.BinPackingErwin(items, h, w, H, W, m, 10, False, self.IncompatibleItems)

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
                self.RemovedItems.append(Item(dx, dy))
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
                self.RemovedItems.append(Item(dx, dy))
                continue

            filteredItemIndices.append(i)

        self.BinPacking.Model.ObjCon = len(self.RemovedItems)

        newItems = []
        for i in filteredItemIndices:
            newItems.append(Item(items[i].Dx, items[i].Dy))

        return newItems, len(filteredItemIndices)

    def DetermineConflicts(self, items, H, W):
        for i, itemI in enumerate(items):
            for j in range(i + 1, len(items)):
                itemJ = items[j]

                if itemI.Dx + itemJ.Dx > W and itemI.Dy + itemJ.Dy > H:
                    self.IncompatibleItems.add(frozenset((i, j)))
                    self.IncompatibleItems.add(frozenset((j, i)))
                    continue

    def Preprocess(self, items, H, W, m):
        items = sorted(items, reverse=True)
        items, m = self.RemoveLargeItems(items, H, W, m)

        self.DetermineConflicts(items, H, W)

        return items, m

    def Run(self, h, w, H, W, m):   
        items = []
        for i in range(len(h)):
            items.append(Item(w[i], h[i]))

        items, m = self.Preprocess(items, H, W, m)

        self.BinPacking.AddItems(items)
        self.BinPacking.AddBins([W] * m, [H] * m)

        self.BinPacking.CreateVariables()
        self.BinPacking.CreateConstraints()

        self.BinPacking.ApplyPreprocessChanges(self.IncompatibleItems)
        
        #self.BinPacking.Model.write('BPP.lp')

        isOptimalCP, lbCP, solverType, rectangles = self.DetermineStartSolution(items, H, W, m)

        if isOptimalCP:
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
    # Single bin 2D-BPP CP model cannot prove feasibility/infeasibility on instances:
    # - 173 (strange behavior: initial num_bool: 1 even though there are no bools in the model)
    #hardInstances = [173, 146, 292, 332]
    for instance in range(332, 500):
    #for instance in hardInstances:
        h, w, H, W, m = ReadBenchmarkData(instance)
        
        solver = BinPackingBranchAndCutSolver()
        rectangles = solver.Run(h, w, H, W, m)

        solver.RetrieveSolutionStatistics()

        # TODO: introduce solution statistics strct
        bestBoundMIP = solver.LB
        upperBoundMIP = solver.UB
        solverType = solver.SolverType
        isOptimalMIP = solver.IsOptimal
        
        PlotSolution(upperBoundMIP * W, H, rectangles)

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