from ortools.sat.python import cp_model
from ortools.sat.python.cp_model import Domain

import math
import pandas
import numpy

from BinPackingData import *

from Model import *
from PlacementPoints import *

import matplotlib

from SymmetryBreaking import SymmetryBreaking

class BinPackingSolverCP:
    def __init__(self, items, binDy, binDx, lowerBoundBins, upperBoundBins, timeLimit = 3600, enableLogging = True, incompatibleItems = None):
        self.items = items
        self.binDx = binDx
        self.binDy = binDy
        self.lowerBoundBins = lowerBoundBins
        self.upperBoundBins = upperBoundBins
        self.timeLimit = timeLimit
        self.enableLogging = enableLogging
        self.incompatibleItems = incompatibleItems
        
        self.IsOptimal = False
        self.LB = -1
        self.UB = -1

        self.EnableNormalPatterns = False

        self.ItemBinAssignments = []
    
    def Solve(self, modelType = 'OneBigBin'):
        if modelType == 'OneBigBin':
            model = OneBigBinModel()
        else:
            raise ValueError("Invalid bin packing model type.")
        
        rectangles = model.Solve(self.items, self.binDy, self.binDx, self.lowerBoundBins, self.upperBoundBins, self.timeLimit, self.enableLogging, self.incompatibleItems)
       
        self.IsOptimal = model.IsOptimal
        self.LB = model.LB
        self.UB = model.UB
        self.ItemBinAssignments = model.ItemBinAssignments

        return rectangles

class OneBigBinModel:
    def __init__(self):
        self.IsOptimal = False
        self.LB = -1
        self.UB = -1

        self.EnableNormalPatterns = False

        self.ItemBinAssignments = []

    def FixIncompatibleItems(self, incompatibleItems, numberOfItems):
        # This is a similar logic as in section 6.2 in Cote, Haouari, Iori (2019). 
        fixItemToBin = [False] * numberOfItems
        fixItemToBin[0] = True

        if incompatibleItems == None or len(incompatibleItems) == 0:
            return fixItemToBin

        for i in range(1, numberOfItems):
            isIncompatible = True
            for j in range(0, i):
                if frozenset((i, j)) not in incompatibleItems:
                    isIncompatible = False
                    return fixItemToBin
            
            if isIncompatible:
                fixItemToBin[i] = True

        return fixItemToBin

    """ For reference, see section 6.1 from Cote, Iori (2018). """
    def CreateReducedBinDomains(self, incompatibleItems, numberOfItems, numberOfBins, fixItemToBin):
        fullDomains = []
        fullDomains.append([0])
        if incompatibleItems == None or len(incompatibleItems) == 0:
            for i in range(1, numberOfItems):
                boundedNumberOfBins = i if i < numberOfBins else numberOfBins - 1
                fullDomains.append([j for j in range(boundedNumberOfBins + 1)])

            return fullDomains

        
        for i in range(1, numberOfItems):
            boundedNumberOfBins = i if i < numberOfBins else numberOfBins - 1
            fullDomains.append([boundedNumberOfBins])

            if fixItemToBin[i]:
                continue

            for j in range(0, numberOfItems):
                if j >= boundedNumberOfBins:
                    break

                if fixItemToBin[j] and frozenset((i, j)) in incompatibleItems:
                    continue
                else:
                    fullDomains[i].append(j)

        return fullDomains

    def CreateBinDependentNormalPatterns(self, incompatibleItems, fixItemToBin, items, numberOfBins, binDx, binDy):
        bin = Bin(binDx, binDy)
        numberOfItems = len(items)
        
        itemSpecificNormalPatternsX = []
        itemSpecificNormalPatternsY = []

        # Determine placement points for specific items depending on the combination with every other compatible item.
        for i, itemI in enumerate(items):
            itemSubset = []
            for j, itemJ in enumerate(items):
                if (fixItemToBin[i] and fixItemToBin[j]) or i == j or frozenset((i, j)) in incompatibleItems:
                    continue

                itemSubset.append(itemJ)

            normalPatternsX, normalPatternsY = PlacementPointGenerator.DetermineNormalPatterns(itemSubset, bin.Dx - itemI.Dx, bin.Dy - itemI.Dy)
            itemSpecificNormalPatternsX.append(normalPatternsX)
            itemSpecificNormalPatternsY.append(normalPatternsY)

        # Determine placement points for items in specific bins and all other compatible items.
        itemBinNormalPatternsX = []
        for i, itemI in enumerate(items):
            itemBinNormalPatternsX.append([])

            if fixItemToBin[i]:
                itemSubset = []
                for j, itemJ in enumerate(items):
                    if fixItemToBin[j] or j <= i or frozenset((i, j)) in incompatibleItems:
                        continue
                    itemSubset.append(itemJ)

                fixedBinNormalPatternsX, fixedBinNormalPatternsY = PlacementPointGenerator.DetermineNormalPatterns(itemSubset, bin.Dx - itemI.Dx, bin.Dy - itemI.Dy, i*bin.Dx)
                itemBinNormalPatternsX[i].extend(fixedBinNormalPatternsX)

                continue

            for b in range(numberOfBins):
                if b > i:
                    break

                itemSubset = []
                for j, itemJ in enumerate(items):
                    if i == j or b > j or frozenset((i, j)) in incompatibleItems:
                        continue
                    
                    itemSubset.append(itemJ)

                binSpecificNormalPatternsX, binSpecificNormalPatternsY = PlacementPointGenerator.DetermineNormalPatterns(itemSubset, bin.Dx - itemI.Dx, bin.Dy - itemI.Dy, b*bin.Dx)
                itemBinNormalPatternsX[i].extend(binSpecificNormalPatternsX)

        return itemSpecificNormalPatternsX, itemSpecificNormalPatternsY, itemBinNormalPatternsX

    def AddIncompatibilityCuts(self, incompatibleItems, fixItemToBin, model, binVariables):
        if incompatibleItems == None:
            return

        for i, j in incompatibleItems:
            if fixItemToBin[i] and fixItemToBin[j]:
                continue
            
            model.Add(binVariables[i] != binVariables[j])

    """ 
    A simpler version of this model can be found at
    https://yetanothermathprogrammingconsultant.blogspot.com/2021/02/2d-bin-packing-with-google-or-tools-cp.html 
    https://yetanothermathprogrammingconsultant.blogspot.com/2021/02/2d-bin-packing.html 
    """
    def Solve(self, items, binDy, binDx, lowerBoundBins, upperBoundBins, timeLimit = 3600, enableLogging = True, incompatibleItems = None):

        n = len(items)

        model = cp_model.CpModel()

        #preprocessing
        fixItemToBin = self.FixIncompatibleItems(incompatibleItems, n)

        binDomains = self.CreateReducedBinDomains(incompatibleItems, n, upperBoundBins, fixItemToBin)

        itemNormalPatternsX, itemNormalPatternsY, globalNormalPatternsX = [], [], []
        if self.EnableNormalPatterns:
            itemNormalPatternsX, itemNormalPatternsY, globalNormalPatternsX = self.CreateBinDependentNormalPatterns(incompatibleItems, fixItemToBin, items, upperBoundBins, binDx, binDy)
        
        # variables

        startVariablesLocalX = [] 
        placedBinVariables = [] # bin numbers
        startVariablesGlobalX = []
        endVariablesGlobalX = []
        startVariablesY = []
        endVariablesY = []
        totalArea = 0.0
        for i, item in enumerate(items):

            totalArea += item.Dx * item.Dy

            itemFeasibleBins = model.NewIntVarFromDomain(Domain.FromValues(binDomains[i]), f'b{i}')
            placedBinVariables.append(itemFeasibleBins)

            if self.EnableNormalPatterns:
                filteredStartLocalX = [p for p in itemNormalPatternsX[i] if (p % binDx) + item.Dx <= binDx]
                xStart = model.NewIntVarFromDomain(Domain.FromValues(filteredStartLocalX),f'x{i}')
                startVariablesLocalX.append(xStart)

                filteredStartX = [p for p in globalNormalPatternsX[i] if (p % binDx) + item.Dx <= binDx]
                filteredEndX = [p + item.Dx for p in filteredStartX]
                globalStartX = model.NewIntVarFromDomain(Domain.FromValues(filteredStartX), f'xb1.{i}')
                globalEndX = model.NewIntVarFromDomain(Domain.FromValues(filteredEndX), f'xb2.{i}')

                startVariablesGlobalX.append(globalStartX)
                endVariablesGlobalX.append(globalEndX)

                filteredStartY = [p for p in itemNormalPatternsY[i] if p + item.Dy <= binDy]
                filteredEndY = [p + item.Dy for p in filteredStartY]
                yStart = model.NewIntVarFromDomain(Domain.FromValues(filteredStartY),f'y1.{i}')
                yEnd = model.NewIntVarFromDomain(Domain.FromValues(filteredEndY),f'y2.{i}')

                startVariablesY.append(yStart)
                endVariablesY.append(yEnd)
            else:
                xStart = model.NewIntVar(0, binDx - item.Dx, f'x{i}')
                startVariablesLocalX.append(xStart)

                if fixItemToBin[i]:
                    # do domain reduction
                    reducedDomainX = SymmetryBreaking.ReducedDomainX(binDx, item)
                    globalStartX = model.NewIntVar(i*binDx, i*binDx + reducedDomainX, f'xb1.{i}')
                    globalEndX = model.NewIntVar(i*binDx + item.Dx, i*binDx + reducedDomainX + item.Dx, f'xb2.{i}')
                    #d = model.NewIntVar(i*W, (i + 1)*W - item.Dx, f'xb1.{i}')
                    #e = model.NewIntVar(i*W + item.Dx, (i + 1)*W, f'xb2.{i}')

                    startVariablesGlobalX.append(globalStartX)
                    endVariablesGlobalX.append(globalEndX)
                    
                    reducedDomainY = SymmetryBreaking.ReducedDomainY(binDy, item)
                    yStart = model.NewIntVar(0, reducedDomainY,f'y1.{i}')
                    yEnd = model.NewIntVar(item.Dy, reducedDomainY + item.Dy,f'y2.{i}')

                    startVariablesY.append(yStart)
                    endVariablesY.append(yEnd)
                else:
                    # TODO: domain reduction for each bin where item i is the biggest placeable item
                    boundedM = i if i < upperBoundBins else upperBoundBins - 1

                    # TODO: apply bin domains to these variables
                    globalStartX = model.NewIntVar(0, (boundedM + 1) * binDx - item.Dx, f'xb1.{i}')
                    globalEndX = model.NewIntVar(item.Dx, (boundedM + 1) * binDx, f'xb2.{i}')
                    
                    startVariablesGlobalX.append(globalStartX)
                    endVariablesGlobalX.append(globalEndX)

                    yStart = model.NewIntVar(0, binDy - item.Dy,f'y1.{i}')
                    yEnd = model.NewIntVar(item.Dy, binDy,f'y2.{i}')

                    startVariablesY.append(yStart)
                    endVariablesY.append(yEnd)

        # interval variables
        intervalX = [model.NewIntervalVar(startVariablesGlobalX[i], items[i].Dx, endVariablesGlobalX[i],f'xival{i}') for i in range(n)]
        intervalY = [model.NewIntervalVar(startVariablesY[i], items[i].Dy, endVariablesY[i],f'yival{i}') for i in range(n)]

        # Also add lifted cuts from MIP via https://developers.google.com/optimization/reference/python/sat/python/cp_model#addforbiddenassignments
        #model.AddForbiddenAssignments([b[1], b[2]], [(0, 0), (1, 1)])
        self.AddIncompatibilityCuts(incompatibleItems, fixItemToBin, model, placedBinVariables)

        lowerBoundAreaBin = math.ceil(float(totalArea) / float(binDy * binDx))
        #lowerBound = min(lowerBoundAreaBin, lowerBoundBin)
        lowerBound = lowerBoundAreaBin

        # objective
        z = model.NewIntVar(lowerBound - 1, upperBoundBins - 1,'z')

        # constraints
        for i, item in enumerate(items):
            model.Add(startVariablesGlobalX[i] == startVariablesLocalX[i] + placedBinVariables[i]*binDx)
            #model.Add(xb2[i] == xb1[i] + item.Dx) # should already be ensured by interval variables

        model.AddNoOverlap2D(intervalX, intervalY)

        #upperBoundBin = m
        #demandsX = [item.Dx for item in items]
        demandsY = [item.Dy for item in items]
        model.AddCumulative(intervalX, demandsY, binDy)
        #model.AddCumulative(intervalY, demandsX, upperBoundBin * binDx)

        model.AddMaxEquality(z, [placedBinVariables[i] for i in range(n)])

        # objective
        model.Minimize(z + 1)    

        # solve model
        solver = cp_model.CpSolver()
        
        solver.parameters.log_search_progress = enableLogging
        solver.parameters.max_time_in_seconds = timeLimit
        solver.parameters.num_search_workers = 8

        #with open(f"Model_{0}.txt","a") as f:
        #    f.write(str(model.Proto()))

        rc = solver.Solve(model)
        #print(f"return code:{rc}")
        #print(f"status:{solver.StatusName()}")
        #print(f"Objective:{solver.ObjectiveValue()}")
        #print(f"Objective:{solver.BestObjectiveBound()}")
        #print(f"Objective:{solver.ResponseProto()}")
        #print(f"Objective:{solver.ResponseStats()}")
        
        status = solver.StatusName()
        if status == 'UNKNOWN':
            raise ValueError("Start solution could not be determined (CP status == UNKNOWN)")

        self.IsOptimal = 1 if solver.StatusName() == 'OPTIMAL' else 0
        self.LB = solver.BestObjectiveBound()
        self.UB = solver.ObjectiveValue()

        self.ItemBinAssignments = [solver.Value(placedBinVariables[i]) for i in range(n)]

        xArray = [solver.Value(startVariablesGlobalX[i]) for i in range(n)]
        yArray = [solver.Value(startVariablesY[i]) for i in range(n)]

        rectangles = ExtractDataForPlot(xArray, yArray, items, binDx, binDy)

        return rectangles

def main():
    #h, w, H, W, m = ReadExampleData()
    items, H, W = ReadBenchmarkData(9)

    solver = BinPackingSolverCP(items, H, W, 1, len(items), 30)
    rectangles = solver.Solve('OneBigBin')
    
    #solver = OneBigBinModel()
    #rectangles = solver.SolveOneBigBinModel(items, H, W, 1, len(items), 30)

    objBoundUB = solver.UB

    PlotSolution(objBoundUB * W, H, rectangles)

if __name__ == "__main__":
    main()
