from ortools.sat.python import cp_model
from ortools.sat.python.cp_model import Domain

import math
import pandas
import numpy

from BinPackingData import *
from BinPacking import BinPacking2D, Bin

import matplotlib

class BinPackingSolverCP():
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

            normalPatternsX, normalPatternsY = BinPacking2D.DetermineNormalPatterns(itemSubset, bin.Dx - itemI.Dx, bin.Dy - itemI.Dy)
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

                fixedBinNormalPatternsX, fixedBinNormalPatternsY = BinPacking2D.DetermineNormalPatterns(itemSubset, bin.Dx - itemI.Dx, bin.Dy - itemI.Dy, i*bin.Dx)
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

                binSpecificNormalPatternsX, binSpecificNormalPatternsY = BinPacking2D.DetermineNormalPatterns(itemSubset, bin.Dx - itemI.Dx, bin.Dy - itemI.Dy, b*bin.Dx)
                itemBinNormalPatternsX[i].extend(binSpecificNormalPatternsX)

        return itemSpecificNormalPatternsX, itemSpecificNormalPatternsY, itemBinNormalPatternsX

    def AddIncompatibilityCuts(self, incompatibleItems, fixItemToBin, model, binVariables):
        if incompatibleItems == None:
            return

        for i, j in incompatibleItems:
            if fixItemToBin[i] and fixItemToBin[j]:
                continue
            
            model.Add(binVariables[i] != binVariables[j])

    """ from https://yetanothermathprogrammingconsultant.blogspot.com/2021/02/2d-bin-packing-with-google-or-tools-cp.html 
    https://yetanothermathprogrammingconsultant.blogspot.com/2021/02/2d-bin-packing.html """
    def BinPackingErwin(self, items, h, w, H, W, lowerBoundBin, m, timeLimit = 3600, enableLogging = True, incompatibleItems = None):

        n = len(items)

        model = cp_model.CpModel()

        #preprocessing
        fixItemToBin = self.FixIncompatibleItems(incompatibleItems, n)

        binDomains = self.CreateReducedBinDomains(incompatibleItems, n, m, fixItemToBin)

        itemNormalPatternsX, itemNormalPatternsY, globalNormalPatternsX = [], [], []
        if self.EnableNormalPatterns:
            itemNormalPatternsX, itemNormalPatternsY, globalNormalPatternsX = self.CreateBinDependentNormalPatterns(incompatibleItems, fixItemToBin, items, m, W, H)


        # variables

        x = [] 
        b = [] # bin numbers
        xb1 = []
        xb2 = []
        y1 = []
        y2 = []
        totalArea = 0.0
        for i, item in enumerate(items):

            totalArea += item.Dx * item.Dy

            f = model.NewIntVarFromDomain(Domain.FromValues(binDomains[i]), f'b{i}')
            b.append(f)

            if self.EnableNormalPatterns:
                filteredStartLocalX = [p for p in itemNormalPatternsX[i] if (p % W) + item.Dx <= W]
                xStart = model.NewIntVarFromDomain(Domain.FromValues(filteredStartLocalX),f'x{i}')
                x.append(xStart)

                filteredStartX = [p for p in globalNormalPatternsX[i] if (p % W) + item.Dx <= W]
                filteredEndX = [p + item.Dx for p in filteredStartX]
                d = model.NewIntVarFromDomain(Domain.FromValues(filteredStartX), f'xb1.{i}')
                e = model.NewIntVarFromDomain(Domain.FromValues(filteredEndX), f'xb2.{i}')

                xb1.append(d)
                xb2.append(e)

                filteredStartY = [p for p in itemNormalPatternsY[i] if p + item.Dy <= H]
                filteredEndY = [p + item.Dy for p in filteredStartY]
                yStart = model.NewIntVarFromDomain(Domain.FromValues(filteredStartY),f'y1.{i}')
                yEnd = model.NewIntVarFromDomain(Domain.FromValues(filteredEndY),f'y2.{i}')

                y1.append(yStart)
                y2.append(yEnd)
            else:
                xStart = model.NewIntVar(0, W - item.Dx, f'x{i}')
                x.append(xStart)

                yStart = model.NewIntVar(0, H - item.Dy,f'y1.{i}')
                yEnd = model.NewIntVar(item.Dy, H,f'y2.{i}')

                y1.append(yStart)
                y2.append(yEnd)

                if fixItemToBin[i]:
                    d = model.NewIntVar(i*W, (i + 1)*W - item.Dx, f'xb1.{i}')
                    e = model.NewIntVar(i*W + item.Dx, (i + 1)*W, f'xb2.{i}')

                    xb1.append(d)
                    xb2.append(e)
                else:
                    boundedM = i if i < m else m - 1

                    # TODO: apply bin domains to these variables
                    d = model.NewIntVar(0, (boundedM + 1) * W - item.Dx, f'xb1.{i}')
                    e = model.NewIntVar(item.Dx, (boundedM + 1) * W, f'xb2.{i}')
                    
                    xb1.append(d)
                    xb2.append(e)


        # interval variables
        xival = [model.NewIntervalVar(xb1[i], items[i].Dx, xb2[i],f'xival{i}') for i in range(n)]
        yival = [model.NewIntervalVar(y1[i], items[i].Dy, y2[i],f'yival{i}') for i in range(n)]

        # Also add lifted cuts from MIP via https://developers.google.com/optimization/reference/python/sat/python/cp_model#addforbiddenassignments
        #model.AddForbiddenAssignments([b[1], b[2]], [(0, 0), (1, 1)])
        self.AddIncompatibilityCuts(incompatibleItems, fixItemToBin, model, b)

        lowerBoundAreaBin = math.ceil(float(totalArea) / float(H * W))
        lowerBound = min(lowerBoundAreaBin, lowerBoundBin)

        # objective
        z = model.NewIntVar(lowerBound - 1, m - 1,'z')

        # constraints
        for i, item in enumerate(items):
            model.Add(xb1[i] == x[i] + b[i]*W)
            model.Add(xb2[i] == xb1[i] + item.Dx) # should already be ensured by interval vars xival

        model.AddNoOverlap2D(xival, yival)

        model.AddMaxEquality(z, [b[i] for i in range(n)])

        # objective
        model.Minimize(z + 1)    

        # solve model
        solver = cp_model.CpSolver()
        # log does not work inside a Jupyter notebook
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

        self.IsOptimal = 1 if solver.StatusName() == 'OPTIMAL' else 0
        self.LB = solver.BestObjectiveBound()
        self.UB = solver.ObjectiveValue()

        self.ItemBinAssignments = [solver.Value(b[i]) for i in range(n)]

        xArray = [solver.Value(xb1[i]) for i in range(n)]
        yArray = [solver.Value(y1[i]) for i in range(n)]

        rectangles = ExtractDataForPlot(xArray, yArray, w, h, W, H)

        return rectangles

def main():
    #h, w, H, W, m = ReadExampleData()
    h, w, H, W, m = ReadBenchmarkData(10)

    solver = BinPackingSolverCP()
    rectangles = solver.BinPackingErwin(h, w, H, W, m)

    objBoundUB = solver.UB

    PlotSolution(objBoundUB * W, H, rectangles)

if __name__ == "__main__":
    main()
