from ortools.sat.python import cp_model
from ortools.sat.python.cp_model import Domain

import math
import pandas
import numpy

from BinPackingData import *

import matplotlib

class BinPackingSolverCP():
    def __init__(self):
        self.IsOptimal = False
        self.LB = -1
        self.UB = -1

        self.ItemBinAssignments = []
    
    def FixIncompatibleItems(self, incompatibleItems, numberOfItems):
        # This is a similar logic as in section 6.2 in Cote, Haouari, Iori (2019). 
        fixItemToBin = [False] * numberOfItems
        fixItemToBin[0] = True

        if incompatibleItems == None or len(incompatibleItems) == 0:
            return 1, fixItemToBin

        for i in range(1, numberOfItems):
            isIncompatible = True
            for j in range(0, i):
                if frozenset((i, j)) not in incompatibleItems:
                    isIncompatible = False
                    return i, fixItemToBin
            
            if isIncompatible:
                fixItemToBin[i] = True

        return numberOfItems - 1, fixItemToBin

    """ For reference, see section 6.1 from Cote, Iori (2018). """
    def CreateReducedBinDomains(self, incompatibleItems, numberOfItems, numberOfBins, fixItemToBin, lastFixedItemIndex):
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

    def AddIncompatibilityCuts(self, incompatibleItems, lastFixedItemIndex, model, binVariables):
        if incompatibleItems == None:
            return

        for i, j in incompatibleItems:
            if i > lastFixedItemIndex and j < lastFixedItemIndex:
                #binVariables[i].RemoveValue(j)
                pass

            if i < lastFixedItemIndex and j < lastFixedItemIndex:
                continue
            
            model.Add(binVariables[i] != binVariables[j])

    """ from https://yetanothermathprogrammingconsultant.blogspot.com/2021/02/2d-bin-packing-with-google-or-tools-cp.html 
    https://yetanothermathprogrammingconsultant.blogspot.com/2021/02/2d-bin-packing.html """
    def BinPackingErwin(self, h, w, H, W, m, timeLimit = 3600, enableLogging = True, incompatibleItems = None):

        n = len(h)

        model = cp_model.CpModel()

        #
        # variables
        #

        lastFixedIndex, fixItemToBin = self.FixIncompatibleItems(incompatibleItems, n)

        binDomains = self.CreateReducedBinDomains(incompatibleItems, n, m, fixItemToBin, lastFixedIndex)

        # x and y
        x = [model.NewIntVar(0, W-w[i],f'x{i}') for i in range(n)]

        # TODO: fix item x ranges to particular bin as in MIP
        b = [] # bin numbers
        xb1 = []
        xb2 = []
        for i in range(n):
            f = model.NewIntVarFromDomain(Domain.FromValues(binDomains[i]), f'b{i}')
            b.append(f)

            """
            d = model.NewIntVar(f*W, (f + 1)*W-w[i], f'xb1.{i}')
            e = model.NewIntVar(f*W + w[i], (f + 1)*W, f'xb2.{i}')

            xb1.append(d)
            xb2.append(e)
            """

            if fixItemToBin[i]:
                d = model.NewIntVar(i*W, (i + 1)*W-w[i], f'xb1.{i}')
                e = model.NewIntVar(i*W + w[i], (i + 1)*W, f'xb2.{i}')

                xb1.append(d)
                xb2.append(e)
            else:
                boundedM = i if i < m else m - 1

                # TODO: apply bin domains to these variables
                d = model.NewIntVar(0, boundedM*W-w[i], f'xb1.{i}')
                e = model.NewIntVar(w[i], boundedM*W, f'xb2.{i}')
                
                xb1.append(d)
                xb2.append(e)

        y1 = [model.NewIntVar(0,H-h[i],f'y1.{i}') for i in range(n)]
        y2 = [model.NewIntVar(h[i],H,f'y2.{i}') for i in range(n)]

        # interval variables
        xival = [model.NewIntervalVar(xb1[i],w[i],xb2[i],f'xival{i}') for i in range(n)]
        yival = [model.NewIntervalVar(y1[i],h[i],y2[i],f'yival{i}') for i in range(n)]

        self.AddIncompatibilityCuts(incompatibleItems, lastFixedIndex, model, b)

        totalArea = numpy.dot(w,h)
        lowerBoundBin = math.ceil(float(totalArea) / float(H * W))

        # objective
        z = model.NewIntVar(lowerBoundBin - 1, m - 1,'z')

        #
        # constraints
        #
        for i in range(n):
            model.Add(xb1[i] == x[i] + b[i]*W)
            model.Add(xb2[i] == xb1[i] + w[i])

        model.AddNoOverlap2D(xival,yival)

        model.AddMaxEquality(z,[b[i] for i in range(n)])

        # objective
        model.Minimize(z + 1)    

        #
        # solve model
        #
        solver = cp_model.CpSolver()
        # log does not work inside a Jupyter notebook
        solver.parameters.log_search_progress = enableLogging
        solver.parameters.max_time_in_seconds = timeLimit
        solver.parameters.num_search_workers = 8
        rc = solver.Solve(model)
        #print(f"return code:{rc}")
        #print(f"status:{solver.StatusName()}")
        #print(f"Objective:{solver.ObjectiveValue()}")
        #print(f"Objective:{solver.BestObjectiveBound()}")
        #print(f"Objective:{solver.ResponseProto()}")
        #print(f"Objective:{solver.ResponseStats()}")

        """
        if rc == 4:
            df = pandas.DataFrame({ 
                'bin' : [solver.Value(b[i]) for i in range(n)],
                'x'   : [solver.Value(x[i]) for i in range(n)],
                'y'   : [solver.Value(y1[i]) for i in range(n)],
                'w'   : w,
                'h'   : h})
            #display(df)
        """
        
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