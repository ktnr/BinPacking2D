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
        elif modelType == 'StripPackOneBigBin':
            model = StripPackOneBigBinModel()
        else:
            raise ValueError("Invalid bin packing model type.")
        
        rectangles = model.Solve(self.items, self.binDy, self.binDx, self.lowerBoundBins, self.upperBoundBins, self.timeLimit, self.enableLogging, self.incompatibleItems)
       
        self.IsOptimal = model.IsOptimal
        self.LB = model.LB
        self.UB = model.UB
        self.ItemBinAssignments = model.ItemBinAssignments

        return rectangles

class PairwiseAssignmentModel:
    def __init__(self):
        self.IsOptimal = False
        self.LB = -1
        self.UB = -1

        self.EnableNormalPatterns = False

        self.ItemBinAssignments = []

    def Solve(self, items, binDy, binDx, lowerBoundBins, upperBoundBins, timeLimit = 3600, enableLogging = True, incompatibleItems = None):
        pass


class StripPackOneBigBinModel:
    def __init__(self):
        self.IsOptimal = False
        self.LB = -1
        self.UB = -1

        self.EnableNormalPatterns = False

        self.ItemBinAssignments = []

    class EndPositionPlacementAbortCallback(cp_model.CpSolverSolutionCallback):
        def __init__(self, binDx):
            cp_model.CpSolverSolutionCallback.__init__(self)
            self.binDx = binDx
            self.LB = 1
            self.UB = binDx * 100

        def on_solution_callback(self):
            """Called at each new solution."""
            lb = self.BestObjectiveBound()
            ub = self.ObjectiveValue()

            self.LB = math.ceil(float(lb) / float(self.binDx))
            self.UB = math.ceil(float(ub) / float(self.binDx))
            if int(self.UB) - int(self.LB) == 0:
                self.StopSearch()

    def Solve(self, items, binDy, binDx, lowerBoundBins, upperBoundBins, timeLimit = 3600, enableLogging = True, incompatibleItems = None):

        n = len(items)

        model = cp_model.CpModel()

        #preprocessing
        fixItemToBin = SymmetryBreaking.FixIncompatibleItems(incompatibleItems, n)

        #binDomains = SymmetryBreaking.CreateReducedBinDomains(incompatibleItems, n, upperBoundBins, fixItemToBin)

        itemNormalPatternsX, itemNormalPatternsY, globalNormalPatternsX = [], [], []
        if self.EnableNormalPatterns:
            itemNormalPatternsX, itemNormalPatternsY, globalNormalPatternsX = self.CreateBinDependentNormalPatterns(incompatibleItems, fixItemToBin, items, upperBoundBins, binDx, binDy)
        
        # variables

        startVariablesGlobalX = []
        #endVariablesGlobalX = []
        startVariablesY = []
        #endVariablesY = []
        totalArea = 0.0
        for i, item in enumerate(items):

            totalArea += item.Dx * item.Dy
            
            if self.EnableNormalPatterns:
                filteredStartX = [p for p in globalNormalPatternsX[i] if (p % binDx) + item.Dx <= binDx]
                filteredEndX = [p + item.Dx for p in filteredStartX]
                globalStartX = model.NewIntVarFromDomain(Domain.FromValues(filteredStartX), f'xb1.{i}')
                #globalEndX = model.NewIntVarFromDomain(Domain.FromValues(filteredEndX), f'xb2.{i}')

                startVariablesGlobalX.append(globalStartX)
                #endVariablesGlobalX.append(globalEndX)

                filteredStartY = [p for p in itemNormalPatternsY[i] if p + item.Dy <= binDy]
                #filteredEndY = [p + item.Dy for p in filteredStartY]
                yStart = model.NewIntVarFromDomain(Domain.FromValues(filteredStartY),f'y1.{i}')
                #yEnd = model.NewIntVarFromDomain(Domain.FromValues(filteredEndY),f'y2.{i}')

                startVariablesY.append(yStart)
                #endVariablesY.append(yEnd)
            else:
                if fixItemToBin[i]:
                    # do domain reduction
                    reducedDomainX = SymmetryBreaking.ReducedDomainX(binDx, item)
                    globalStartX = model.NewIntVar(i*binDx, i*binDx + reducedDomainX, f'xb1.{i}')
                    #globalEndX = model.NewIntVar(i*binDx + item.Dx, i*binDx + reducedDomainX + item.Dx, f'xb2.{i}')
                    #d = model.NewIntVar(i*W, (i + 1)*W - item.Dx, f'xb1.{i}')
                    #e = model.NewIntVar(i*W + item.Dx, (i + 1)*W, f'xb2.{i}')

                    startVariablesGlobalX.append(globalStartX)
                    #endVariablesGlobalX.append(globalEndX)
                    
                    reducedDomainY = SymmetryBreaking.ReducedDomainY(binDy, item)
                    yStart = model.NewIntVar(0, reducedDomainY,f'y1.{i}')
                    #yEnd = model.NewIntVar(item.Dy, reducedDomainY + item.Dy,f'y2.{i}')

                    startVariablesY.append(yStart)
                    #endVariablesY.append(yEnd)
                else:
                    # TODO: domain reduction for each bin where item i is the biggest placeable item
                    boundedM = i if i < upperBoundBins else upperBoundBins - 1

                    # TODO: apply bin domains to these variables
                    globalStartX = model.NewIntVarFromDomain(cp_model.Domain.FromIntervals([[binDx * m, binDx * (m + 1) - item.Dx] for m in range(boundedM + 1)]), f'xb1.{i}')
                    #globalEndX = model.NewIntVarFromDomain(cp_model.Domain.FromIntervals([[(binDx * m) + item.Dx, binDx * (m + 1)] for m in range(boundedM + 1)]), f'xb1.{i}')
                    
                    startVariablesGlobalX.append(globalStartX)
                    #endVariablesGlobalX.append(globalEndX)

                    yStart = model.NewIntVar(0, binDy - item.Dy,f'y1.{i}')
                    yEnd = model.NewIntVar(item.Dy, binDy,f'y2.{i}')

                    startVariablesY.append(yStart)
                    #endVariablesY.append(yEnd)

        # interval variables
        intervalX = [model.NewFixedSizeIntervalVar(startVariablesGlobalX[i], items[i].Dx,f'xival{i}') for i in range(n)]
        intervalY = [model.NewFixedSizeIntervalVar(startVariablesY[i], items[i].Dy,f'yival{i}') for i in range(n)]

        model.AddNoOverlap2D(intervalX, intervalY)

        #upperBoundBin = m
        #demandsX = [item.Dx for item in items]
        demandsY = [item.Dy for item in items]
        model.AddCumulative(intervalX, demandsY, binDy)
        #model.AddCumulative(intervalY, demandsX, upperBoundBin * binDx)

        # objective
        lowerBoundAreaBin = math.ceil(float(totalArea) / float(binDy * binDx))
        z = model.NewIntVar(lowerBoundAreaBin * binDx, (upperBoundBins + 1) * binDx,'z')
        model.AddMaxEquality(z, [startVariablesGlobalX[i] + item.Dx for i, item in enumerate(items)])

        # objective
        model.Minimize(z)    

        # solve model
        solver = cp_model.CpSolver()
        
        solver.parameters.log_search_progress = enableLogging
        solver.parameters.max_time_in_seconds = timeLimit
        solver.parameters.num_search_workers = 8
        solver.parameters.use_timetable_edge_finding_in_cumulative_constraint = True
        #solver.parameters.use_cumulative_in_no_overlap_2d = True
        #solver.parameters.use_precedences_in_disjunctive_constraint = True
        #solver.parameters.use_disjunctive_constraint_in_cumulative_constraint = True
        #solver.parameters.use_overload_checker_in_cumulative_constraint = True
        #solver.parameters.optimize_with_lb_tree_search = True
        #solver.parameters.binary_search_num_conflicts  = 8
        #solver.parameters.stop_after_first_solution = True
        #solver.parameters.symmetry_level = 2
        #solver.parameters.cp_model_presolve = False

        #with open(f"Model_{0}.txt","a") as f:
        #    f.write(str(model.Proto()))

        abortCallback = self.EndPositionPlacementAbortCallback(binDx)
        rc = solver.Solve(model, abortCallback)
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
        self.LB = math.ceil(float(solver.BestObjectiveBound()) / float(binDx))
        self.UB = math.ceil(float(solver.ObjectiveValue()) / float(binDx))

        #self.ItemBinAssignments = [solver.Value(placedBinVariables[i]) for i in range(n)]

        xArray = [solver.Value(startVariablesGlobalX[i]) for i in range(n)]
        yArray = [solver.Value(startVariablesY[i]) for i in range(n)]

        rectangles = ExtractDataForPlot(xArray, yArray, items, binDx, binDy)

        return rectangles


class OneBigBinModel:
    def __init__(self):
        self.IsOptimal = False
        self.LB = -1
        self.UB = -1

        self.EnableNormalPatterns = False

        self.ItemBinAssignments = []

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
        fixItemToBin = SymmetryBreaking.FixIncompatibleItems(incompatibleItems, n)

        binDomains = SymmetryBreaking.CreateReducedBinDomains(incompatibleItems, n, upperBoundBins, fixItemToBin)

        itemNormalPatternsX, itemNormalPatternsY, globalNormalPatternsX = [], [], []
        if self.EnableNormalPatterns:
            itemNormalPatternsX, itemNormalPatternsY, globalNormalPatternsX = SymmetryBreaking.CreateBinDependentNormalPatterns(incompatibleItems, fixItemToBin, items, upperBoundBins, binDx, binDy)
        
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
                    globalStartX = model.NewIntVarFromDomain(cp_model.Domain.FromIntervals([[binDx * m, binDx * (m + 1) - item.Dx] for m in range(boundedM + 1)]), f'xb1.{i}')
                    globalEndX = model.NewIntVarFromDomain(cp_model.Domain.FromIntervals([[(binDx * m) + item.Dx, binDx * (m + 1)] for m in range(boundedM + 1)]), f'xb1.{i}')
                    #globalStartX = model.NewIntVar(0, (boundedM + 1) * binDx - item.Dx, f'xb1.{i}')
                    #globalEndX = model.NewIntVar(item.Dx, (boundedM + 1) * binDx, f'xb2.{i}')
                    
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
        #solver.parameters.use_cumulative_in_no_overlap_2d = True
        #solver.parameters.use_disjunctive_constraint_in_cumulative_constraint = True
        #solver.parameters.use_timetable_edge_finding_in_cumulative_constraint = True
        #solver.parameters.use_overload_checker_in_cumulative_constraint = True
        #solver.parameters.max_time_in_seconds = 300
        #solver.parameters.optimize_with_lb_tree_search = True
        #solver.parameters.binary_search_num_conflicts  = 99
        #solver.parameters.stop_after_first_solution = True
        #solver.parameters.symmetry_level = 1
        #solver.parameters.cp_model_presolve = False

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
    items, H, W = ReadBenchmarkData(173)

    solver = BinPackingSolverCP(items, H, W, 1, len(items), 16*3600)
    rectangles = solver.Solve('StripPackOneBigBin')
    
    #solver = OneBigBinModel()
    #rectangles = solver.SolveOneBigBinModel(items, H, W, 1, len(items), 30)

    objBoundUB = solver.UB

    PlotSolution(objBoundUB * W, H, rectangles)

if __name__ == "__main__":
    main()
