from ortools.sat.python import cp_model
from ortools.sat.python.cp_model import Domain

import sys
import math
import pandas
import numpy

from Preprocess import Preprocess
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
        bin = Bin(self.binDx, self.binDy)
        if modelType == 'OneBigBin':
            model = OneBigBinModel()
        elif modelType == 'StripPackOneBigBin':
            model = StripPackOneBigBinModel(self.items, bin)
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
    def __init__(self, items, bin):
        self.IsOptimal = False
        self.LB = -1
        self.UB = -1

        self.EnableNormalPatterns = False

        self.ItemBinAssignments = []

        self.solver = cp_model.CpSolver()
        self.model = cp_model.CpModel()

        self.startVariablesGlobalX = []
        self.startVariablesY = []
        self.intervalX = []
        self.intervalY = []
        self.loadingLength = None

        self.itemArea = 0.0
        self.lowerBoundAreaBin = 0.0
        #self.upperBoundBins = sys.maxsize

        self.preprocess = Preprocess(items, bin)
        #self.globalNormalPatternsX = []
        #self.itemNormalPatternsX = []
        #self.itemNormalPatternsY = []
        #self.fixItemToBin = []

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

    def CreateVariables(self, items, binDx, binDy):
        self.CreatePlacementVariables(items, binDx, binDy)
        self.CreateIntervalVariables(items)

    def CreatePlacementVariables(self, items, binDx, binDy):
        for i, item in enumerate(items):

            self.itemArea += item.Dx * item.Dy
            
            if self.EnableNormalPatterns:
                filteredStartX = [p for p in self.preprocess.GlobalNormalPatternsX[i] if (p % binDx) + item.Dx <= binDx]
                globalStartX = self.model.NewIntVarFromDomain(Domain.FromValues(filteredStartX), f'xb1.{i}')

                self.startVariablesGlobalX.append(globalStartX)

                filteredStartY = [p for p in self.preprocess.ItemNormalPatternsY[i] if p + item.Dy <= binDy]
                yStart = self.model.NewIntVarFromDomain(Domain.FromValues(filteredStartY),f'y1.{i}')

                self.startVariablesY.append(yStart)
            else:
                if self.preprocess.FixItemToBin[i]:
                    # do domain reduction
                    reducedDomainX = SymmetryBreaking.ReducedDomainX(binDx, item)
                    globalStartX = self.model.NewIntVar(i*binDx, i*binDx + reducedDomainX, f'xb1.{i}')

                    self.startVariablesGlobalX.append(globalStartX)
                    
                    reducedDomainY = SymmetryBreaking.ReducedDomainY(binDy, item)
                    yStart = self.model.NewIntVar(0, reducedDomainY,f'y1.{i}')

                    self.startVariablesY.append(yStart)
                else:
                    # TODO: domain reduction for each bin where item i is the biggest placeable item
                    boundedM = i if i < self.preprocess.UpperBoundsBin else self.preprocess.UpperBoundsBin - 1

                    # TODO: apply bin domains to these variables
                    globalStartX = self.model.NewIntVarFromDomain(cp_model.Domain.FromIntervals([[binDx * m, binDx * (m + 1) - item.Dx] for m in range(boundedM + 1)]), f'xb1.{i}')
                    
                    self.startVariablesGlobalX.append(globalStartX)

                    yStart = self.model.NewIntVar(0, binDy - item.Dy,f'y1.{i}')

                    self.startVariablesY.append(yStart)

    def CreateIntervalVariables(self, items):
        self.intervalX = [self.model.NewFixedSizeIntervalVar(self.startVariablesGlobalX[i], items[i].Dx,f'xival{i}') for i in range(len(items))]
        self.intervalY = [self.model.NewFixedSizeIntervalVar(self.startVariablesY[i], items[i].Dy,f'yival{i}') for i in range(len(items))]

    def CreateConstraints(self, items, binDx, binDy):
        self.model.AddNoOverlap2D(self.intervalX, self.intervalY)

        demandsY = [item.Dy for item in items]
        self.model.AddCumulative(self.intervalX, demandsY, binDy)
        #demandsX = [item.Dx for item in items]
        #model.AddCumulative(intervalY, demandsX, upperBoundBin * binDx)

    def CreateObjective(self, items, binDx, binDy):
        lowerBoundAreaBin = math.ceil(float(self.itemArea) / float(binDy * binDx))

        self.loadingLength = self.model.NewIntVar(lowerBoundAreaBin * binDx, (self.preprocess.UpperBoundsBin + 1) * binDx,'z')

        self.model.AddMaxEquality(self.loadingLength, [self.startVariablesGlobalX[i] + item.Dx for i, item in enumerate(items)])

        self.model.Minimize(self.loadingLength)    

    def SetParameters(self, solver, enableLogging, timeLimit):
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

    def Solve(self, items, binDy, binDx, lowerBoundBins, upperBoundBins, timeLimit = 3600, enableLogging = True, incompatibleItems = None):

        self.preprocess.Run()

        newItems = self.preprocess.ProcessedItems

        """
        self.upperBoundBins = upperBoundBins

        #preprocessing
        self.fixItemToBin = SymmetryBreaking.FixIncompatibleItems(incompatibleItems, n)

        #binDomains = SymmetryBreaking.CreateReducedBinDomains(incompatibleItems, n, upperBoundBins, fixItemToBin)

        if self.EnableNormalPatterns:
            self.itemNormalPatternsX, self.itemNormalPatternsY, self.globalNormalPatternsX = SymmetryBreaking.CreateBinDependentNormalPatterns(incompatibleItems, self.fixItemToBin, items, upperBoundBins, binDx, binDy)
        """

        self.CreateVariables(newItems, binDx, binDy)
        self.CreateConstraints(newItems, binDx, binDy)
        self.CreateObjective(newItems, binDx, binDy)

        self.SetParameters(self.solver, enableLogging, timeLimit)
        # solve model

        #with open(f"Model_{0}.txt","a") as f:
        #    f.write(str(model.Proto()))

        abortCallback = self.EndPositionPlacementAbortCallback(binDx)
        rc = self.solver.Solve(self.model, abortCallback)
        #print(f"return code:{rc}")
        #print(f"status:{solver.StatusName()}")
        #print(f"Objective:{solver.ObjectiveValue()}")
        #print(f"Objective:{solver.BestObjectiveBound()}")
        #print(f"Objective:{solver.ResponseProto()}")
        #print(f"Objective:{solver.ResponseStats()}")
        
        status = self.solver.StatusName()
        if status == 'UNKNOWN':
            raise ValueError("Start solution could not be determined (CP status == UNKNOWN)")

        self.IsOptimal = 1 if self.solver.StatusName() == 'OPTIMAL' else 0
        self.LB = math.ceil(float(self.solver.BestObjectiveBound()) / float(binDx))
        self.UB = math.ceil(float(self.solver.ObjectiveValue()) / float(binDx))

        #self.ItemBinAssignments = [solver.Value(placedBinVariables[i]) for i in range(n)]

        xArray = [self.solver.Value(self.startVariablesGlobalX[i]) for i in range(len(items))]
        yArray = [self.solver.Value(self.startVariablesY[i]) for i in range(len(items))]

        rectangles = ExtractDataForPlot(xArray, yArray, newItems, binDx, binDy)

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
        fixItemToBin = SymmetryBreaking.DetermineFixedItems(incompatibleItems, n)

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
