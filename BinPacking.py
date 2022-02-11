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
    def __init__(
        self, 
        items, 
        binDy, 
        binDx, 
        lowerBoundBins, 
        upperBoundBins, 
        placementPointStrategy = PlacementPointStrategy.UnitDiscretization, 
        timeLimit = 3600, 
        enableLogging = True,
        preprocess = None):

        self.items = items
        self.binDx = binDx
        self.binDy = binDy
        self.lowerBoundBins = lowerBoundBins
        self.upperBoundBins = upperBoundBins
        self.placementPointStrategy = placementPointStrategy
        self.timeLimit = timeLimit
        self.enableLogging = enableLogging

        self.preprocess = preprocess
        
        self.IsOptimal = False
        self.LB = -1
        self.UB = -1

        self.ItemBinAssignments = []
    
    def Solve(self, modelType = 'OneBigBin'):
        bin = Bin(self.binDx, self.binDy)

        if self.preprocess == None:
            self.preprocess = Preprocess(self.items, bin, self.placementPointStrategy)
            self.preprocess.Run()

        if modelType == 'OneBigBin':
            model = OneBigBinModel(self.items, bin, self.preprocess, self.placementPointStrategy)
        elif modelType == 'StripPackOneBigBin':
            model = StripPackOneBigBinModel(self.items, bin, self.preprocess, self.placementPointStrategy)
        else:
            raise ValueError("Invalid bin packing model type.")
        
        rectangles = model.Solve(self.items, bin, self.lowerBoundBins, self.upperBoundBins, self.timeLimit, self.enableLogging)
       
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

        self.ItemBinAssignments = []

    def Solve(self, items, bin, lowerBoundBins, upperBoundBins, timeLimit = 3600, enableLogging = True):
        pass

class StripPackOneBigBinModel:
    def __init__(self, items, bin, preprocess, placementPointStrategy = PlacementPointStrategy.UnitDiscretization):
        self.IsOptimal = False
        self.LB = -1
        self.UB = -1

        self.placementPointStrategy = placementPointStrategy

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

        if preprocess == None:
            self.preprocess = Preprocess(items, bin, placementPointStrategy)
        else:
            self.preprocess = preprocess

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
            
            filteredStartX = [p for p in self.preprocess.GlobalPlacementPatternsX[i] if (p % binDx) + item.Dx <= binDx]
            globalStartX = self.model.NewIntVarFromDomain(Domain.FromValues(filteredStartX), f'xb1.{i}')

            self.startVariablesGlobalX.append(globalStartX)

            filteredStartY = [p for p in self.preprocess.ItemPlacementPatternsY[i] if p + item.Dy <= binDy]
            yStart = self.model.NewIntVarFromDomain(Domain.FromValues(filteredStartY),f'y1.{i}')

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

        self.loadingLength = self.model.NewIntVar((lowerBoundAreaBin - 1) * binDx + 1, (self.preprocess.UpperBoundsBin + 1) * binDx,'z')

        self.model.AddMaxEquality(self.loadingLength, [self.intervalX[i].EndExpr() for i, item in enumerate(items)])

        self.model.Minimize(self.loadingLength)    

    def SetParameters(self, solver, enableLogging, timeLimit):
        solver.parameters.log_search_progress = enableLogging
        solver.parameters.max_time_in_seconds = timeLimit
        solver.parameters.num_search_workers = 8
        #solver.parameters.use_timetable_edge_finding_in_cumulative_constraint = True
        #solver.parameters.use_cumulative_in_no_overlap_2d = True
        #solver.parameters.use_precedences_in_disjunctive_constraint = True
        #solver.parameters.use_disjunctive_constraint_in_cumulative_constraint = True
        #solver.parameters.use_overload_checker_in_cumulative_constraint = True
        #solver.parameters.optimize_with_lb_tree_search = True
        #solver.parameters.binary_search_num_conflicts  = 8
        #solver.parameters.stop_after_first_solution = True
        #solver.parameters.symmetry_level = 2
        #solver.parameters.cp_model_presolve = False

    def SolveModel(self):
        #with open(f"Model_{0}.txt","a") as f:
        #    f.write(str(model.Proto()))

        binDx = self.preprocess.Bin.Dx
        abortCallback = self.EndPositionPlacementAbortCallback(binDx)
        rc = self.solver.Solve(self.model, abortCallback)
        #print(f"return code:{rc}")
        #print(f"status:{solver.StatusName()}")
        #print(f"Objective:{solver.ObjectiveValue()}")
        #print(f"Objective:{solver.BestObjectiveBound()}")
        #print(f"Objective:{solver.ResponseProto()}")
        #print(f"Objective:{solver.ResponseStats()}")

    def ExtractSolution(self, items, bin):
        status = self.solver.StatusName()
        if status == 'UNKNOWN':
            raise ValueError("Start solution could not be determined (CP status == UNKNOWN)")

        self.IsOptimal = 1 if self.solver.StatusName() == 'OPTIMAL' else 0
        self.LB = math.ceil(float(self.solver.BestObjectiveBound()) / float(bin.Dx))
        self.UB = math.ceil(float(self.solver.ObjectiveValue()) / float(bin.Dx))

        self.ItemBinAssignments = []
        for i in range(len(items)):
            binId = math.floor(self.startVariablesGlobalX[i] / float(bin.Dx))
            self.ItemBinAssignments.append(binId)

        xArray = [self.solver.Value(self.startVariablesGlobalX[i]) for i in range(len(items))]
        yArray = [self.solver.Value(self.startVariablesY[i]) for i in range(len(items))]

        rectangles = ExtractDataForPlot(xArray, yArray, items, bin.Dx, bin.Dy)

        return rectangles

    def Solve(self, items, bin, lowerBoundBins, upperBoundBins, timeLimit = 3600, enableLogging = True):

        self.preprocess.Run()

        newItems = self.preprocess.ProcessedItems

        self.CreateVariables(newItems, bin.Dx, bin.Dy)
        self.CreateConstraints(newItems, bin.Dx, bin.Dy)
        self.CreateObjective(newItems, bin.Dx, bin.Dy)

        self.SetParameters(self.solver, enableLogging, timeLimit)
        self.SolveModel()

        rectangles = self.ExtractSolution(newItems, bin)

        return rectangles

class OneBigBinModel(StripPackOneBigBinModel):
    def __init__(self, items, bin, preprocess, placementPointStrategy = PlacementPointStrategy.UnitDiscretization):
        super().__init__(items, bin, preprocess, placementPointStrategy)
        
        self.startVariablesLocalX = []
        self.binCountVariables = []
        self.placedBinVariables = []
        self.ItemBinAssignments = []

    def CreateVariables(self, items, binDx, binDy):
        super().CreateVariables(items, binDx, binDy)

        self.CreateLocalStartVariablesX(items, binDx)
        self.CreateItemBinAssignmentVariables(items)
        self.CreateBinCountVariables(binDx, binDy)

    def CreateLocalStartVariablesX(self, items, binDx):
        itemNormalPatternsX = self.preprocess.ItemPlacementPatternsX
        for i, item in enumerate(items):
            filteredStartLocalX = [p for p in itemNormalPatternsX[i] if (p % binDx) + item.Dx <= binDx]
            xStart = self.model.NewIntVarFromDomain(Domain.FromValues(filteredStartLocalX),f'x{i}')
            self.startVariablesLocalX.append(xStart)

    def CreateBinCountVariables(self, binDx, binDy):
        lowerBoundAreaBin = math.ceil(float(self.itemArea) / float(binDx * binDy))
        #lowerBound = min(lowerBoundAreaBin, lowerBoundBin)
        lowerBound = lowerBoundAreaBin

        self.binCountVariables = self.model.NewIntVar(lowerBound - 1, self.preprocess.UpperBoundsBin - 1,'z')

    def CreateItemBinAssignmentVariables(self, items):
        binDomains = self.preprocess.BinDomains
        for i, item in enumerate(items):
            itemFeasibleBins = self.model.NewIntVarFromDomain(Domain.FromValues(binDomains[i]), f'b{i}')
            self.placedBinVariables.append(itemFeasibleBins)

    def CreateMaximumActiveBinConstraints(self, items):
        self.model.AddMaxEquality(self.binCountVariables, [self.placedBinVariables[i] for i in range(len(items))])

    def AddIncompatibilityCuts(self, incompatibleItems, fixItemToBin, model, binVariables):
        if incompatibleItems == None:
            return

        for i, j in incompatibleItems:
            if fixItemToBin[i] and fixItemToBin[j]:
                continue
            
            model.Add(binVariables[i] != binVariables[j])

    def CreateIntervalSynchronizationConstraints(self, items, binDx):
        for i, item in enumerate(items):
            self.model.Add(self.startVariablesGlobalX[i] == self.startVariablesLocalX[i] + self.placedBinVariables[i] * binDx)

    def CreateConstraints(self, items, binDx, binDy):
        super().CreateConstraints(items, binDx, binDy)

        self.CreateIntervalSynchronizationConstraints(items, binDx)
        self.CreateMaximumActiveBinConstraints(items)
        self.AddIncompatibilityCuts(self.preprocess.IncompatibleItems, self.preprocess.FixItemToBin, self.model, self.placedBinVariables)

    def CreateObjective(self, items, binDx, binDy):
        self.model.Minimize(self.binCountVariables + 1)

    def SolveModel(self):
        rc = self.solver.Solve(self.model)  

    def ExtractSolution(self, items, bin):
        status = self.solver.StatusName()
        if status == 'UNKNOWN':
            raise ValueError("Start solution could not be determined (CP status == UNKNOWN)")

        self.IsOptimal = 1 if self.solver.StatusName() == 'OPTIMAL' else 0
        self.LB = self.solver.BestObjectiveBound()
        self.UB = self.solver.ObjectiveValue()

        self.ItemBinAssignments = [self.solver.Value(self.placedBinVariables[i]) for i in range(len(items))]

        xArray = [self.solver.Value(self.startVariablesGlobalX[i]) for i in range(len(items))]
        yArray = [self.solver.Value(self.startVariablesY[i]) for i in range(len(items))]

        rectangles = ExtractDataForPlot(xArray, yArray, items, bin.Dx, bin.Dy)

        return rectangles

def main():
    items, H, W = ReadBenchmarkData(173)

    solver = BinPackingSolverCP(items, H, W, 1, len(items), PlacementPointStrategy.NormalPatterns, 16*3600)
    rectangles = solver.Solve('OneBigBin')

    objBoundUB = solver.UB

    PlotSolution(objBoundUB * W, H, rectangles)

if __name__ == "__main__":
    main()
