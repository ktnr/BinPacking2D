from ortools.sat.python import cp_model
from ortools.sat.python.cp_model import Domain

import sys
import math
import pandas
import numpy

from Preprocess import PreprocessBinPacking
from BinPackingData import *
from PlacementPoints import PlacementPointStrategy

from Model import *

import matplotlib

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
            self.preprocess = PreprocessBinPacking(self.items, bin, self.placementPointStrategy)
            self.preprocess.Run()

        if modelType == 'OneBigBin':
            model = OneBigBinModel(self.items, bin, self.preprocess, self.placementPointStrategy)
        elif modelType == 'StripPackOneBigBin':
            model = StripPackOneBigBinModel(self.items, bin, self.preprocess, self.placementPointStrategy)
        elif modelType == 'PairwiseAssignment':
            model = PairwiseAssignmentModel(self.items, bin, self.preprocess, self.placementPointStrategy)
        else:
            raise ValueError("Invalid bin packing model type.")
        
        rectangles = model.Solve(self.items, bin, self.lowerBoundBins, self.upperBoundBins, self.timeLimit, self.enableLogging)
       
        self.IsOptimal = model.IsOptimal
        self.LB = model.LB
        self.UB = model.UB
        self.ItemBinAssignments = model.ItemBinAssignments

        return rectangles

class BinPackBaseModel:
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

        self.itemArea = 0.0
        self.lowerBoundAreaBin = 0.0

        if preprocess == None:
            self.preprocess = PreprocessBinPacking(items, bin, placementPointStrategy)
        else:
            self.preprocess = preprocess

    def RetrieveLowerBound(self):
        return self.solver.BestObjectiveBound()

    def RetrieveUpperBound(self):
        return self.solver.ObjectiveValue()

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
    
    def CreateObjective(self, items, binDx, binDy):
        raise ValueError("CreateObjective() not implemented in BinPackBaseModel class.")

    def SolveModel(self):
        rc = self.solver.Solve(self.model)  

    def ExtractSolution(self, items, bin):
        status = self.solver.StatusName()
        if status == 'UNKNOWN':
            raise ValueError("Start solution could not be determined (CP status == UNKNOWN)")

        self.IsOptimal = 1 if self.solver.StatusName() == 'OPTIMAL' else 0
        self.LB = self.RetrieveLowerBound()
        self.UB = self.RetrieveUpperBound()

        self.DetermineItemBinAssignments(items, bin)

        xArray = [self.solver.Value(self.startVariablesGlobalX[i]) for i in range(len(items))]
        yArray = [self.solver.Value(self.startVariablesY[i]) for i in range(len(items))]

        rectangles = ExtractDataForPlot(xArray, yArray, items, bin.Dx, bin.Dy)

        return rectangles

    def DetermineItemBinAssignments(self, items, bin):
        self.ItemBinAssignments = []
        for i in range(len(items)):
            binId = math.floor(self.solver.Value(self.startVariablesGlobalX[i]) / float(bin.Dx))
            self.ItemBinAssignments.append(binId)

""" A simpler version of this model can be found at https://yetanothermathprogrammingconsultant.blogspot.com/2021/02/2d-bin-packing-with-google-or-tools-cp.html. """
class OneBigBinModel(BinPackBaseModel):
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

""" A minimal bin packing model without bin assignment variables that uses a strip packing search function. """
class StripPackOneBigBinModel(BinPackBaseModel):
    def __init__(self, items, bin, preprocess, placementPointStrategy = PlacementPointStrategy.UnitDiscretization):
        super().__init__(items, bin, preprocess, placementPointStrategy)

        self.loadingLength = None

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

    def CreateObjective(self, items, binDx, binDy):
        lowerBoundAreaBin = math.ceil(float(self.itemArea) / float(binDy * binDx))

        self.loadingLength = self.model.NewIntVar((lowerBoundAreaBin - 1) * binDx + 1, (self.preprocess.UpperBoundsBin + 1) * binDx,'z')

        self.model.AddMaxEquality(self.loadingLength, [self.intervalX[i].EndExpr() for i, item in enumerate(items)])

        self.model.Minimize(self.loadingLength)    

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

    def RetrieveLowerBound(self):
        return math.ceil(float(self.solver.BestObjectiveBound()) / float(self.preprocess.Bin.Dx))

    def RetrieveUpperBound(self):
        return math.ceil(float(self.solver.ObjectiveValue()) / float(self.preprocess.Bin.Dx))

""" Similar to the model of Aliaksei Vaidzelevich, cp. https://yetanothermathprogrammingconsultant.blogspot.com/2021/02/2d-bin-packing-with-google-or-tools-cp.html. 
    It is the same OneBigBin but with a different item bin assignment formulation. """
class PairwiseAssignmentModel(BinPackBaseModel):
    def __init__(self, items, bin, preprocess, placementPointStrategy = PlacementPointStrategy.UnitDiscretization):
        super().__init__(items, bin, preprocess, placementPointStrategy)

        self.itemBinAssignmentVariables = []
        self.binCountVariables = None

    def CreateVariables(self, items, binDx, binDy):
        self.CreatePlacementVariables(items, binDx, binDy)
        self.CreateIntervalVariables(items)
        self.CreateItemBinAssignmentVariables(items, self.preprocess.UpperBoundsBin)
        self.CreateBinCountVariables(binDx, binDy)

    def CreateItemBinAssignmentVariables(self, items, numberOfBins):
        self.itemBinAssignmentVariables = [[self.model.NewBoolVar(f'lit[{i}][{j}]') for j in range(numberOfBins)] for i in range(len(items))]

    def CreateBinCountVariables(self, binDx, binDy):
        lowerBoundAreaBin = math.ceil(float(self.itemArea) / float(binDx * binDy))
        #lowerBound = min(lowerBoundAreaBin, lowerBoundBin)
        lowerBound = lowerBoundAreaBin

        self.binCountVariables = self.model.NewIntVar(lowerBound, self.preprocess.UpperBoundsBin,'z')

    def CreateConstraints(self, items, binDx, binDy):
        self.model.AddNoOverlap2D(self.intervalX, self.intervalY)

        demandsY = [item.Dy for item in items]
        self.model.AddCumulative(self.intervalX, demandsY, binDy)

        self.CreateItemBinAssignmentConstraints(items, binDx, self.preprocess.UpperBoundsBin)

    def CreateItemBinAssignmentConstraints(self, items, binDx, numberOfBins):
        for i, item in enumerate(items):
            self.model.Add(sum(self.itemBinAssignmentVariables[i]) == 1)
            self.model.Add(sum(j * binDx * self.itemBinAssignmentVariables[i][j] for j in range(numberOfBins)) <= self.intervalX[i].StartExpr())
            self.model.Add(sum((j + 1) * binDx * self.itemBinAssignmentVariables[i][j] for j in range(numberOfBins)) >= self.intervalX[i].EndExpr())

    def CreateObjective(self, items, binDx, binDy):
        for i in range(len(items)):
            for j in range(self.preprocess.UpperBoundsBin):
                self.model.Add(self.binCountVariables >= j + 1).OnlyEnforceIf(self.itemBinAssignmentVariables[i][j])

        self.model.Minimize(self.binCountVariables)

def main():
    path = 'data/input/BPP/CLASS'
    instanceId = 9
    fileName = path 
    items, H, W = ReadBenchmarkData(path, str(instanceId) + '.json')

    solver = BinPackingSolverCP(items, H, W, 1, len(items), PlacementPointStrategy.NormalPatterns, 16*3600)
    rectangles = solver.Solve('OneBigBin')

    objBoundUB = solver.UB

    PlotSolution(objBoundUB * W, H, rectangles)

if __name__ == "__main__":
    main()
