from ortools.sat.python import cp_model
from ortools.sat.python.cp_model import Domain

from Preprocess import Preprocess
from PlacementPoints import *
from SymmetryBreaking import *

from HelperIO import Converter

import math

class OrthogonalPackingSolver:
    def __init__(
        self, 
        items, 
        bin,
        placementPointStrategy = PlacementPointStrategy.UnitDiscretization, 
        timeLimit = 0, 
        enableLogging = False):

        self.items = items
        self.bin = bin
        self.placementPointStrategy = placementPointStrategy
        self.timeLimit = timeLimit
        self.enableLogging = enableLogging

        self.PositionsX = []
        self.PositionsY = []
    
    def Solve(self, instanceName = '9999', modelType = 'BaseModel'):
        if modelType == 'BaseModel':
            model = OrthogonalPacking2D(self.items, self.bin, self.placementPointStrategy)
        elif modelType == 'BranchAndCut':
            model = OrthogonalPackingRelaxed2D(self.items, self.bin, self.placementPointStrategy)
        else:
            raise ValueError("Invalid bin packing model type.")
        
        isFeasible = model.Solve(instanceName)
        
        numberOfItems = len(self.items)
        if numberOfItems == 1:
            self.PositionsX = [0.0]
            self.PositionsY = [0.0]
        else:
            self.PositionsX = model.StartPositionsX #[model.Solver.Value(model.StartX[i]) for i in range(numberOfItems)]
            self.PositionsY = model.StartPositionsY

        return isFeasible

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

class OrthogonalPackingBase2D:
    def __init__(self, items, bin, placementPointStrategy, positionsX = []):
        self.Model = cp_model.CpModel()
        self.Solver = None

        self.Items = items
        self.Bin = bin

        self.placementPointStrategy = placementPointStrategy
        self.fixedPositionsX = positionsX

        self.StartPositionsX = []
        self.StartPositionsY = []

        self.StartX = []
        self.EndX = []
        self.IntervalX = []

        self.StartY = []
        self.EndY = []
        self.IntervalY = []

    def Solve(self, instanceId):
        if len(self.Items) == 1:
            if self.Items[0].Dx <= self.Bin.Dx and self.Items[0].Dy <= self.Bin.Dy:
                return True
            else:
                return False

        # https://github.com/google/or-tools/blob/2cb85b4eead4c38e1c54b48044f92087cf165bce/ortools/sat/sat_parameters.proto
        self.Solver = cp_model.CpSolver()
        self.SetParameters()

        self.CreateVariables()
        self.CreateConstraints()
        self.FixVariablesX()

        """
        numberOfItems = len(self.Items)
        with open(f"Model_{numberOfItems}.txt","a") as f:
            f.write(str(self.Model.Proto()))
        """

        #status = self.Solver.Solve(self.Model)
        status = self.SolveModel()

        if status == cp_model.FEASIBLE or status == cp_model.OPTIMAL:
            self.ExtractSolution()
            return True
        elif status == cp_model.INFEASIBLE:
            return False
        elif status == cp_model.MODEL_INVALID:
            raise ValueError("Model invalid")
        elif status == cp_model.UNKNOWN:
            #Converter.ConvertFromSubproblem(self.Items, self.Bin, self.Solver.parameters.max_time_in_seconds, instanceId)
            #return False
            raise ValueError("Unkown status type")
        else:
            raise ValueError("Invalid status type")

        print(status)
        print(self.Solver.StatusName())

class OrthogonalPackingRelaxed2D(OrthogonalPackingBase2D):
    def __init__(self, items, bin, placementPointStrategy):
        super().__init__(items, bin, placementPointStrategy)

        self.AbortCallback = None

    class OrthogonalPackingFixedPositionCallback(cp_model.CpSolverSolutionCallback):
        def __init__(self, items, bin, placementPointStrategy, startVariablesX, startPositionsY, solver):
            cp_model.CpSolverSolutionCallback.__init__(self)
            self.items = items
            self.bin = bin
            self.placementPointStrategy = placementPointStrategy
            self.startVariablesX = startVariablesX
            self.StartPositionsY = startPositionsY
            self.solver = solver

        def on_solution_callback(self):
            if len(self.startVariablesX) != len(self.items):
                raise ValueError("Inconsistent variable count.")

            xPositions =  [self.Value(variable) for variable in self.startVariablesX]

            orthogonalPacking = OrthogonalPacking2D(self.items, self.bin, self.placementPointStrategy, xPositions)
            isFeasible = orthogonalPacking.Solve(0)

            if isFeasible:
                self.StartPositionsY = [orthogonalPacking.Solver.Value(variable) for variable in orthogonalPacking.StartY]
                self.StopSearch()
    
    def CreateVariables(self):
        binDx = self.Bin.Dx

        reducedItemIndex = SymmetryBreaking.DetermineMaximumItemIndexDx(self.Items)

        reducedItem = self.Items[reducedItemIndex]
        reducedDomainThresholdX = SymmetryBreaking.ReducedDomainX(binDx, reducedItem)

        for i, item in enumerate(self.Items):
        
            filteredItems = [itemJ for j, itemJ in enumerate(self.Items) if i != j]

            placementPointsX, placementPointsY = PlacementPointGenerator.CreatePlacementPatterns(self.placementPointStrategy, item, filteredItems, self.Bin)

            if i == reducedItemIndex:
                placementPointsStartX = [p for p in placementPointsX if p + item.Dx <= binDx and p <= reducedDomainThresholdX]
            else:
                placementPointsStartX = [p for p in placementPointsX if p + item.Dx <= binDx] # unnecessary, is satisfied by construction

            x1 = self.Model.NewIntVarFromDomain(Domain.FromValues(placementPointsStartX), f'x1.{i}')

            self.StartX.append(x1)

            intervalX = self.Model.NewFixedSizeIntervalVar(x1, item.Dx, f'xI{i}')

            self.IntervalX.append(intervalX)
    
    def CreateConstraints(self):
        self.CreateCumulativeConstraints()

    def CreateCumulativeConstraints(self):
        demandsY = [item.Dy for item in self.Items]

        self.Model.AddCumulative(self.IntervalX, demandsY, self.Bin.Dy)
    
    def ExtractSolution(self):
        self.StartPositionsX = [self.Solver.Value(variable) for variable in self.StartX]
        self.StartPositionsY = self.AbortCallback.StartPositionsY

    def FixVariablesX(self):
        return

    def SetParameters(self):
        self.Solver.parameters.log_search_progress = False 
        self.Solver.parameters.num_search_workers = 1
        self.Solver.parameters.enumerate_all_solutions = True
        #self.Solver.parameters.use_cumulative_in_no_overlap_2d = True
        #self.Solver.parameters.use_disjunctive_constraint_in_cumulative_constraint = True
        #self.Solver.parameters.use_timetable_edge_finding_in_cumulative_constraint = True
        #self.Solver.parameters.use_overload_checker_in_cumulative_constraint = True
        #self.Solver.parameters.max_time_in_seconds = 300
        #self.Solver.parameters.optimize_with_lb_tree_search = True;
        #self.Solver.parameters.binary_search_num_conflicts  = 99;
        #self.Solver.parameters.stop_after_first_solution = True;
        #self.Solver.parameters.symmetry_level = 1;
        #self.Solver.parameters.cp_model_presolve = False

    def SolveModel(self):
        self.AbortCallback = self.OrthogonalPackingFixedPositionCallback(self.Items, self.Bin, self.placementPointStrategy, self.StartX, self.StartPositionsY, self.Solver)
        status = self.Solver.Solve(self.Model, self.AbortCallback)
        #status = self.Solver.StatusName()

        return status

# https://www.xiang.dev/cp-sat/
class OrthogonalPacking2D(OrthogonalPackingBase2D):
    def __init__(self, items, bin, placementPointStrategy, positionsX = []):
        super().__init__(items, bin, placementPointStrategy, positionsX)

    def CreateVariables(self):
        binDx = self.Bin.Dx
        binDy = self.Bin.Dy

        reducedItemIndex = SymmetryBreaking.DetermineMaximumItemIndexDx(self.Items)

        reducedItem = self.Items[reducedItemIndex]
        reducedDomainThresholdX = SymmetryBreaking.ReducedDomainX(binDx, reducedItem)
        reducedDomainThresholdY = SymmetryBreaking.ReducedDomainY(binDy, reducedItem)

        # Consider modifying domains according to Cote, Iori (2018): The meet-in-the-middle principle for cutting and packing problems.
        for i, item in enumerate(self.Items):
        
            filteredItems = [itemJ for j, itemJ in enumerate(self.Items) if i != j]

            placementPointsX, placementPointsY = PlacementPointGenerator.CreatePlacementPatterns(self.placementPointStrategy, item, filteredItems, self.Bin)

            if i == reducedItemIndex:
                placementPointsStartX = [p for p in placementPointsX if p + item.Dx <= binDx and p <= reducedDomainThresholdX]
            else:
                placementPointsStartX = [p for p in placementPointsX if p + item.Dx <= binDx] # unnecessary, is satisfied by construction

            x1 = self.Model.NewIntVarFromDomain(Domain.FromValues(placementPointsStartX), f'x1.{i}')

            #placementPointsEndX = [p + item.Dx for p in placementPointsStartX]
            #x2 = self.Model.NewIntVarFromDomain(Domain.FromValues(placementPointsEndX), f'x2.{i}')

            self.StartX.append(x1)
            #self.EndX.append(x2)

            if i == reducedItemIndex:
                placementPointsStartY = [p for p in placementPointsY if p + item.Dy <= binDy and p <= reducedDomainThresholdY]
            else:
                placementPointsStartY = [p for p in placementPointsY if p + item.Dy <= binDy]  # is satisfied by construction

            y1 = self.Model.NewIntVarFromDomain(Domain.FromValues(placementPointsStartY), f'y1.{i}')
            #placementPointsEndY = [p + item.Dy for p in placementPointsStartY]
            #y2 = self.Model.NewIntVarFromDomain(Domain.FromValues(placementPointsEndY), f'y2.{i}')

            self.StartY.append(y1)
            #self.EndY.append(y2)

            intervalX = self.Model.NewFixedSizeIntervalVar(x1, item.Dx, f'xI{i}')
            intervalY = self.Model.NewFixedSizeIntervalVar(y1, item.Dy, f'yI{i}')

            self.IntervalX.append(intervalX)
            self.IntervalY.append(intervalY)

    def FixVariablesX(self):
        for i, placement in enumerate(self.fixedPositionsX):
            #self.StartX[i].FixVariable(placement) # will be available in future releases
            self.Model.Add(self.StartX[i] == placement)

    def CreateConstraints(self):
        self.CreateNoOverlapConstraints()
        self.CreateCumulativeConstraints()

    def CreateNoOverlapConstraints(self):
        self.Model.AddNoOverlap2D(self.IntervalX, self.IntervalY)

    def CreateCumulativeConstraints(self):
        demandsX = [item.Dx for item in self.Items]
        demandsY = [item.Dy for item in self.Items]

        self.Model.AddCumulative(self.IntervalY, demandsX, self.Bin.Dx)
        self.Model.AddCumulative(self.IntervalX, demandsY, self.Bin.Dy)

    def ExtractSolution(self):
        self.StartPositionsX = [self.Solver.Value(variable) for variable in self.StartX]
        self.StartPositionsY = [self.Solver.Value(variable) for variable in self.StartY]

    def SetParameters(self):
        self.Solver.parameters.log_search_progress = False 
        self.Solver.parameters.num_search_workers = 8
        #self.Solver.parameters.use_cumulative_in_no_overlap_2d = True
        #self.Solver.parameters.use_disjunctive_constraint_in_cumulative_constraint = True
        #self.Solver.parameters.use_timetable_edge_finding_in_cumulative_constraint = True
        #self.Solver.parameters.use_overload_checker_in_cumulative_constraint = True
        #self.Solver.parameters.max_time_in_seconds = 300
        #self.Solver.parameters.optimize_with_lb_tree_search = True;
        #self.Solver.parameters.binary_search_num_conflicts  = 99;
        #self.Solver.parameters.stop_after_first_solution = True;
        #self.Solver.parameters.symmetry_level = 1;
        #self.Solver.parameters.cp_model_presolve = False

    def SolveModel(self):
        return self.Solver.Solve(self.Model)

