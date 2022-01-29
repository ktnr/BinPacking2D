from ortools.sat.python import cp_model
from ortools.sat.python.cp_model import Domain

from PlacementPoints import *

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
class OrthogonalPacking2D:
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
                placementPointsX, placementPointsY = PlacementPointGenerator.DetermineNormalPatterns(filteredItems, binDx - item.Dx, binDy - item.Dy)
            elif placementPointStrategy == PlacementPointStrategy.MeetInTheMiddlePatterns:
                placementPointsX, placementPointsY = PlacementPointGenerator.DetermineMeetInTheMiddlePatterns(filteredItems, item, binDx, binDy)
                raise ValueError("Meet-in-the-middle patterns might not be accurate, see #2.")
            elif placementPointStrategy == PlacementPointStrategy.MinimalMeetInTheMiddlePatterns:
                placementPointsX, placementPointsY = PlacementPointGenerator.DetermineMinimalMeetInTheMiddlePatterns(filteredItems, item, binDx, binDy)
                raise ValueError("Minimal meet-in-the-middle patterns might not be accurate, see #2.")
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