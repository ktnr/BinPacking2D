from ortools.sat.python import cp_model

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
    
    """ from https://yetanothermathprogrammingconsultant.blogspot.com/2021/02/2d-bin-packing-with-google-or-tools-cp.html 
    https://yetanothermathprogrammingconsultant.blogspot.com/2021/02/2d-bin-packing.html """
    def BinPackingErwin(self, h, w, H, W, m, timeLimit = 3600, enableLogging = True):

        n = len(h)

        model = cp_model.CpModel()

        #
        # variables
        #

        # x and y
        x = [model.NewIntVar(0,W-w[i],f'x{i}') for i in range(n)]

        xb1 = [model.NewIntVar(0,m*W-w[i],f'xb1.{i}') for i in range(n)]
        xb2 = [model.NewIntVar(w[i],m*W,f'xb2.{i}') for i in range(n)]

        y1 = [model.NewIntVar(0,H-h[i],f'y1.{i}') for i in range(n)]
        y2 = [model.NewIntVar(h[i],H,f'y2.{i}') for i in range(n)]

        # interval variables
        xival = [model.NewIntervalVar(xb1[i],w[i],xb2[i],f'xival{i}') for i in range(n)]
        yival = [model.NewIntervalVar(y1[i],h[i],y2[i],f'yival{i}') for i in range(n)]

        # bin numbers
        b = [model.NewIntVar(0,m-1,f'b{i}') for i in range(n)]

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