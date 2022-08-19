import re
import os, sys
import time
from os import path

# https://towardsdatascience.com/understanding-python-imports-init-py-and-pythonpath-once-and-for-all-4c5249ab6355
fpathParent = os.path.join(os.path.dirname(__file__), '..')
fpathSolver = os.path.join(os.path.dirname(__file__), '../packingSolver')
sys.path.append(fpathParent)
sys.path.append(fpathSolver)
#print(sys.path)

#from packingSolver import *
from packingSolver.OrthogonalPacking import OrthogonalPackingSolver, ParametersPackingCP
from packingSolver.PlacementPoints import PlacementPointStrategy
from packingSolver.BinPackingData import *

def SetSearchParameters(parameters):
    SetDefault(parameters)
    #SetBasic(parameters)
    #SetFull(parameters)

def SetDefault(parameters):
    pass

def SetBasic(parameters):
    parameters.EnableDisjunctiveConstraintsInCumulative = True
    parameters.EnableTimeTableEdgeFindingInCumulative = True
    parameters.EnableTimeTablingInNoOverlap = True
    parameters.EnableEnergeticReasoningInNoOverlap = True

def SetFull(parameters):
    #parameters.UseCombinedNoOverlap = False
    parameters.EnableDisjunctiveConstraintsInCumulative = True
    parameters.EnableTimeTableEdgeFindingInCumulative = True
    parameters.EnableTimeTablingInNoOverlap = True
    parameters.EnableOverloadCheckingInCumulative = True
    parameters.EnableEnergeticReasoningInNoOverlap = True

runtimesMaxCJCM = {
    "E00N23": 70.0,
    "E00X23": 160.0
}

runtimesMinCJCM = {
    "E00N23": 1.14,
    "E00X23": 6.0
}

instanceFilter = [r'E00N23', r'E00X23']
combinedRegex = "(" + ")|(".join(instanceFilter) + ")"
path = 'data/input/OPP/CJCM08'
for root, dirs, files in os.walk(path):
    for fileName in files:
        if '.json' not in fileName:
            continue

        if len(instanceFilter) > 0 and not re.match(combinedRegex, fileName):
            continue
        
        items, H, W = ReadBenchmarkData(path, fileName)
        bin = Bin(W, H)

        instanceName = fileName.split('.')[0]
        #print(f'{fileName} ', end='', flush=True)
        
        t1 = time.time()

        parametersCP = ParametersPackingCP()
        SetSearchParameters(parametersCP)

        solver = OrthogonalPackingSolver(items, bin, PlacementPointStrategy.StandardUnitDiscretization)
        isFeasible = solver.Solve(True, instanceName, 'BaseModel', parametersCP)
            
        t2 = time.time()
        elapsedTime = t2 - t1

        feasibility = "feasible" if isFeasible else "infeasible"
        print(f'{instanceName} is {feasibility} in {elapsedTime:.2f}s vs. CJCM08 in max = {runtimesMaxCJCM[instanceName]}s and min = {runtimesMinCJCM[instanceName]}s')
        
        if isFeasible:
            #rectangles = ExtractDataForPlot(solver.PositionsX, solver.PositionsY, items, W, H)
            #PlotSolution(W, H, rectangles)
            pass