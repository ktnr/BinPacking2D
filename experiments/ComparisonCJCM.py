import re
import os, sys
from os import path

# https://towardsdatascience.com/understanding-python-imports-init-py-and-pythonpath-once-and-for-all-4c5249ab6355
fpathParent = os.path.join(os.path.dirname(__file__), '..')
fpathSolver = os.path.join(os.path.dirname(__file__), '../packingSolver')
sys.path.append(fpathParent)
sys.path.append(fpathSolver)
#print(sys.path)

#from packingSolver import *
from packingSolver.OrthogonalPacking import OrthogonalPackingSolver
from packingSolver.PlacementPoints import PlacementPointStrategy
from packingSolver.BinPackingData import *

instanceFilter = [r'E02F17', r'E04F15', r'E00N23']
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

        print(f'{fileName} ', end='', flush=True)
        
        solver = OrthogonalPackingSolver(items, bin, PlacementPointStrategy.NormalPatterns)
        isFeasible = solver.Solve()

        if isFeasible:
            rectangles = ExtractDataForPlot(solver.PositionsX, solver.PositionsY, items, W, H)
            PlotSolution(W, H, rectangles)
            print(f'is feasible (#items = {len(items)})')
        else:
            print(f'is infeasible (#items = {len(items)})')