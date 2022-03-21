import re
import os, sys
import time

from os import path

import pandas as pd
import json

# https://towardsdatascience.com/understanding-python-imports-init-py-and-pythonpath-once-and-for-all-4c5249ab6355
fpathParent = os.path.join(os.path.dirname(__file__), '..')
fpathSolver = os.path.join(os.path.dirname(__file__), '../packingSolver')
sys.path.append(fpathParent)
sys.path.append(fpathSolver)
#print(sys.path)

#from packingSolver import *
from packingSolver.BinPacking import BinPackingSolverCP
from packingSolver.PlacementPoints import PlacementPointStrategy
from packingSolver.BinPackingData import *
from packingSolver.Preprocess import PreprocessBinPacking
from packingSolver.BranchAndCut import BinPackingMip

solutions = []

timeLimit = 180
path = 'data/input/BPP/CLASS'
for instance in range(1, 501):
#for instance in hardInstances:
    currentInstanceId = instance
    items, H, W = ReadBenchmarkData(path, str(instance) + '.json')
    bin = Bin(W, H)

    preprocess = PreprocessBinPacking(items, bin, PlacementPointStrategy.NormalPatterns)
    preprocess.Run()
    newItems = preprocess.ProcessedItems
    numberOfItems = len(newItems)

    lowerBoundBin = BinPackingMip.DetermineLowerBoundBinPacking(newItems, bin, numberOfItems, preprocess, True)
    
    for binPackingModel in ['OneBigBin', 'PairwiseAssignment', 'StripPackOneBigBin']:
        solverCP = BinPackingSolverCP(items, H, W, lowerBoundBin, len(items), preprocess.PlacementPointStrategy, timeLimit, False, preprocess)
        
        t1 = time.time()
        
        rectangles = solverCP.Solve('OneBigBin')

        t2 = time.time()
        elapsedTime = t2 - t1

        solution = {'Instance': instance, 'ModelType': binPackingModel, 'LB': solverCP.LB, 'UB': solverCP.UB, 'Runtime': elapsedTime}
        solutions.append(solution)

        print(f'Instance {instance}: Bounds = [{solverCP.LB}, {solverCP.UB}] in {elapsedTime:.2f}s with {binPackingModel}')

df = pd.DataFrame(solutions)

result = df.to_json()
with open("experiments/BinPackingBenchmarks.json", "w") as outfile:
    outfile.write(result)