import sys, os

fpathParent = os.path.join(os.path.dirname(__file__), '..')
fpathSolver = os.path.join(os.path.dirname(__file__), '../packingSolver')
sys.path.append(fpathParent)
sys.path.append(fpathSolver)

from packingSolver.PlacementPoints import *
from packingSolver.BinPackingData import *
from packingSolver.OrthogonalPacking import *

import re

import unittest
import timeit

class TestPackingProblems(unittest.TestCase):
    #def __init__(self, *args, **kwargs):
    #    super(TestPackingProblems, self).__init__(*args, **kwargs)
    def setUp(self):
        self.IsFeasibleDict = {
            'E00N10.json': False,
            'E00N15.json': False,
            #'E00N23.json': False, # takes too long
            #'E00X23.json': False, # takes way too long
            'E02F17.json': True,
            'E02F20.json': True,
            'E02F22.json': True,
            'E02N20.json': False,
            'E03N10.json': False,
            'E03N15.json': False,
            'E03N16.json': False,
            'E03N17.json': False,
            #'E03X18.json': ,
            'E04F15.json': True,
            'E04F17.json': True,
            'E04F19.json': True,
            'E04F20.json': True,
            'E04N15.json': False,
            'E04N17.json': False,
            'E04N18.json': False,
            'E05F15.json': True,
            'E05F18.json': True,
            'E05F20.json': True,
            'E05N15.json': False,
            'E05N17.json': False,
            #'E05X15.json': ,
            'E07F15.json': True,
            'E07N10.json': False,
            'E07N15.json': False,
            #'E07X15.json': ,
            'E08F15.json': True,
            'E08N15.json': False,
            'E10N10.json': False,
            'E10N15.json': False,
            #'E10X15.json': ,
            'E13N10.json': False,
            'E13N15.json': False,
            #'E13X15.json': ,
            'E15N10.json': False,
            'E15N15.json': False,
            #'E20X15.json': ,
            'E20F15.json': True
        }

    def TestOPP(self, path, fileName, placementPoints):
        items, H, W = ReadBenchmarkData(path, fileName)
        bin = Bin(W, H)

        solver = OrthogonalPackingSolver(items, bin, placementPoints)
        isFeasible = solver.Solve()

        return isFeasible

    def test_all_instances(self):
        path = 'data/input/OPP/CJCM08'
        for root, dirs, files in os.walk(path):
            for fileName in files:
                if fileName not in self.IsFeasibleDict:
                    continue

                for placementPointStrategy in PlacementPointStrategy:
                    isExpectedFeasible = self.IsFeasibleDict[fileName]
                    isFeasible = self.TestOPP(path, fileName, placementPointStrategy)
                    
                    isExpectedFeasibleString = 'feasible' if isExpectedFeasible else 'infeasible'
                    isFeasibleString = 'feasible' if isFeasible else 'infeasible'

                    message = f'Instance {fileName} with {placementPointStrategy.name} should be {isExpectedFeasibleString} but is {isFeasibleString}'

                    self.assertEqual(isFeasible, isExpectedFeasible, message)

def test_performance(benchmark):
    test = TestPackingProblems()
    test.setUp()
    result = benchmark(test.test_all_instances)

def main():
    unittest.main()

if __name__ == "__main__":
    main()