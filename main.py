import json

from BinPackingData import *
from BranchAndCut import BinPackingBranchAndCutSolver

def main():
    solutions = {}
    # Single bin 2D-BPP CP model takes ages to prove feasibility/infeasibility on instances: 173, 262, 292, 297, 298, 317
    # Meet in the middle produces erroneous 1D KP instances for instances 11
    # Postsolve error: 51, 45, 125/126, 31, 26 or 27, 311
    # Double count: 21, 22, 26, 27, 31, 32
    # Negative lifting coefficient: 120, 123, 
    #hardInstances = [226, 232, 242, 242, 244, 245, 247, 248, 249, 261, 292, 313, 314, 332, 173, 187, 188, 191, 197, 142, 143, 145, 146, 149]
    #mediumInstance = [149, 174]

    for instance in range(120, 501):
    #for instance in hardInstances:
        currentInstanceId = instance
        items, H, W = ReadBenchmarkData(instance)
        
        solver = BinPackingBranchAndCutSolver(instance)
        rectangles = solver.Run(items, H, W)

        solver.RetrieveSolutionStatistics()

        # TODO: introduce solution statistics struct
        bestBoundMIP = solver.LB
        upperBoundMIP = solver.UB
        solverType = solver.SolverType
        isOptimalMIP = solver.IsOptimal
        
        #PlotSolution(upperBoundMIP * W, H, rectangles)

        if isOptimalMIP:
            print(f'Instance {instance}: Optimal solution = {int(bestBoundMIP)} found by {solverType} (#items = {len(items)})')
        else:
            print(f'Instance {instance}: No optimal solution found, [lb, ub] = [{bestBoundMIP}, {upperBoundMIP}] (#items = {len(items)})')
    
        solutions[instance] = {'LB': bestBoundMIP, 'UB': upperBoundMIP, 'Solver': solverType}

    solutionsJson = json.dumps(solutions, indent = 4)
    with open("Solutions.json", "w") as outfile:
        outfile.write(solutionsJson)

if __name__ == "__main__":
    main()