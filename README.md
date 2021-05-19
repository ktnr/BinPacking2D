# Exact solutions to two-dimensional bin packing problems

This is a partial replication of [Côté, Haouari, Iori (2019): A Primal Decomposition Algorithm for the Two-dimensional Bin Packing Problem. arXiv:1909.06835 [math.OC]](https://arxiv.org/abs/1909.06835) resp. [Côté, Haouari, Iori (2021): Combinatorial Benders Decomposition for the Two-Dimensional Bin Packing Problem. INFORMS Journal on Computing.](https://doi.org/10.1287/ijoc.2020.1014)

In contrast to the original paper, this implementation uses constraint programming to solve subproblems and generate initial solutions, instead of mixed-integer programms (MIPs) or specialized heuristics. Not all preprocessing routines and lower bounding procedures are implemented.

## Algorithmic outline

The two-dimensional bin packing problem (2D-BPP) is decomposed into a master problem and several subproblems. The master problem is a (one-dimensional) BPP, which assigns items to bins. Feasibility of each bin assignment is then checked in a separate two-dimensional knapsack subproblem.

The master problem is modeled as a MIP and solved using [Gurobi](https://www.gurobi.com/). Callbacks are set at integer nodes, where the subproblems are generated and solved by constraint programming (CP) using [google or-tools](https://developers.google.com/optimization/cp/cp_solver). When an infeasible subproblem is found, a cut of the form

<img src="https://render.githubusercontent.com/render/math?math=\sum_{j \in \mathcal{C}} x_{ij} \leq |\mathcal{C}| - 1, \quad i \in \mathcal{B}">

is added to exclude the infeasible assignment of items <img src="https://render.githubusercontent.com/render/math?math=j \in \mathcal{C}"> to bins <img src="https://render.githubusercontent.com/render/math?math=i \in \mathcal{B}">.

- Start solutions are generated by a solving a 2D-BPP problem via CP; a simpler version of this model can be found [here](http://yetanothermathprogrammingconsultant.blogspot.com/2021/02/2d-bin-packing-with-google-or-tools-cp.html)
- Preprocessing routines include: filtering of large items, determining incompatible items, fixing and restricting item assignment through a conflict graph based on item incompatility, different placement point strategies (normal patterns, meet-in-the-middle) 
- Improvements to the decomposition approach include: cut strengthening by finding reduced infeasible sets (decrease lhs and rhs of infeasibility cuts), cut lifting by finding valid item displacements for infeasible sets (increasing lhs while keeping rhs constant)

## Usage
The entry point is the `main()`method of `BinPacking.py`. Reading, writing and displaying data is handled by `BinPackingData.py`. By default the benchmark sets of [J. O. Berkey and P. Y. Wang, “Two-Dimensional Finite Bin-Packing Algorithms,” J. Oper. Res. Soc., vol. 38, no. 5, p. 423, May 1987 and Martello and D. Vigo, “Exact Solution of the Two-Dimensional Finite Bin Packing Problem,” Manage. Sci., vol. 44, no. April 2015, pp. 388–399, 1998](https://github.com/Oscar-Oliveira/OR-Datasets/tree/master/Cutting-and-Packing/2D/Datasets/CLASS) are included.

The high level algorithm is implemented in the `BinPackingBranchAndCutSolver` class.
```
itemHeights, itemWidths, binHeight, bindWidth, numberOfItems = ReadBenchmarkData(instance)
        
solver = BinPackingBranchAndCutSolver()
rectangles = solver.Run(itemHeights, itemWidths, binHeight, bindWidth, numberOfItems)
```

Algorithmic features can be parametrized.
```
def SetCallbackData(self):
    self.Model._EnableLifting = False
    self.Model._EnableCutStrengthening = True
    self.Model._InfeasibleSuproblemCutThreshold = 1

    self.Model._EnablePreprocessLifting = False
    self.Model._PlacementPointStrategy = PlacementPointStrategy.MinimalMeetInTheMiddlePatterns
```

## Results
The output indicates whether the solution is optimal and if the solution was found during the constructive phase by constraint programming (solving a 2D-BPP with google's CP-SAT solver) or by branch-and-cut.
```
Instance 1: Optimal solution = 8 found by CP (#items = 20)
Instance 2: Optimal solution = 5 found by CP (#items = 20)
Instance 3: Optimal solution = 9 found by CP (#items = 20)
Instance 4: Optimal solution = 6 found by CP (#items = 20)
Instance 5: Optimal solution = 6 found by CP (#items = 20)
Instance 6: Optimal solution = 9 found by CP (#items = 20)
Instance 7: Optimal solution = 6 found by CP (#items = 20)
Instance 8: Optimal solution = 6 found by CP (#items = 20)
Instance 9: Optimal solution = 8 found by CP (#items = 20)
Instance 10: Optimal solution = 8 found by CP (#items = 20)
Instance 11: Optimal solution = 10 found by B&C (#items = 40)
```

## Evaluation
The algorithm produces optimal solutions for a majority of the 500 benchmark instances in less than 20 minutes. It has difficulty in proving feasibility/infeasibility of 2D-Knapsack subproblems for instances with many small items (e.g. google's CP-SAT solver takes more than 1000 seconds to solve a single 2D-Knapsack problem for benchmark instance 173). 

By far the most impactful algorithmic components are
- symmetry breaking constraints (24) and (25) of the original paper (arXiv version),
- removing large items and fixing conflicting items to bins (27) in accordance with (24) and (25), 
- and excluding pairwise incompatible items from bins with fixed items (26),
- no-good cuts of type (4),
- strong constraint programming formulations for the start solution (2D-BPP) and subproblems (2D-Knapsack) by reducing decision variable domains as much as possible, i.e. reducing available placement points through the techniques mentioned above and placement point patterns such as the meet-in-the-middle patterns.

The benefit of producing no-good cuts (29) from reduced feasible sets is only marginal. The benefit of lifting these cuts (30) is also only marginal, mainly due to the numerous 2D-Knapsack problems that must be solved. Hence, speeding up the solution of the 2D-Knapsack problem of Algorithm 2 (currently solved as 2D-Knapsack with google's CP-SAT and a time limit of 1 second) might increase the impact of lifting. 

Components with the greatest potential to improve solution times:
- an efficient algorithm to prove feasibility/infeasibility of problematic 2D-Knapsack subproblems with numerous small items
- a variant of the BKRS lower bound (23)
