# Exact solutions to two-dimensional bin packing problems

This is a partial replication of [Côté, Haouari, Iori (2019): A Primal Decomposition Algorithm for the Two-dimensional Bin Packing Problem. arXiv:1909.06835 [math.OC]](https://arxiv.org/abs/1909.06835) resp. [Côté, Haouari, Iori (2021): Combinatorial Benders Decomposition for the Two-Dimensional Bin Packing Problem. INFORMS Journal on Computing.](https://doi.org/10.1287/ijoc.2020.1014). 

In contrast to the original paper, this implementation uses Constraint Programming to solve subproblems and generate initial solutions. Not all preprocessing routines and lower bounding procedures are implemented.
