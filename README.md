# Logic-Based Discrete-Steepest Descent Algorithm (L-DSDA)
This repository contains the code for the LD-SDA. The code here receives a Generalized Disjunctive Programming (GDP) problem modeled through [Pyomo.GDP](https://pyomo.readthedocs.io/en/latest/modeling_extensions/gdp/modeling.html) and performs the *external variable* reformulation such that the Discrete-Steepest Descent Algorithm, including a neighbor and a line search method in an integer lattice described by those external variables.
This algorithm is particularly well-suited for GDP problems with ordered sets and knapsack equality constraints, i.e., where from an ordered set of Boolean variables only a given number can be *True* at the same time.
This kind of problem structure arises very often in process superstructure optimization.

The implementation of the method can be found [here](gdp/dsda). To execute the problems and recreate the resuts run:

``python main_*.py``

Where * represents the problem you would like to solve.
This repository contains several examples.
- [CSTR in series network](gdp/cstr) : Example adjusted from [^1] of a series of CSTR reactors where the volume is minimized depending a parametric number of total possible reacotrs to install and the position of the reflux in the network. Has an analytical solution when the number of reactors tend to be infinite as the Plug Flow Reactor solution.
- [Distillation column](gdp/column) : Example from [^2] avaliable [here](https://github.com/grossmann-group/gdplib/tree/master/gdplib/gdp_col) of a distillation column design problem to determine the optimal number of trays and feed stage.
- [Small batch scheduling](gdp/smallbatch) : Example from [^3] avaliable [here](https://www.gams.com/33/gamslib_ml/libhtml/gamslib_batchdes.html) of a small batch scheduling process involving certain number of parallel unit operations.
- [Catalytic distillation column design](catalytic_distillation) : Equilibrium based example from [^4] and rate-based example from [^5] implemented in GAMS both for the traditional MINLP-based D-SDA and the Logic-based D-SDA.
References:
[^1]: Liñán, D. A., Bernal, D. E., Ricardez-Sandoval, L. A., & Gómez, J. M. (2020a). Optimal design of superstructures for placing units and streams with multiple and ordered available locations. Part I: A new mathematical framework. Computers & Chemical Engineering, 137, 106794. [link](https://doi.org/10.1016/j.compchemeng.2020.106794)
[^2]: Jaffer H. Ghouse, Qi Chen, Miguel A. Zamarripa, Andrew Lee, Anthony P. Burgard, Ignacio E. Grossmann, David C. Miller. A comparative study between GDP and NLP formulations for conceptual design of distillation columns. Computer Aided Chemical Engineering. Volume 44, 2018, Pages 865-870. [link](https://doi.org/10.1016/B978-0-444-64241-7.50139-7)
[^3]: Kocis, G R, and Grossmann, I E, Global Optimization of Nonconvex MINLP. Problems in Process Synthesis. Independent Engineering Chemical Research 27 (1988), 1407-1421. [link](https://doi.org/10.1021/ie00080a013)
[^4]: Liñán, D. A., Bernal, D. E., Ricardez-Sandoval, L. A., & Gómez, J. M. (2020). Optimal design of superstructures for placing units and streams with multiple and ordered available locations. Part II: Rigorous design of catalytic distillation columns. Computers & Chemical Engineering, 139, 106845. [link](https://doi.org/10.1016/j.compchemeng.2020.106845)
[^5]: Liñán, D. A., Bernal, D. E., Gómez, J. M., & Ricardez-Sandoval, L. A. (2021). Optimal synthesis and design of catalytic distillation columns: A rate-based modeling approach. Chemical Engineering Science, 231, 116294. [link](https://doi.org/10.1016/j.ces.2020.116294)
