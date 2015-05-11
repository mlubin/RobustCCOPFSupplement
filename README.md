# RobustCCOPFSupplement
Repository containing supplementary data and code for "A Robust Approach to Chance Constrained Optimal Power Flow with Renewable Generation" by Lubin, Dvorkin, and Backhaus

*Installation instructions:*

The optimization model was implemented by using the [JuMPChance](https://github.com/mlubin/JuMPChance.jl) extension to [JuMP](https://github.com/JuliaOpt/JuMP.jl) in the [Julia](http://julialang.org/downloads/) programming language.
Additionaly, we used Gurobi 6.0 in our numerical experiments. [Gurobi](http://www.gurobi.com/) is a commercial solver which must be installed and licensed separately (one may easily use a different solver if Gurobi is not available, see the JuMP documentation).

The experiments require Julia 0.3 or later, and the following Julia packages:
- [JuMP](https://github.com/JuliaOpt/JuMP.jl) 0.9.0
- [JuMPChance](https://github.com/mlubin/JuMPChance.jl) 0.1.1
- [Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl) 0.1.26
- [MAT](https://github.com/simonster/MAT.jl)

You should force the use of particular versions of these Julia packages with
```
julia> Pkg.pin("JuMP", v"0.9.0")
julia> Pkg.pin("JuMPChance", v"0.1.1")
julia> Pkg.pin("Gurobi", v"0.1.26")
```

*Running the code:*

The code for the experiments is contained in the ``codejl`` directory. Note that it is under 700 lines of code. The file ``input.jl`` contains routines and data structures for processing the input, and the file ``ccopfmodel_simulation.jl`` (heavily commented) contains the main simulation logic and optimization model.

You can run the robust chance-constrained OPF model by entering the ``run`` directory and executing:
```
julia ../codejl/ccopfmodel_simulation.jl bpa-season1-gamma1.0.dat
```

By default, this will run a small sample (10 hours) from the first seasion used in the experiment. See line 220 of ``ccopfmodel_simulation.jl`` to adjust this. The output of the simulation is a ``.mat`` file (``ccopf_bpa_simulation_results_season1_gamma1.0.mat``) which can be opened directly in MATLAB or loaded in Julia via the MAT package. Contained in the ``.mat`` file are the optimal objective values, solution times, solution status, and optimal values for the decision variables. See lines 323 and below of ``ccopfmodel_simulation.jl`` for more details on the output format.

The file ``bpa-season1-gamma1.0.dat`` specifies all of the input paths and parameters for the simulation. In particular, one can modify the ``robust_budget`` to chance the budget of uncertainty Γ for the parameters of the Gaussian distributions, as discussed in the paper.

Further documentation is available on request.
