This directory contains logfiles from LP solvers and scripts to evaluate them.
The information from the logfiles has been the source for the numerical
comparisons in [1,2].

The files in the results/ directory have been generated by the following Julia
commands:

julia> include("scripts/tblwriter.jl");

# Writes table with problem dimensions and runtimes for each of the three test sets.
julia> tblwriter.Results("testsets/diverse.tbl","logfiles/diverse/","results/diverse.tbl",false);
julia> tblwriter.Results("testsets/srd.tbl","logfiles/srd/","results/srd.tbl");
julia> tblwriter.Results("testsets/nug.tbl","logfiles/nug/","results/nug.tbl");

# Writes table comparing the geometric mean of the runtimes on the diverse test set.
julia> tblwriter.Comparison("testsets/diverse.tbl","logfiles/diverse/","results/diverse_comparison.tbl");

# Writes table with number of basis updates for the srd and nug models.
julia> tblwriter.BasisUpdates("testsets/srd.tbl","logfiles/srd/","results/srd_updates.tbl");
julia> tblwriter.BasisUpdates("testsets/nug.tbl","logfiles/nug/","results/nug_updates.tbl");

julia> include("scripts/matwriter.jl");

# Writes Matlab binary file with runtimes for different parts of the IPX algorithm.
julia> matwriter.IPXRuntimes("testsets/diverse.tbl","logfiles/diverse/","results/diverse_runtimes.mat");


The figure showing the fraction of IPX runtime spend in diiferent parts of the
algorithm has been generated by the following Matlab commands:

>> addpath('scripts');
>> load results/diverse_runtimes.mat
>> plot_runtime_breakdown


[1] L. Schork, J. Gondzio, "Implementation of an Interior Point Method with
    Basis Preconditioning", University of Edinburgh, School of Mathematics,
    Technical Report ERGO-18-014 (2018)

[2] L. Schork, "Basis Preconditioning in Interior Point Methods", PhD thesis,
    University of Edinburgh (2018)
