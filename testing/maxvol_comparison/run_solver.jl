function RunSolver(name::String, group::String)
    model = Mathprog.LPget(name, group)
    solver = ipx.LPSolver()
    ipx.SetParameter(solver, :debug, 2)
    ipx.SetParameter(solver, :crossover, 0)
    ipx.SetParameter(solver, :ipm_drop_primal, 0.0)
    ipx.SetParameter(solver, :ipm_drop_primal, 0.0)
    volume_tol = [10.0; 4.0; 2.0; 1.1]

    # Maxvolume heuristic.
    ipx.SetParameter(solver, :update_heuristic, 1)
    for v in volume_tol
        ipx.SetParameter(solver, :volume_tol, v)
        logfile = string(name, "_heuristic_", @sprintf("%05.1f", v), ".log")
        logfile = joinpath("logs", logfile)
        ipx.Solve(solver, model, logfile)
    end

    # True maxvolume basis.
    ipx.SetParameter(solver, :update_heuristic, 0)
    ipx.SetParameter(solver, :maxpasses, -1)
    for v in volume_tol
        ipx.SetParameter(solver, :volume_tol, v)
        logfile = string(name, "_sequential_", @sprintf("%05.1f", v), ".log")
        logfile = joinpath("logs", logfile)
        ipx.Solve(solver, model, logfile)
    end
end
