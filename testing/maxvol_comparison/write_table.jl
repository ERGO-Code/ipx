"""
  WriteTable(probname_method)

  Writes table from all logfiles whose filename starts with probname_method.
"""
function WriteTable(pattern::String)
    tblfile = joinpath("results", string(pattern, ".tbl"))
    tblfile = open(tblfile, "w")
    @printf(tblfile, " %12s", "maxvol_tol")
    @printf(tblfile, " %12s", "updates_ipm")
    @printf(tblfile, " %12s", "time_maxvol")
    @printf(tblfile, " %12s", "kktiter2")
    @printf(tblfile, " %12s", "time_cr2")
    @printf(tblfile, " %12s", "time_ipm2")
    @printf(tblfile, "\n")
    files = readdir("logs")
    for fname in files
        fbase, fext = splitext(fname)
        r = Regex(string("^", pattern))
        if fext == ".log" && match(r, fname) != nothing
            ipxinfo = ipx.ParseInfo(joinpath("logs", fname))
            maxvol_tol = match(r"_([.\d]+)$", fbase)[1]
            maxvol_tol = parse(Float64, maxvol_tol)
            @printf(tblfile, " %12.1f", maxvol_tol)
            @printf(tblfile, " %12d", ipxinfo.updates_ipm)
            @printf(tblfile, " %12.2f", ipxinfo.time_maxvol)
            @printf(tblfile, " %12d", ipxinfo.kktiter2)
            @printf(tblfile, " %12.2f", ipxinfo.time_cr2)
            @printf(tblfile, " %12.2f", ipxinfo.time_ipm2)
            @printf(tblfile, "\n")
        end
    end
    close(tblfile)
end
