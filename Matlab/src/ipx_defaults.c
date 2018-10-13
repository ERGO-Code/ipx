#include "mex.h"
#include "ipx_c.h"

static void set_double(mxArray *pm, const char *fieldname, double value)
{
    mxSetField(pm, 0, fieldname, mxCreateDoubleScalar(value));
}

static void set_string(mxArray *pm, const char *fieldname, const char *value)
{
    if (!value)
        mxSetField(pm, 0, fieldname, mxCreateString(""));
    else
        mxSetField(pm, 0, fieldname, mxCreateString(value));
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    const char *field_names[] = {
        "display",
        "logfile",
        "print_interval",
        "time_limit",
        "dualize",
        "scale",
        "residual_tol",
        "optimality_tol",
        "drop_primal",
        "drop_dual",
        "maxiter",
        "kkt_tol",
        "precond_dense_cols",
        "crash_basis",
        "dependency_tol",
        "volume_tol",
        "update_heuristic",
        "maxpasses",
        "maxskip_updates",
        "slices",
        "lu_kernel",
        "lu_pivottol",
        "crossover",
        "crossover_start",
        "pfeastol",
        "dfeastol",
        "debug",
        "switchiter",
        "stop_at_switch"
    };

    if (nrhs != 0)
        mexErrMsgTxt("No input argument required.");
    if (nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");
    if (nlhs == 0)
        return;

    mwSize dims[2] = {1, 1};
    int num_fields = sizeof(field_names)/sizeof(field_names[0]);
    plhs[0] = mxCreateStructArray(2, dims, num_fields, field_names);
    struct ipx_parameters defaults = ipx_default_parameters();
    set_double(plhs[0], "display",          defaults.display);
    set_string(plhs[0], "logfile",          defaults.logfile);
    set_double(plhs[0], "print_interval",   defaults.print_interval);
    set_double(plhs[0], "time_limit",       defaults.time_limit);
    set_double(plhs[0], "dualize",          defaults.dualize);
    set_double(plhs[0], "scale",            defaults.scale);
    set_double(plhs[0], "residual_tol",     defaults.residual_tol);
    set_double(plhs[0], "optimality_tol",   defaults.optimality_tol);
    set_double(plhs[0], "drop_primal",      defaults.drop_primal);
    set_double(plhs[0], "drop_dual",        defaults.drop_dual);
    set_double(plhs[0], "maxiter",          defaults.maxiter);
    set_double(plhs[0], "kkt_tol",          defaults.kkt_tol);
    set_double(plhs[0], "precond_dense_cols", defaults.precond_dense_cols);
    set_double(plhs[0], "crash_basis",      defaults.crash_basis);
    set_double(plhs[0], "dependency_tol",   defaults.dependency_tol);
    set_double(plhs[0], "volume_tol",       defaults.volume_tol);
    set_double(plhs[0], "update_heuristic", defaults.update_heuristic);
    set_double(plhs[0], "maxpasses",        defaults.maxpasses);
    set_double(plhs[0], "maxskip_updates",  defaults.maxskip_updates);
    set_double(plhs[0], "slices",           defaults.slices);
    set_double(plhs[0], "lu_kernel",        defaults.lu_kernel);
    set_double(plhs[0], "lu_pivottol",      defaults.lu_pivottol);
    set_double(plhs[0], "crossover",        defaults.crossover);
    set_double(plhs[0], "crossover_start",  defaults.crossover_start);
    set_double(plhs[0], "pfeastol",         defaults.pfeastol);
    set_double(plhs[0], "dfeastol",         defaults.dfeastol);
    set_double(plhs[0], "debug",            defaults.debug);
    set_double(plhs[0], "switchiter",       defaults.switchiter);
    set_double(plhs[0], "stop_at_switch",   defaults.stop_at_switch);
}
