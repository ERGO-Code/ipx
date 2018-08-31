#ifndef IPX_PARAMETERS_H_
#define IPX_PARAMETERS_H_

#include "ipx_config.h"

struct ipx_parameters {
    /* Solver control */
    ipxint debug;
    ipxint display;
    const char* logfile;
    double print_interval;
    double time_limit;

    /* Preprocessing */
    ipxint dualize;
    ipxint scale;

    /* Interior point method */
    double residual_tol;
    double optimality_tol;
    double drop_primal;
    double drop_dual;
    ipxint maxiter;

    /* Linear solver */
    double kkt_tol;
    ipxint precond_dense_cols;
    ipxint switchiter;
    ipxint stop_at_switch;

    /* Basis construction in IPM */
    ipxint crash_basis;
    double dependency_tol;
    double volume_tol;
    ipxint update_heuristic;
    ipxint maxpasses;
    ipxint maxskip_updates;
    ipxint slices;

    /* LU factorization */
    ipxint lu_kernel;
    double lu_pivottol;

    /* Crossover */
    ipxint crossover;
    double crossover_start;
    double pfeastol;
    double dfeastol;

    #ifdef __cplusplus
    ipx_parameters() {
        debug = 0;
        display = 1;
        logfile = nullptr;
        print_interval = 5.0;
        time_limit = -1.0;
        dualize = -1;
        scale = 1;
        residual_tol = 1e-6;
        optimality_tol = 1e-8;
        drop_primal = 1e-9;
        drop_dual = 1e-9;
        maxiter = 300;
        kkt_tol = 0.3;
        precond_dense_cols = 1;
        switchiter = -1;
        stop_at_switch = 0;
        crash_basis = 1;
        dependency_tol = 1e-7;
        volume_tol = 2.0;
        update_heuristic = 1;
        maxpasses = -1;
        maxskip_updates = 10;
        slices = -1;
        lu_kernel = 0;
        lu_pivottol = 0.0625;
        crossover = 1;
        crossover_start = 1e-8;
        pfeastol = 1e-7;
        dfeastol = 1e-7;
    }
    #endif
};

#endif  /* IPX_PARAMETERS_H_ */
