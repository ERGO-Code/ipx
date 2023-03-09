// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "lp_solver.h"
#include <algorithm>
#include <cassert>
#include <vector>
#include <utility>
#include "crossover.h"
#include "info.h"
#include "kkt_solver_basis.h"
#include "kkt_solver_diag.h"
#include "starting_basis.h"
#include "utils.h"
#include "version.h"

namespace ipx {

Int LpSolver::LoadModel(Int num_var, const double* obj, const double* lb,
                        const double* ub, Int num_constr, const Int* Ap,
                        const Int* Ai, const double* Ax, const double* rhs,
                        const char* constr_type) {
    ClearModel();
    Int errflag = user_model_.Load(control_, num_constr, num_var, Ap, Ai, Ax,
                                   rhs, constr_type, obj, lb, ub);
    if (errflag)
        return errflag;
    user_model_.GetInfo(&info_);
    return 0;
}

Int LpSolver::LoadIPMStartingPoint(const double* x, const double* xl,
                                   const double* xu, const double* slack,
                                   const double* y, const double* zl,
                                   const double* zu) {
    ClearIPMStartingPoint();
    Int errflag = user_model_.CheckInteriorPoint(x, xl, xu, slack, y, zl, zu);
    if (errflag)
        return errflag;

    const Int num_var = user_model_.num_var();
    const Int num_constr = user_model_.num_constr();
    ipm_start_.reset(new InteriorSolution(num_var, num_constr));
    std::copy_n(x, num_var, std::begin(ipm_start_->x));
    std::copy_n(xl, num_var, std::begin(ipm_start_->xl));
    std::copy_n(xu, num_var, std::begin(ipm_start_->xu));
    std::copy_n(slack, num_constr, std::begin(ipm_start_->slack));
    std::copy_n(y, num_constr, std::begin(ipm_start_->y));
    std::copy_n(zl, num_var, std::begin(ipm_start_->zl));
    std::copy_n(zu, num_var, std::begin(ipm_start_->zu));
    return 0;
}

Int LpSolver::Solve() {
    ClearSolution();
    if (user_model_.empty())
        return info_.status = IPX_STATUS_no_model;
    control_.ResetTimer();
    control_.OpenLogfile();
    control_.Log() << "IPX version " << Version{} << '\n';
    try {
        Presolve();
        if (info_.status == IPX_STATUS_not_run) {
            InteriorPointSolve();
            if ((info_.status_ipm == IPX_STATUS_optimal ||
                 info_.status_ipm == IPX_STATUS_imprecise) &&
                control_.crossover())
                RunCrossover();
        }
        if (basis_) {
            info_.ftran_sparse = basis_->frac_ftran_sparse();
            info_.btran_sparse = basis_->frac_btran_sparse();
            info_.time_lu_invert = basis_->time_factorize();
            info_.time_lu_update = basis_->time_update();
            info_.time_ftran = basis_->time_ftran();
            info_.time_btran = basis_->time_btran();
            info_.mean_fill = basis_->mean_fill();
            info_.max_fill = basis_->max_fill();
        }
        if (info_.status_ipm == IPX_STATUS_primal_infeas ||
            info_.status_ipm == IPX_STATUS_dual_infeas ||
            info_.status_crossover == IPX_STATUS_primal_infeas ||
            info_.status_crossover == IPX_STATUS_dual_infeas) {
            // When IPM or crossover detect the model to be infeasible
            // (currently only the former is implemented), then the problem is
            // solved.
            info_.status = IPX_STATUS_solved;
        } else {
            Int method_status = control_.crossover() ?
                info_.status_crossover : info_.status_ipm;
            if (method_status == IPX_STATUS_optimal ||
                method_status == IPX_STATUS_imprecise)
                info_.status = IPX_STATUS_solved;
            else
                info_.status = IPX_STATUS_stopped;
        }
        PrintSummary();
    }
    catch (const std::bad_alloc&) {
        control_.Log() << " out of memory\n";
        info_.status = IPX_STATUS_out_of_memory;
    }
    catch (const std::exception& e) {
        control_.Log() << " internal error: " << e.what() << '\n';
        info_.status = IPX_STATUS_internal_error;
    }
    info_.time_total = control_.Elapsed();
    control_.Debug(2) << info_;
    control_.CloseLogfile();
    return info_.status;
}

Info LpSolver::GetInfo() const {
    return info_;
}

namespace {
    void copy_if_not_null(const Vector& src, double* dest) {
        if (dest)
            std::copy(std::begin(src), std::end(src), dest);
    }
    void copy_if_not_null(const std::vector<Int>& src, Int* dest) {
        if (dest)
            std::copy(std::begin(src), std::end(src), dest);
    }
} // namespace

Int LpSolver::GetInteriorSolution(InteriorSolution& solution) const {
    if (!interior_solution_)
        return -1;
    solution = *interior_solution_;
    return 0;
}

Int LpSolver::GetInteriorSolution(double* x, double* xl, double* xu,
                                  double* slack, double* y, double* zl,
                                  double* zu) const {
    if (!interior_solution_)
        return -1;
    copy_if_not_null(interior_solution_->x, x);
    copy_if_not_null(interior_solution_->xl, xl);
    copy_if_not_null(interior_solution_->xu, xu);
    copy_if_not_null(interior_solution_->slack, slack);
    copy_if_not_null(interior_solution_->y, y);
    copy_if_not_null(interior_solution_->zl, zl);
    copy_if_not_null(interior_solution_->zu, zu);
    return 0;
}

Int LpSolver::GetBasicSolution(BasicSolution& solution) const {
    if (!basic_solution_)
        return -1;
    solution = *basic_solution_;
    return 0;
}

Int LpSolver::GetBasicSolution(double* x, double* slack, double* y, double* z,
                               Int* cbasis, Int* vbasis) const {
    if (!basic_solution_)
        return -1;
    copy_if_not_null(basic_solution_->x, x);
    copy_if_not_null(basic_solution_->slack, slack);
    copy_if_not_null(basic_solution_->y, y);
    copy_if_not_null(basic_solution_->z, z);
    copy_if_not_null(basic_solution_->cbasis, cbasis);
    copy_if_not_null(basic_solution_->vbasis, vbasis);
    return 0;
}

Parameters LpSolver::GetParameters() const {
    return control_.parameters();
}

void LpSolver::SetParameters(Parameters new_parameters) {
    control_.parameters(new_parameters);
}

Int LpSolver::ReadParameters(const char* filename) {
    return control_.ReadParameters(filename);
}

Int LpSolver::WriteParameters(const char* filename) {
    return control_.WriteParameters(filename);
}

void LpSolver::ClearModel() {
    user_model_.clear();
    model_.clear();
    presolver_.clear();
    ClearSolution();
    ClearIPMStartingPoint();
}

void LpSolver::ClearIPMStartingPoint() {
    ipm_start_.reset(nullptr);
}

Int LpSolver::GetIterate(double* x, double* y, double* zl, double* zu,
                         double* xl, double* xu) {
    if (!iterate_)
        return -1;
    copy_if_not_null(iterate_->x(), x);
    copy_if_not_null(iterate_->y(), y);
    copy_if_not_null(iterate_->zl(), zl);
    copy_if_not_null(iterate_->zu(), zu);
    copy_if_not_null(iterate_->xl(), xl);
    copy_if_not_null(iterate_->xu(), xu);
    return 0;
}

// Returns a vector of basic statuses that is consistent with the basis and
// the bounds from the model.
static std::vector<Int> BuildBasicStatuses(const Basis& basis) {
    const Model& model = basis.model();
    const Int m = model.rows();
    const Int n = model.cols();
    const Vector& lb = model.lb();
    const Vector& ub = model.ub();
    std::vector<Int> basic_statuses(n+m);
    for (Int j = 0; j < n+m; j++) {
        if (basis.IsBasic(j)) {
            basic_statuses[j] = IPX_basic;
        } else if (std::isfinite(lb[j])) {
            basic_statuses[j] = IPX_nonbasic_lb;
        } else if (std::isfinite(ub[j])) {
            basic_statuses[j] = IPX_nonbasic_ub;
        } else {
            basic_statuses[j] = IPX_superbasic;
        }
    }
    return basic_statuses;
}

Int LpSolver::GetBasis(Int* cbasis, Int* vbasis) {
    if (!basis_)
        return -1;
    if (basic_solution_) {
        // crossover provides basic statuses
        copy_if_not_null(basic_solution_->cbasis, cbasis);
        copy_if_not_null(basic_solution_->vbasis, vbasis);
    } else {
        presolver_.PostsolveBasis(BuildBasicStatuses(*basis_), cbasis, vbasis);
    }
    return 0;
}

Int LpSolver::GetKKTMatrix(Int* AIp, Int* AIi, double* AIx, double* g) {
    if (!iterate_)
        return -1;
    if (AIp && AIi && AIx) {
        const SparseMatrix& AI = model_.AI();
        std::copy_n(AI.colptr(), AI.cols()+1, AIp);
        std::copy_n(AI.rowidx(), AI.entries(), AIi);
        std::copy_n(AI.values(), AI.entries(), AIx);
    }
    if (g) {
        Int m = model_.rows();
        Int n = model_.cols();
        for (Int j = 0; j < n+m; j++) {
            switch (iterate_->StateOf(j)) {
            case Iterate::State::fixed:
                g[j] = INFINITY;
                break;
            case Iterate::State::free:
                g[j] = 0.0;
                break;
            case Iterate::State::barrier:
                g[j] = iterate_->zl(j)/iterate_->xl(j) +
                    iterate_->zu(j)/iterate_->xu(j);
                assert(std::isfinite(g[j]));
                assert(g[j] > 0.0);
                break;
            default:
                assert(0);
            }
        }
    }
    return 0;
}

Int LpSolver::SymbolicInvert(Int* rowcounts, Int* colcounts) {
    if (!basis_)
        return -1;
    basis_->SymbolicInvert(rowcounts, colcounts);
    return 0;
}

void LpSolver::ClearSolution() {
    iterate_.reset(nullptr);
    basis_.reset(nullptr);
    simplex_iterate_.reset(nullptr);
    interior_solution_.reset(nullptr);
    basic_solution_.reset(nullptr);
    info_ = Info();
    // Restore info entries that belong to model.
    user_model_.GetInfo(&info_);
}

void LpSolver::Presolve() {
    Int errflag = presolver_.PresolveModel(control_);
    assert(errflag == 0);
    model_.GetInfo(&info_);
    presolver_.GetInfo(&info_);
}

void LpSolver::InteriorPointSolve() {
    control_.Log() << "Interior Point Solve\n";

    // Allocate new iterate and set tolerances for IPM termination test.
    iterate_.reset(new Iterate(model_));
    iterate_->feasibility_tol(control_.ipm_feasibility_tol());
    iterate_->optimality_tol(control_.ipm_optimality_tol());
    if (control_.crossover())
        iterate_->crossover_start(control_.crossover_start());

    RunIPM();

    iterate_->Postprocess();

    interior_solution_.reset(
        new InteriorSolution(user_model_.num_var(), user_model_.num_constr()));
    presolver_.PostsolveInteriorSolution(*iterate_, *interior_solution_);
    user_model_.EvaluateInteriorPoint(*interior_solution_, &info_);

    // Declare status_ipm "imprecise" if the IPM terminated optimal but the
    // solution after postprocessing/postsolve does not satisfy tolerances.
    if (info_.status_ipm == IPX_STATUS_optimal) {
        if (std::abs(info_.rel_objgap) > control_.ipm_optimality_tol() ||
            info_.rel_presidual > control_.ipm_feasibility_tol() ||
            info_.rel_dresidual > control_.ipm_feasibility_tol())
            info_.status_ipm = IPX_STATUS_imprecise;
    }
}

void LpSolver::RunIPM() {
    IPM ipm(control_);

    if (ipm_start_ && !model_.dualized()) {
        control_.Log() << " Using starting point provided by user."
            " Skipping initial iterations.\n";
        LoadStartingPoint(ipm);
        if (info_.status_ipm != IPX_STATUS_not_run)
            return;
    }
    else {
        if (ipm_start_)
            control_.Log() << " Ignoring starting point provided by user"
                " because presolver dualized model.\n";
        ComputeStartingPoint(ipm);
        if (info_.status_ipm != IPX_STATUS_not_run)
            return;
        RunInitialIPM(ipm);
        if (info_.status_ipm != IPX_STATUS_not_run)
            return;
    }
    BuildStartingBasis();
    if (info_.status_ipm != IPX_STATUS_not_run)
        return;
    RunMainIPM(ipm);
}

void LpSolver::LoadStartingPoint(IPM& ipm) {
    const Int m = model_.rows();
    const Int n = model_.cols();
    Vector x(n+m), xl(n+m), xu(n+m), y(m), zl(n+m), zu(n+m);

    Int errflag = presolver_.PresolveIPMStartingPoint(*ipm_start_,
                                                      x, xl, xu, y, zl, zu);
    assert(errflag == 0);

    ipm.LoadStartingPoint(x, xl, xu, y, zl, zu, iterate_.get(), &info_);
}

void LpSolver::ComputeStartingPoint(IPM& ipm) {
    Timer timer;
    KKTSolverDiag kkt(control_, model_);

    // If the starting point procedure fails, then iterate_ remains as
    // initialized by the constructor, which is a valid state for
    // postprocessing/postsolving.
    ipm.ComputeStartingPoint(&kkt, iterate_.get(), &info_);
    info_.time_ipm1 += timer.Elapsed();
}

void LpSolver::RunInitialIPM(IPM& ipm) {
    Timer timer;
    KKTSolverDiag kkt(control_, model_);

    Int switchiter = control_.switchiter();
    if (switchiter < 0) {
        // Switch iteration not specified by user. Run as long as KKT solver
        // converges within min(500,10+m/20) iterations.
        Int m = model_.rows();
        kkt.maxiter(std::min(500l, (long) (10+m/20) ));
        ipm.maxiter(control_.ipm_maxiter());
    } else {
        ipm.maxiter(std::min(switchiter, control_.ipm_maxiter()));
    }
    ipm.Driver(&kkt, iterate_.get(), &info_);
    switch (info_.status_ipm) {
    case IPX_STATUS_optimal:
        // If the IPM reached its termination criterion in the initial
        // iterations (happens rarely), we still call the IPM again with basis
        // preconditioning. In exact arithmetic it would terminate without an
        // additional iteration. A starting basis is then available for
        // crossover.
        info_.status_ipm = IPX_STATUS_not_run;
        break;
    case IPX_STATUS_no_progress:
        info_.status_ipm = IPX_STATUS_not_run;
        break;
    case IPX_STATUS_failed:
        info_.status_ipm = IPX_STATUS_not_run;
        info_.errflag = 0;
        break;
    case IPX_STATUS_iter_limit:
        if (info_.iter < control_.ipm_maxiter()) // stopped at switchiter
            info_.status_ipm = IPX_STATUS_not_run;
    }
    info_.time_ipm1 += timer.Elapsed();
}

void LpSolver::BuildStartingBasis() {
    if (control_.stop_at_switch() < 0) {
        info_.status_ipm = IPX_STATUS_debug;
        return;
    }
    basis_.reset(new Basis(control_, model_));
    control_.Log() << " Constructing starting basis...\n";
    StartingBasis(iterate_.get(), basis_.get(), &info_);
    if (info_.errflag == IPX_ERROR_interrupt_time) {
        info_.errflag = 0;
        info_.status_ipm = IPX_STATUS_time_limit;
        return;
    } else if (info_.errflag) {
        info_.status_ipm = IPX_STATUS_failed;
        return;
    }
    if (model_.dualized()) {
        std::swap(info_.dependent_rows, info_.dependent_cols);
        std::swap(info_.rows_inconsistent, info_.cols_inconsistent);
    }
    if (control_.stop_at_switch() > 0) {
        info_.status_ipm = IPX_STATUS_debug;
        return;
    }
    if (info_.rows_inconsistent) {
        info_.status_ipm = IPX_STATUS_primal_infeas;
        return;
    }
    if (info_.cols_inconsistent) {
        info_.status_ipm = IPX_STATUS_dual_infeas;
        return;
    }
}

void LpSolver::RunMainIPM(IPM& ipm) {
    KKTSolverBasis kkt(control_, *basis_);
    Timer timer;
    ipm.maxiter(control_.ipm_maxiter());
    ipm.Driver(&kkt, iterate_.get(), &info_);
    info_.time_ipm2 = timer.Elapsed();
}

void LpSolver::RunCrossover() {
    control_.Log() << "Crossover\n";
    assert(basis_);
    const Int m = model_.rows();
    const Int n = model_.cols();
    const Vector& lb = model_.lb();
    const Vector& ub = model_.ub();

    // Construct a complementary primal-dual point from the final IPM iterate.
    // This usually increases the residuals to Ax=b and A'y+z=c.
    simplex_iterate_.reset(new SimplexIterate(model_));
    iterate_->DropToComplementarity(*simplex_iterate_);

    // Run crossover. Perform dual pushes in increasing order and primal pushes
    // in decreasing order of the scaling factors from the final IPM iterate.
    {
        Vector weights(n+m);
        for (Int j = 0; j < n+m; j++)
            weights[j] = iterate_->ScalingFactor(j);
        Crossover crossover(control_);
        crossover.PushAll(basis_.get(), *simplex_iterate_, &weights[0], &info_);
        info_.time_crossover =
            crossover.time_primal() + crossover.time_dual();
        info_.updates_crossover =
            crossover.primal_pivots() + crossover.dual_pivots();
        if (info_.status_crossover != IPX_STATUS_optimal) {
            // Crossover failed. Discard solution.
            simplex_iterate_.reset(nullptr);
            return;
        }
    }

    // Recompute vertex solution and set basic statuses.
    basis_->ComputeBasicSolution(*simplex_iterate_);
    std::vector<Int> basic_statuses(n+m);
    for (Int j = 0; j < basic_statuses.size(); j++) {
        if (basis_->IsBasic(j)) {
            basic_statuses[j] = IPX_basic;
        } else {
            if (lb[j] == ub[j])
                basic_statuses[j] = simplex_iterate_->z[j] >= 0.0 ?
                    IPX_nonbasic_lb : IPX_nonbasic_ub;
            else if (simplex_iterate_->x[j] == lb[j])
                basic_statuses[j] = IPX_nonbasic_lb;
            else if (simplex_iterate_->x[j] == ub[j])
                basic_statuses[j] = IPX_nonbasic_ub;
            else
                basic_statuses[j] = IPX_superbasic;
        }
    }
    control_.Debug()
        << Textline("Bound violation of basic solution:")
        << sci2(PrimalInfeasibility(model_, simplex_iterate_->x)) << '\n'
        << Textline("Dual sign violation of basic solution:")
        << sci2(DualInfeasibility(
                    model_, simplex_iterate_->x, simplex_iterate_->z)) << '\n';
    control_.Debug()
        << Textline("Minimum singular value of basis matrix:")
        << sci2(basis_->MinSingularValue()) << '\n';

    basic_solution_.reset(
        new BasicSolution(user_model_.num_var(), user_model_.num_constr()));
    presolver_.PostsolveBasicSolution(*simplex_iterate_, basic_statuses,
                                      *basic_solution_);
    presolver_.PostsolveBasis(basic_statuses, &basic_solution_->cbasis[0],
                              &basic_solution_->vbasis[0]);
    user_model_.EvaluateBasicPoint(*basic_solution_, &info_);

    // Declare crossover status "imprecise" if the vertex solution defined by
    // the final basis does not satisfy tolerances.
    if (info_.primal_infeas > control_.pfeasibility_tol() ||
        info_.dual_infeas > control_.dfeasibility_tol())
        info_.status_crossover = IPX_STATUS_imprecise;
}

void LpSolver::PrintSummary() {
    control_.Log() << "Summary\n"
                   << Textline("Runtime:") << fix2(control_.Elapsed()) << "s\n"
                   << Textline("Status interior point solve:")
                   << StatusString(info_.status_ipm) << '\n'
                   << Textline("Status crossover:")
                   << StatusString(info_.status_crossover) << '\n';
    if (info_.status_ipm == IPX_STATUS_optimal ||
        info_.status_ipm == IPX_STATUS_imprecise) {
        control_.Log()
            << Textline("objective value:") << sci8(info_.pobjval) << '\n'
            << Textline("interior solution primal residual (abs/rel):")
            << sci2(info_.abs_presidual) << " / " << sci2(info_.rel_presidual)
            << '\n'
            << Textline("interior solution dual residual (abs/rel):")
            << sci2(info_.abs_dresidual) << " / " << sci2(info_.rel_dresidual)
            << '\n'
            << Textline("interior solution objective gap (abs/rel):")
            << sci2(info_.pobjval-info_.dobjval) << " / "
            << sci2(info_.rel_objgap)  << '\n';
    }
    if (info_.status_crossover == IPX_STATUS_optimal ||
        info_.status_crossover == IPX_STATUS_imprecise) {
        control_.Log()
            << Textline("basic solution primal infeasibility:")
            << sci2(info_.primal_infeas) << '\n'
            << Textline("basic solution dual infeasibility:")
            << sci2(info_.dual_infeas) << '\n';
    }
}

}  // namespace ipx
