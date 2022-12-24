// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "control.h"
#include <algorithm>
#include <cctype>
#include <iostream>
#include <iterator>

namespace ipx {

Control::Control() {
    // When failbit is set, the stream evaluates to false.
    dummy_.setstate(std::ios::failbit);
    MakeStream();
}

Int Control::InterruptCheck() const {
    if (parameters_.time_limit >= 0.0 &&
        parameters_.time_limit < timer_.Elapsed())
        return IPX_ERROR_interrupt_time;
    return 0;
}

std::ostream& Control::Log() const {
    return output_;
}

std::ostream& Control::IntervalLog() const {
    if (parameters_.print_interval >= 0.0 &&
        interval_.Elapsed() >= parameters_.print_interval) {
        interval_.Reset();
        return output_;
    } else {
        return dummy_;
    }
}

std::ostream& Control::Debug(Int level) const {
    if (parameters_.debug >= level)
        return output_;
    else
        return dummy_;
}

void Control::ResetPrintInterval() const {
    interval_.Reset();
}

double Control::Elapsed() const {
    return timer_.Elapsed();
}

const Parameters& Control::parameters() const {
    return parameters_;
}

void Control::parameters(const Parameters& new_parameters) {
    parameters_ = new_parameters;
    MakeStream();
}

namespace {

bool isspace(char c) {
    return std::isspace(c);
}

// Extracts and trims two substrings from s, separated by whitespace.
// Returns true when exactly two substrings existed.
bool split_two(const std::string& s, std::string& s1, std::string& s2) {
    s1.clear();
    s2.clear();
    const auto beg1 = std::find_if_not(s.begin(), s.end(), isspace);
    if (beg1 == s.end())
        return false;
    const auto end1 = std::find_if(beg1 + 1, s.end(), isspace);
    const auto beg2 = std::find_if_not(end1, s.end(), isspace);
    if (end1 == s.end())
        return false;
    const auto end2 = std::find_if(beg2 + 1, s.end(), isspace);
    if (!std::all_of(end2, s.end(), isspace))
        return false;
    std::copy(beg1, end1, std::back_inserter(s1));
    std::copy(beg2, end2, std::back_inserter(s2));
    return true;
}

// Parses str as Int value. Returns true on success.
bool parse(const std::string& str, Int& result) {
    try {
        std::size_t pos{0};
        result = std::stoi(str, &pos);
        if (pos != str.size())
            throw std::invalid_argument{""};
    }
    catch (std::invalid_argument const&) {
        return false;
    }
    catch (std::out_of_range const&) {
        return false;
    }
    return true;
}

// Parses str as double value. Returns true on success.
bool parse(const std::string& str, double& result) {
    try {
        std::size_t pos{0};
        result = std::stod(str, &pos);
        if (pos != str.size())
            throw std::invalid_argument{""};
    }
    catch (std::invalid_argument const&) {
        return false;
    }
    catch (std::out_of_range const&) {
        return false;
    }
    return true;
}

// Sets parameter of name key to value val. Returns true on success.
bool read(const std::string& key, const std::string& val, Parameters& param) {
    if (key == "display")
        return parse(val, param.display);
    else if (key == "print_interval")
        return parse(val, param.print_interval);
    else if (key == "time_limit")
        return parse(val, param.time_limit);
    else if (key == "dualize")
        return parse(val, param.dualize);
    else if (key == "scale")
        return parse(val, param.scale);
    else if (key == "ipm_maxiter")
        return parse(val, param.ipm_maxiter);
    else if (key == "ipm_feasibility_tol")
        return parse(val, param.ipm_feasibility_tol);
    else if (key == "ipm_optimality_tol")
        return parse(val, param.ipm_optimality_tol);
    else if (key == "ipm_drop_primal")
        return parse(val, param.ipm_drop_primal);
    else if (key == "ipm_drop_dual")
        return parse(val, param.ipm_drop_dual);
    else if (key == "kkt_tol")
        return parse(val, param.kkt_tol);
    else if (key == "precond_dense_cols")
        return parse(val, param.precond_dense_cols);
    else if (key == "crash_basis")
        return parse(val, param.crash_basis);
    else if (key == "dependency_tol")
        return parse(val, param.dependency_tol);
    else if (key == "volume_tol")
        return parse(val, param.volume_tol);
    else if (key == "rows_per_slice")
        return parse(val, param.rows_per_slice);
    else if (key == "maxskip_updates")
        return parse(val, param.maxskip_updates);
    else if (key == "lu_kernel")
        return parse(val, param.lu_kernel);
    else if (key == "lu_pivottol")
        return parse(val, param.lu_pivottol);
    else if (key == "crossover")
        return parse(val, param.crossover);
    else if (key == "crossover_start")
        return parse(val, param.crossover_start);
    else if (key == "pfeasibility_tol")
        return parse(val, param.pfeasibility_tol);
    else if (key == "dfeasibility_tol")
        return parse(val, param.dfeasibility_tol);
    else if (key == "debug")
        return parse(val, param.debug);
    else if (key == "switchiter")
        return parse(val, param.switchiter);
    else if (key == "stop_at_switch")
        return parse(val, param.stop_at_switch);
    else if (key == "update_heuristic")
        return parse(val, param.update_heuristic);
    else if (key == "maxpasses")
        return parse(val, param.maxpasses);
    else
        return false;
}

template <typename T>
void write(std::ostream& os, const char* name, T value) {
    std::ostringstream s;
    s.setf(std::ios_base::left);
    s.width(32);
    s << name;
    os << s.str() << value << '\n';
}

} // namespace

Int Control::ReadParameters(const char* filename) {
    if (!filename)
        return -1;
    std::ifstream file{filename};
    if (!file.is_open())
        return -1;
    // Only overwrite parameters_ when reading the whole file succeeds.
    Parameters param = parameters_;
    std::string line;
    for (int linenb = 1; std::getline(file, line); linenb++) {
        // Remove everything after the first '#'.
        std::string::size_type pos = line.find('#');
        if (pos != std::string::npos)
            line.resize(pos);
        if (std::all_of(line.cbegin(), line.cend(), isspace))
            continue;
        std::string key, val;
        bool ok = split_two(line, key, val);
        if (ok) {
            ok = read(key, val, param);
        }
        if (!ok) {
            Log() << "Failed to parse line " << linenb << " of "
                  << filename << ":\n" << line << '\n';
            return -1;
        }
    }
    if (file.bad())
        return -1;
    parameters_ = param;
    return 0;
}

Int Control::WriteParameters(const char* filename) {
    if (!filename)
        return -1;
    std::ofstream file{filename};
    if (!file.is_open())
        return -1;
    // TODO: Use std::to_chars() for formatting double values when
    // this is supported by clang/gcc. Then 1e-6 is printed as such,
    // and not as 9.9999999999999995e-07.
    file.precision(17);
    write(file, "display", parameters_.display);
    write(file, "print_interval", parameters_.print_interval);
    write(file, "time_limit", parameters_.time_limit);
    write(file, "dualize", parameters_.dualize);
    write(file, "scale", parameters_.scale);
    write(file, "ipm_maxiter", parameters_.ipm_maxiter);
    write(file, "ipm_feasibility_tol", parameters_.ipm_feasibility_tol);
    write(file, "ipm_optimality_tol", parameters_.ipm_optimality_tol);
    write(file, "ipm_drop_primal", parameters_.ipm_drop_primal);
    write(file, "ipm_drop_dual", parameters_.ipm_drop_dual);
    write(file, "kkt_tol", parameters_.kkt_tol);
    write(file, "precond_dense_cols", parameters_.precond_dense_cols);
    write(file, "crash_basis", parameters_.crash_basis);
    write(file, "dependency_tol", parameters_.dependency_tol);
    write(file, "volume_tol", parameters_.volume_tol);
    write(file, "rows_per_slice", parameters_.rows_per_slice);
    write(file, "maxskip_updates", parameters_.maxskip_updates);
    write(file, "lu_kernel", parameters_.lu_kernel);
    write(file, "lu_pivottol", parameters_.lu_pivottol);
    write(file, "crossover", parameters_.crossover);
    write(file, "crossover_start", parameters_.crossover_start);
    write(file, "pfeasibility_tol", parameters_.pfeasibility_tol);
    write(file, "dfeasibility_tol", parameters_.dfeasibility_tol);
    write(file, "debug", parameters_.debug);
    write(file, "switchiter", parameters_.switchiter);
    write(file, "stop_at_switch", parameters_.stop_at_switch);
    write(file, "update_heuristic", parameters_.update_heuristic);
    write(file, "maxpasses", parameters_.maxpasses);
    return file.fail() ? -1 : 0;
}

void Control::OpenLogfile() {
    logfile_.close();
    const char* filename = parameters_.logfile;
    if (filename && filename[0])
        logfile_.open(filename, std::ios_base::out | std::ios_base::app);
    MakeStream();
}

void Control::CloseLogfile() {
    logfile_.close();
    MakeStream();
}

void Control::ResetTimer() {
    timer_.Reset();
}

void Control::MakeStream() {
    output_.clear();
    if (parameters_.display)
        output_.add(std::cout);
    if (logfile_.is_open())
        output_.add(logfile_);
}

std::string Format(Int i, int width) {
    std::ostringstream s;
    s.width(width);
    s << i;
    return s.str();
}

std::string Format(const char* c, int width) {
    std::ostringstream s;
    s.width(width);
    s << c;
    return s.str();
}

std::string Format(double d, int width, int prec,
                   std::ios_base::fmtflags floatfield) {
    std::ostringstream s;
    s.precision(prec);
    s.width(width);
    s.setf(floatfield, std::ios_base::floatfield);
    s << d;
    return s.str();
}

}  // namespace ipx
