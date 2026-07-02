/* spt2str_wasm.cpp
 * Wasm-ready core refactored from spt2str (Qt/SQLite desktop version).
 * - Removed: sqlite3 (LoadDataFromDB / SaveResultsToDB), hardcoded main(db, layers)
 * - Added:   CSV string input, JSON string output, Embind binding (analyzeCSV)
 * - Unchanged: ALL calculation logic (physics-based inverse for c / phi).
 *
 * Build (native test):  g++ -O2 -std=c++17 spt2str_wasm.cpp -o spt_test && ./spt_test
 * Build (wasm):         emcc spt2str_wasm.cpp -o spt_core.js --bind \
 *                            -s MODULARIZE=1 -s EXPORT_ES6=1 -O2
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cstdlib>

#ifdef __EMSCRIPTEN__
#include <emscripten/bind.h>
#endif

// ---------------- Constants (unchanged) ----------------
constexpr double pa = 100.0;
constexpr double pi = 3.14159265358979;
constexpr double d2r = 0.01745329251994; // degree -> radian
constexpr double m2k = 1000.0;
constexpr double qc2E = 5.0;
constexpr double gw = 10.0;
constexpr double Di_shoe = 0.0349;
constexpr double Do_shoe = 0.0510;
constexpr double Di_tube = 0.0381;
constexpr double Do_tube = 0.0510;
constexpr double Wh = 0.635;   // hammer weight (kN)
constexpr double Hdrop = 0.76; // drop height
constexpr double eff = 0.60;   // hammer efficiency
constexpr double disp = 0.30;  // sampler displacement

// ---------------- Data structures ----------------
struct SPTData {
    double depth;
    double spt;
    int type; // 1: Cohesion (Clay), 0/other: Friction (Sand)
};

struct ProfilePoint {
    double depth, c, f, g, E;
};

struct AnalysisResult {
    std::vector<ProfilePoint> layers;  // layer summary
    std::vector<ProfilePoint> profile; // 0.1 m slice profile (for plotting)
};

// ---------------- Analyzer ----------------
class SPTAnalyzer {
private:
    std::vector<SPTData> raw_spt_data;
    std::vector<SPTData> sliced_spt_data;
    std::vector<double> c_z, f_z, E_z, g_z;
    std::vector<double> layer_boundaries;

    struct LayerResult {
        double c = 0.0, f = 0.0, g = 0.0, E = 0.0;
        int count = 0;
    };
    std::vector<LayerResult> layer_results;

    struct Stratum {
        double start_depth;
        double end_depth;
        int type;
        int allocated_sublayers;
    };

    double ge = 0.0;

    // --- Input adapter (replaces LoadDataFromDB): keeps the z=0 dummy logic ---
    void SetRawData(const std::vector<SPTData>& data) {
        raw_spt_data.clear();
        if (!data.empty()) {
            SPTData zero_point = data.front();
            zero_point.depth = 0.0;
            raw_spt_data.push_back(zero_point);
        }
        for (const auto& d : data) raw_spt_data.push_back(d);
    }

    double Heaviside(double x) const { return (x > 0.0) ? 1.0 : 0.0; }

    double GetProfileValue(int propType, double z) const {
        double sum = 0.0;
        if (raw_spt_data.empty()) return 0.0;

        for (size_t i = 0; i < raw_spt_data.size() - 1; ++i) {
            double x1 = raw_spt_data[i].depth;
            double x2 = raw_spt_data[i + 1].depth;
            double val = (propType == 1) ? raw_spt_data[i + 1].spt
                                         : (double)raw_spt_data[i + 1].type;
            double h1 = Heaviside(z - x1);
            double h2 = Heaviside(z - x2);
            sum += val * (h1 - h2);
        }
        if (!raw_spt_data.empty()) {
            const auto& last = raw_spt_data.back();
            double x1 = last.depth;
            double x2 = last.depth + 1000.0;
            double val = (propType == 1) ? last.spt : (double)last.type;
            double h1 = Heaviside(z - x1);
            double h2 = Heaviside(z - x2);
            sum += val * (h1 - h2);
        }
        return sum;
    }

    void GenerateSlicedData() {
        if (raw_spt_data.empty()) return;
        double zs = raw_spt_data.front().depth;
        double ze = raw_spt_data.back().depth;
        double z = zs;

        while (z <= ze + 0.1001) {
            SPTData slice;
            slice.depth = z;
            slice.spt = GetProfileValue(1, z);
            slice.type = static_cast<int>(GetProfileValue(2, z));
            sliced_spt_data.push_back(slice);
            z += 0.1;
        }
        if (!sliced_spt_data.empty()) {
            sliced_spt_data.front().spt = raw_spt_data.front().spt;
            sliced_spt_data.front().type = raw_spt_data.front().type;
            sliced_spt_data.back().spt = raw_spt_data.back().spt;
            sliced_spt_data.back().type = raw_spt_data.back().type;
        }
    }

    // --- Physics / Mechanics (unchanged) ---
    double Vside(double z, double cohesion, double friction, double D) {
        double K = 1.0 - std::sin(friction * d2r);
        double Vs_cohesion = pi * D * 0.54 * cohesion;
        double Vs_friction = pi * D * ge * z * K * std::tan(friction * d2r * 2.0 / 3.0);
        return Vs_cohesion + Vs_friction;
    }

    double Vbase(double z, double cohesion, double friction) {
        double Aeff = pi * (std::pow(Do_shoe, 2.0) - std::pow(Di_shoe, 2.0)) / 4.0;
        double f_val = (friction <= 0.0) ? 1.e-6 : friction;
        double tan_f = std::tan(f_val * d2r);
        double sin_f = std::sin(f_val * d2r);

        double Nq = std::pow(std::tan(pi / 4.0 + f_val * d2r / 2.0), 2.0) * std::exp(pi * tan_f);
        double Nc = (Nq - 1.0) / tan_f;
        double Ng = 1.5 * (Nq - 1.0) * tan_f;

        double B = (Do_shoe - Di_shoe) / 2.0;
        double L = pi * (Di_shoe + Do_shoe) / 2.0;
        double sq = 1.0 + B / L * sin_f;
        double sc = 1.0 + (Nq / Nc) * B / L;
        double sg = 1.0 - 0.4 * B / L;

        double zi = 0.3;
        double zo = z + 0.3;
        double ki = (zi / B <= 1.0) ? zi / B : std::atan(zi / B);
        double ko = (zo / B <= 1.0) ? zo / B : std::atan(zo / B);

        double term = 2.0 * tan_f * std::pow(1.0 - sin_f, 2.0);
        double dqi = 1.0 + term * ki;
        double dqo = 1.0 + term * ko;
        double dci = 1.0 + 0.4 * ki;
        double dco = 1.0 + 0.4 * ko;
        double dg = 1.0;

        double qult_qi = ge * zi * Nq * sq * dqi;
        double qult_qo = ge * zo * Nq * sq * dqo;
        double qult_ci = cohesion * Nc * sc * dci;
        double qult_co = cohesion * Nc * sc * dco;
        double qult_g  = 0.5 * B * ge * Ng * sg * dg;

        double qult_q = (qult_qi + qult_qo) / 2.0;
        double qult_c = (qult_ci + qult_co) / 2.0;
        return Aeff * (qult_g + qult_q + qult_c);
    }

    double Q_Cohesion(double z, double cohesion) {
        double Qp  = 0.5 * disp * (Vbase(z + 0.15, cohesion, 0.0) + Vbase(z + 0.45, cohesion, 0.0));
        double Qs1 = 0.5 * disp * (Vside(z + 0.15, cohesion, 0.0, Do_tube) + Vside(z + 0.45, cohesion, 0.0, Do_tube));
        double Qs2 = 0.5 * disp * (Vside(z + 0.15, cohesion, 0.0, Di_tube) + Vside(z + 0.45, cohesion, 0.0, Di_tube));
        return Qp + Qs1 + Qs2;
    }

    double Q_Friction(double z, double friction) {
        double Qp  = 0.5 * disp * (Vbase(z + 0.15, 0.0, friction) + Vbase(z + 0.45, 0.0, friction));
        double Qs1 = 0.5 * disp * (Vside(z + 0.15, 0.0, friction, Do_tube) + Vside(z + 0.45, 0.0, friction, Do_tube));
        double Qs2 = 0.5 * disp * (Vside(z + 0.15, 0.0, friction, Di_tube) + Vside(z + 0.45, 0.0, friction, Di_tube));
        return Qp + Qs1 + Qs2;
    }

    double SolveForParameter(bool isCohesion, double N_value) {
        double targetEnergy = N_value * Wh * Hdrop * eff;
        double z_opt = 20.0;

        auto func = [&](double x) -> double {
            double work = isCohesion ? Q_Cohesion(z_opt, x) : Q_Friction(z_opt, x);
            return work - targetEnergy;
        };

        double a = 1e-6;
        double b = isCohesion ? 400.0 : 60.0;
        double fa = func(a);
        double fb = func(b);
        if (fa * fb > 0.0) return b;

        double r = a;
        int side = 0;
        for (int n = 0; n < 100; ++n) {
            r = (fa * b - fb * a) / (fa - fb);
            if (std::abs(b - a) / std::abs(b + a) < 1e-6) break;
            double fr = func(r);
            if (fr * fb > 0.0) {
                b = r; fb = fr;
                if (side == -1) fa /= 2.0;
                side = -1;
            } else if (fa * fr > 0.0) {
                a = r; fa = fr;
                if (side == 1) fb /= 2.0;
                side = 1;
            } else break;
        }
        return r;
    }

    void CalculateStrengthParameters() {
        for (const auto& slice : sliced_spt_data) {
            double N = slice.spt;
            bool isClay = (slice.type == 1);

            if (isClay) {
                double unit_w = 1.93e-4 * std::pow(N, 3) - 1.72e-2 * std::pow(N, 2) + 4.95e-1 * N + 14.4;
                ge = unit_w - gw;
                g_z.push_back(unit_w);
                E_z.push_back(2.0 * 1220.0 * std::pow(N, 0.62) * 1.4);
                c_z.push_back(SolveForParameter(true, N));
                f_z.push_back(0.0);
            } else {
                double unit_w = 6.43e-8 * std::pow(N, 5) - 1.95e-5 * std::pow(N, 4) + 1.75e-3 * std::pow(N, 3) - 6.57e-2 * std::pow(N, 2) + 1.14 * N + 10.8;
                ge = unit_w - gw;
                g_z.push_back(unit_w);
                E_z.push_back(650.0 * std::pow(N, 0.94) * 1.35);
                c_z.push_back(0.0);
                f_z.push_back(SolveForParameter(false, N));
            }
        }
    }

    void DetermineLayerBoundaries(int n_layers) {
        if (sliced_spt_data.empty()) return;
        double H = sliced_spt_data.back().depth;
        layer_boundaries.clear();

        std::vector<Stratum> strata;
        {
            Stratum current;
            current.start_depth = sliced_spt_data[0].depth;
            current.type = sliced_spt_data[0].type;
            for (size_t i = 1; i < sliced_spt_data.size(); ++i) {
                if (sliced_spt_data[i].type != current.type) {
                    current.end_depth = sliced_spt_data[i].depth;
                    strata.push_back(current);
                    current.start_depth = sliced_spt_data[i].depth;
                    current.type = sliced_spt_data[i].type;
                }
            }
            current.end_depth = sliced_spt_data.back().depth;
            strata.push_back(current);
        }

        int total_allocated = 0;
        double total_thickness = H - sliced_spt_data[0].depth;
        if (total_thickness <= 0) total_thickness = H;

        for (auto& s : strata) {
            double thickness = s.end_depth - s.start_depth;
            int count = static_cast<int>(std::round((thickness / total_thickness) * n_layers));
            if (count < 1) count = 1;
            s.allocated_sublayers = count;
            total_allocated += count;
        }

        while (total_allocated != n_layers) {
            auto it = std::max_element(strata.begin(), strata.end(),
                [](const Stratum& a, const Stratum& b) {
                    return (a.end_depth - a.start_depth) < (b.end_depth - b.start_depth);
                });
            if (total_allocated < n_layers) {
                it->allocated_sublayers++;
                total_allocated++;
            } else {
                if (it->allocated_sublayers > 1) {
                    it->allocated_sublayers--;
                    total_allocated--;
                } else break;
            }
        }

        for (const auto& s : strata) {
            double thickness = s.end_depth - s.start_depth;
            double interval = thickness / s.allocated_sublayers;
            for (int i = 1; i <= s.allocated_sublayers; ++i) {
                double b = s.start_depth + i * interval;
                if (i == s.allocated_sublayers) b = s.end_depth;
                layer_boundaries.push_back(b);
            }
        }
        layer_results.resize(layer_boundaries.size());
    }

    void CalculateLayerParameters() {
        if (sliced_spt_data.empty()) return;

        for (size_t i = 0; i < sliced_spt_data.size(); ++i) {
            double depth = sliced_spt_data[i].depth;
            int layer_idx = -1;
            for (size_t j = 0; j < layer_boundaries.size(); ++j) {
                double prev_bound = (j == 0) ? 0.0 : layer_boundaries[j - 1];
                if (depth > prev_bound && depth <= layer_boundaries[j] + 1e-6) {
                    layer_idx = j;
                    break;
                }
            }
            if (layer_idx != -1 && layer_idx < (int)layer_results.size()) {
                layer_results[layer_idx].count++;
                layer_results[layer_idx].c += c_z[i];
                layer_results[layer_idx].f += f_z[i];
                layer_results[layer_idx].g += g_z[i];
                layer_results[layer_idx].E += E_z[i];
            }
        }
        for (auto& lr : layer_results) {
            if (lr.count > 0) {
                lr.c /= lr.count; lr.f /= lr.count; lr.g /= lr.count; lr.E /= lr.count;
            }
        }
        for (size_t i = 0; i < layer_results.size(); ++i) {
            if (layer_results[i].count == 0) {
                if (i > 0 && i < layer_results.size() - 1) {
                    layer_results[i].c = (layer_results[i - 1].c + layer_results[i + 1].c) / 2.0;
                    layer_results[i].f = (layer_results[i - 1].f + layer_results[i + 1].f) / 2.0;
                    layer_results[i].g = (layer_results[i - 1].g + layer_results[i + 1].g) / 2.0;
                    layer_results[i].E = (layer_results[i - 1].E + layer_results[i + 1].E) / 2.0;
                }
            }
        }
    }

    // --- Output adapter (replaces SaveResultsToDB): build struct instead of DB write ---
    AnalysisResult BuildResult() const {
        AnalysisResult res;
        for (size_t i = 0; i < layer_boundaries.size(); ++i) {
            ProfilePoint p{ layer_boundaries[i], layer_results[i].c, layer_results[i].f,
                            layer_results[i].g, layer_results[i].E };
            res.layers.push_back(p);
        }
        for (size_t i = 0; i < sliced_spt_data.size(); ++i) {
            ProfilePoint p{ sliced_spt_data[i].depth, c_z[i], f_z[i], g_z[i], E_z[i] };
            res.profile.push_back(p);
        }
        return res;
    }

public:
    AnalysisResult Analyze(const std::vector<SPTData>& data, int num_layers) {
        SetRawData(data);
        GenerateSlicedData();
        CalculateStrengthParameters();
        DetermineLayerBoundaries(num_layers);
        CalculateLayerParameters();
        return BuildResult();
    }
};

// ---------------- I/O helpers (web boundary) ----------------
static std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// Non-throwing numeric parse (Emscripten disables exception catching by default,
// so we must not rely on std::stod throwing on the header row).
static bool toNumber(const std::string& s, double& out) {
    std::string t = trim(s);
    if (t.empty()) return false;
    char* end = nullptr;
    out = std::strtod(t.c_str(), &end);
    return end != t.c_str() && *end == '\0';
}

std::vector<SPTData> parseCSV(const std::string& csv) {
    std::vector<SPTData> out;
    std::stringstream ss(csv);
    std::string line;
    while (std::getline(ss, line)) {
        line = trim(line);
        if (line.empty()) continue;
        std::stringstream ls(line);
        std::string a, b, c;
        if (!std::getline(ls, a, ',')) continue;
        if (!std::getline(ls, b, ',')) continue;
        if (!std::getline(ls, c, ',')) continue;
        double da, db, dc;
        if (!toNumber(a, da) || !toNumber(b, db) || !toNumber(c, dc))
            continue; // skip header or malformed rows (no exceptions)
        SPTData d;
        d.depth = da;
        d.spt   = db;
        d.type  = (int)std::lround(dc);
        out.push_back(d);
    }
    return out;
}

std::string resultToJSON(const AnalysisResult& r) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(3);
    auto safe = [](double v) { return std::isfinite(v) ? v : 0.0; };
    auto emit = [&](const std::vector<ProfilePoint>& v) {
        ss << "[";
        for (size_t i = 0; i < v.size(); ++i) {
            const auto& p = v[i];
            if (i) ss << ",";
            ss << "{\"depth\":" << safe(p.depth) << ",\"c\":" << safe(p.c) << ",\"f\":" << safe(p.f)
               << ",\"g\":" << safe(p.g) << ",\"E\":" << safe(p.E) << "}";
        }
        ss << "]";
    };
    ss << "{\"layers\":";
    emit(r.layers);
    ss << ",\"profile\":";
    emit(r.profile);
    ss << "}";
    return ss.str();
}

// Single entry point exposed to JS: CSV in, JSON out.
std::string analyzeCSV(const std::string& csv, int num_layers) {
    std::vector<SPTData> data = parseCSV(csv);
    SPTAnalyzer analyzer;
    AnalysisResult res = analyzer.Analyze(data, num_layers);
    return resultToJSON(res);
}

#ifdef __EMSCRIPTEN__
EMSCRIPTEN_BINDINGS(spt_module) {
    emscripten::function("analyzeCSV", &analyzeCSV);
}
#else
int main() {
    std::string csv =
        "depth,SPT,type\n"
        "1.0,4,1\n"
        "2.0,6,1\n"
        "3.0,9,1\n"
        "4.0,15,0\n"
        "6.0,22,0\n"
        "8.0,31,0\n"
        "10.0,44,0\n";
    std::cout << analyzeCSV(csv, 8) << std::endl;
    return 0;
}
#endif
