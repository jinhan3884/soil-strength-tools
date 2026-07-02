/* cpt2str_wasm.cpp
 * Wasm-ready core refactored from cpt2str_ver7 (SQLite desktop version).
 * - Removed: sqlite3 (Read_in_situ_Data / Decide_Design_Parameter), DEIntegrator.h (unused),
 *            command-line main().
 * - Added:   CSV string input (depth,qc,fs,u2), JSON string output, Embind (analyzeCSV),
 *            global-state reset so the function is safe to call repeatedly.
 * - Unchanged: ALL calculation logic (Robertson-based, false-position friction solve).
 *
 * Build (native test):  g++ -O2 -std=c++17 cpt2str_wasm.cpp -o cpt_test && ./cpt_test
 * Build (wasm):         emcc cpt2str_wasm.cpp -o cpt_core.js --bind \
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

using namespace std;

// ---------------- Constants (unchanged) ----------------
const double   pi = 3.14159265358979;
const double  d2r = 0.01745329251994; // degree -> radian
const double qc2E = 5.;                // empirical factor
const double   gw = 10.;
const double   pa = 100.;
double  m2k_qc, m2k_fs, m2k_u2;        // unit conversion

// ---------------- Global state (as in original) ----------------
int reading;
int data_count;
int num_of_layer;
int number_of_data_in_layer[8];
double c_layer[8], f_layer[8], E_layer[8], g_layer[8];
double dp_average[8], qc_average[8], fs_average[8], u2_average[8];

class CPT {
public:
    double depth, qc, fs, u2;
    CPT() {};
};
vector<CPT> result;
vector<double> Layer;
vector<double> x;

// ---------------- Output structures (shared with SPT shell) ----------------
struct ProfilePoint { double depth, c, f, g, E; int sbt = 0; };
struct AnalysisResult {
    vector<ProfilePoint> layers;
    vector<ProfilePoint> profile;
};

// Robertson (2010) non-normalized SBT index (I_SBT) -> zone number.
// I_SBT = sqrt[(3.47 - log10(qc/pa))^2 + (log10(Rf) + 1.22)^2], Rf in %, pa = 100 kPa.
// Index applies to zones 2..7; zones 1/8/9 need separate criteria (returned as their
// neighbours here, flagged in UI as index-not-applicable cases).
static int classifySBT(double qc_kPa, double Rf_percent) {
    if (qc_kPa <= 0.0 || Rf_percent <= 0.0) return 0;
    double Isbt = sqrt( pow(3.47 - log10(qc_kPa / pa), 2.0)
                      + pow(log10(Rf_percent) + 1.22, 2.0) );
    if (Isbt > 3.60) return 2; // organic soils
    if (Isbt > 2.95) return 3; // clays: clay to silty clay
    if (Isbt > 2.60) return 4; // silt mixtures
    if (Isbt > 2.05) return 5; // sand mixtures
    if (Isbt > 1.31) return 6; // sands: clean to silty sand
    return 7;                  // gravelly sand to dense sand
}

// prototypes
double Falsi_Method(double (*)(double), double, double);
double Q(double);
void Find_Strength_Parameter();
double Find_Friction_Angle();
void Calculate_Average();
void Auto_Layer_Division();

// ---------------- Physics (unchanged) ----------------
double Q(double x)
{
    if (x <= 0.) x = 1.e-3;
    double Rf = fs_average[reading] * m2k_fs / (qc_average[reading] * m2k_qc); if (Rf <= 0.) Rf = 1e-3;
    double gt = (0.27 * log10(Rf * 100.) + 0.36 * log10(qc_average[reading] * m2k_qc / pa) + 1.236) * gw; // Robertson
    double sv_t = dp_average[reading] * gt;
    double sv_e = dp_average[reading] * (gt - gw); if (sv_e < 0.) sv_e = 1.0e-8;
    double k0 = 1. - sin(x * d2r);
    double sh_t = k0 * sv_t;
    double sh_e = k0 * sv_e;
    double smean_t = (sv_t + 2. * sh_t) / 3.;
    double smean_e = (sv_e + 2. * sh_e) / 3.;
    double qt = qc_average[reading] * m2k_qc;
    double shc_e = 7.89e-4 * pow(((qt - smean_t) / smean_e), 1.44) * sh_e; if (qt <= smean_t) shc_e = sh_e;
    double Nq = pow(tan(pi / 4. + x * d2r / 2.), 2.) * exp(pi * tan(x * d2r));
    double Nc = (Nq - 1.) / tan(x * d2r);
    double Ng = 2. * (Nq + 1.) * tan(x * d2r);
    double Nu = 6. * tan(x * d2r) * (1. + tan(x * d2r));
    double c = fs_average[reading] * m2k_fs - shc_e * tan(2. * x * d2r / 3.); if (c < 0.) c = 0.;
    double q = (gt - gw) * dp_average[reading];
    double f = qt - Nu * u2_average[reading] * m2k_u2 - c * Nc - q * Nq - 0.5 * gt * 0.0357 * Ng;
    return f;
}

void Find_Strength_Parameter()
{
    double qc_total = 0., fs_total = 0., u2_total = 0.;
    for (int i = 0; i < data_count; i++) {
        qc_total += result[i].qc;
        fs_total += result[i].fs;
        u2_total += result[i].u2;
    }

    double qc_avg = qc_total / data_count;
    double fs_avg = fs_total / data_count;
    double u2_avg = u2_total / data_count;

    double RF_avg = fs_avg / qc_avg * 100.;
    double Bq_avg = u2_avg / qc_avg * 100.;

    // unit conversion consideration
    m2k_qc = 1000.;
    if (RF_avg > 8.0) m2k_fs = 1.; else m2k_fs = 1000.;
    if (Bq_avg > 1.4) m2k_u2 = 1.; else m2k_u2 = 1000.;

    for (int i = 0; i < num_of_layer; i++) {
        reading = i;
        double friction = Find_Friction_Angle();
        if (friction < 0.) friction = 0.001;
        double Rf = fs_average[reading] * m2k_fs / (qc_average[reading] * m2k_qc); if (Rf <= 0.) Rf = 1.e-3;
        double gt = (0.27 * log10(Rf * 100.) + 0.36 * log10(qc_average[reading] * m2k_qc / pa) + 1.236) * gw;
        if (gt < 0.) gt = 1e-8;
        double sv_t = dp_average[reading] * gt;
        double sv_e = dp_average[reading] * (gt - gw); if (sv_e < 0.) sv_e = 1.0e-8;
        double k0 = 1. - sin(friction * d2r);
        double sh_t = k0 * sv_t;
        double sh_e = k0 * sv_e;
        double smean_t = (sv_t + 2. * sh_t) / 3.;
        double smean_e = (sv_e + 2. * sh_e) / 3.;
        double qt = qc_average[reading] * m2k_qc;
        double shc_e = 7.89e-4 * pow(((qt - smean_t) / smean_e), 1.44) * sh_e; if (qt <= smean_t) shc_e = sh_e;
        double cohesion = fs_average[i] * m2k_fs - shc_e * tan(friction * d2r * 2. / 3.);
        if (cohesion < 0.) cohesion = 1e-8;
        if (qc_average[reading] == 1e-8 && fs_average[reading] == 1e-8) { cohesion = 1e-8; friction = 1e-8; }
        c_layer[i] = cohesion;
        f_layer[i] = friction;
        g_layer[i] = gt;
    }
}

double Find_Friction_Angle()
{
    double a = 1.e-3;
    double b = 60.;
    double f = Falsi_Method(&Q, a, b);
    return f;
}

double Falsi_Method(double (*func)(double), double a, double b)
{
    double r = a, fr;
    int n, side = 0;
    int m = 100;
    double e = 1.e-6;
    double fa = func(a);
    double fb = func(b);
    for (n = 0; n < m; n++) {
        r = (fa * b - fb * a) / (fa - fb);
        if (fabs(b - a) / fabs(b + a) < e) break;
        fr = func(r);
        if (fr * fb > 0.) {
            b = r; fb = fr;
            if (side == -1) fa /= 2.;
            side = -1;
        } else if (fa * fr > 0.) {
            a = r; fa = fr;
            if (side == +1) fb /= 2.;
            side = +1;
        } else break;
    }
    return r;
}

void Calculate_Average()
{
    for (int i = 0; i < data_count; i++) {
        for (int j = 0; j < num_of_layer; j++) {
            if (result[i].depth > Layer[j] && result[i].depth <= Layer[j + 1]) {
                int layer = j;
                number_of_data_in_layer[layer]++;
                dp_average[layer] += result[i].depth;
                qc_average[layer] += result[i].qc;
                fs_average[layer] += result[i].fs;
                u2_average[layer] += result[i].u2;
            }
        }
    }

    for (int j = 0; j < num_of_layer; j++) {
        if (number_of_data_in_layer[j] == 0) { // guard empty layer (would be nan otherwise)
            dp_average[j] = 1.0e-8; qc_average[j] = 1.0e-8;
            fs_average[j] = 1.0e-8; u2_average[j] = 1.0e-8;
            continue;
        }
        dp_average[j] /= number_of_data_in_layer[j]; if (dp_average[j] <= 0.) dp_average[j] = 1.0e-8;
        qc_average[j] /= number_of_data_in_layer[j]; if (qc_average[j] <= 0.) qc_average[j] = 1.0e-8;
        fs_average[j] /= number_of_data_in_layer[j]; if (fs_average[j] <= 0.) fs_average[j] = 1.0e-8;
        u2_average[j] /= number_of_data_in_layer[j]; if (u2_average[j] <= 0.) u2_average[j] = 1.0e-8;
    }
}

void Auto_Layer_Division()
{
    num_of_layer = 8;
    double H = result[data_count - 1].depth;
    Layer.push_back(0.0 * H);
    Layer.push_back(0.1 * H);
    Layer.push_back(0.2 * H);
    Layer.push_back(0.3 * H);
    Layer.push_back(0.4 * H);
    Layer.push_back(0.5 * H);
    Layer.push_back(0.6 * H);
    Layer.push_back(0.8 * H);
    Layer.push_back(1.0 * H);
}

// ---------------- Web boundary ----------------
static void ResetState() {
    reading = 0; data_count = 0; num_of_layer = 0;
    m2k_qc = m2k_fs = m2k_u2 = 0.;
    for (int i = 0; i < 8; i++) {
        number_of_data_in_layer[i] = 0;
        c_layer[i] = f_layer[i] = E_layer[i] = g_layer[i] = 0.;
        dp_average[i] = qc_average[i] = fs_average[i] = u2_average[i] = 0.;
    }
    result.clear();
    Layer.clear();
    x.clear();
}

static string trim(const string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

static bool toNumber(const string& s, double& out) {
    string t = trim(s);
    if (t.empty()) return false;
    char* end = nullptr;
    out = strtod(t.c_str(), &end);
    return end != t.c_str() && *end == '\0';
}

static void parseCSV(const string& csv) {
    stringstream ss(csv);
    string line;
    while (getline(ss, line)) {
        line = trim(line);
        if (line.empty()) continue;
        stringstream ls(line);
        string a, b, c, d;
        if (!getline(ls, a, ',')) continue;
        if (!getline(ls, b, ',')) continue;
        if (!getline(ls, c, ',')) continue;
        if (!getline(ls, d, ',')) continue;
        double da, db, dc, dd;
        if (!toNumber(a, da) || !toNumber(b, db) || !toNumber(c, dc) || !toNumber(d, dd))
            continue; // skip header / malformed
        CPT p;
        p.depth = fabs(da);
        p.qc    = fabs(db);
        p.fs    = fabs(dc);
        p.u2    = dd; // u2 may be negative
        result.push_back(p);
    }
    data_count = (int)result.size();
}

static AnalysisResult BuildResult() {
    AnalysisResult res;
    for (int i = 0; i < num_of_layer; i++) {
        double E = qc_average[i] * m2k_qc * qc2E;
        double qc_kPa = qc_average[i] * m2k_qc;
        double Rf_pct = (qc_average[i] > 0.0)
            ? (fs_average[i] * m2k_fs) / (qc_average[i] * m2k_qc) * 100.0 : 0.0;
        int sbt = classifySBT(qc_kPa, Rf_pct);
        ProfilePoint lp{ Layer[i + 1], c_layer[i], f_layer[i], g_layer[i], E, sbt };
        res.layers.push_back(lp);
        // stepped profile: top & bottom of each layer carry the layer value
        res.profile.push_back(ProfilePoint{ Layer[i],     c_layer[i], f_layer[i], g_layer[i], E, sbt });
        res.profile.push_back(ProfilePoint{ Layer[i + 1], c_layer[i], f_layer[i], g_layer[i], E, sbt });
    }
    return res;
}

static string resultToJSON(const AnalysisResult& r) {
    stringstream ss;
    ss << fixed << setprecision(3);
    auto safe = [](double v) { return isfinite(v) ? v : 0.0; };
    auto emit = [&](const vector<ProfilePoint>& v) {
        ss << "[";
        for (size_t i = 0; i < v.size(); ++i) {
            const auto& p = v[i];
            if (i) ss << ",";
            ss << "{\"depth\":" << safe(p.depth) << ",\"c\":" << safe(p.c) << ",\"f\":" << safe(p.f)
               << ",\"g\":" << safe(p.g) << ",\"E\":" << safe(p.E) << ",\"sbt\":" << p.sbt << "}";
        }
        ss << "]";
    };
    ss << "{\"layers\":"; emit(r.layers);
    ss << ",\"profile\":"; emit(r.profile);
    ss << "}";
    return ss.str();
}

// Single entry point exposed to JS. num_layers is accepted for API parity with the
// SPT engine, but CPT uses the fixed 8-layer (0.1H x6 + 0.2H x2) scheme by design.
string analyzeCSV(const string& csv, int num_layers) {
    (void)num_layers;
    ResetState();
    parseCSV(csv);
    if (data_count < 2) return "{\"layers\":[],\"profile\":[]}";
    Auto_Layer_Division();
    Calculate_Average();
    Find_Strength_Parameter();
    return resultToJSON(BuildResult());
}

#ifdef __EMSCRIPTEN__
EMSCRIPTEN_BINDINGS(cpt_module) {
    emscripten::function("analyzeCSV", &analyzeCSV);
}
#else
int main() {
    string csv =
        "depth,qc,fs,u2\n"
        "1.0,2.1,0.025,0.10\n"
        "2.0,3.4,0.035,0.15\n"
        "3.0,5.2,0.050,0.20\n"
        "4.0,6.8,0.061,0.22\n"
        "5.0,8.1,0.070,0.25\n"
        "6.0,10.5,0.082,0.28\n"
        "7.0,12.0,0.090,0.30\n"
        "8.0,14.3,0.101,0.33\n"
        "9.0,16.0,0.110,0.35\n"
        "10.0,18.5,0.120,0.38\n";
    cout << analyzeCSV(csv, 8) << endl;
    return 0;
}
#endif
