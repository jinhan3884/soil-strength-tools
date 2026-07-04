/* cpt2str_wasm.cpp  (v2 — slice-based, aligned with the SPT engine)
 * Reconstructed so CPT works like SPT:
 *   raw CSV (depth,qc,fs,u2, possibly 1000s of rows)
 *     -> resample to 0.1 m slices by LINEAR INTERPOLATION
 *     -> compute strength parameters (c, phi, gamma, E, SBT) at EACH slice
 *     -> profile  = per-slice parameters   (drawn as the graph)
 *     -> layers   = slices grouped into num_layers EQUAL layers, averaged (the table)
 * Calculation physics (Motaghedi & Eslami Q() + false-position solve) is unchanged;
 * only the unit it operates on changed from "8 fixed layers" to "0.1 m slices".
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

// ---------------- Constants ----------------
const double   pi = 3.14159265358979;
const double  d2r = 0.01745329251994; // degree -> radian
const double qc2E = 5.;                // empirical factor
const double   gw = 10.;
const double   pa = 100.;
double  m2k_qc, m2k_fs, m2k_u2;        // unit conversion

// ---------------- State ----------------
class CPT { public: double depth, qc, fs, u2; CPT() {} };
vector<CPT> raw;                       // parsed raw data

// current slice values that Q() operates on
double cur_qc, cur_fs, cur_u2, cur_dp;

// 0.1 m sliced data + per-slice results
vector<double> s_depth, s_qc, s_fs, s_u2;
vector<double> s_c, s_f, s_g, s_E;
vector<int>    s_sbt;

// uniform layer boundaries
vector<double> Layer;
int num_of_layer;

// ---------------- Output structures (shared with SPT shell) ----------------
struct ProfilePoint { double depth, c, f, g, E; int sbt = 0; };
struct AnalysisResult {
    vector<ProfilePoint> layers;
    vector<ProfilePoint> profile;
};

// Robertson (2010) non-normalized SBT index -> zone number (applies to zones 2..7).
static int classifySBT(double qc_kPa, double Rf_percent) {
    if (qc_kPa <= 0.0 || Rf_percent <= 0.0) return 0;
    double Isbt = sqrt( pow(3.47 - log10(qc_kPa / pa), 2.0)
                      + pow(log10(Rf_percent) + 1.22, 2.0) );
    if (Isbt > 3.60) return 2;
    if (Isbt > 2.95) return 3;
    if (Isbt > 2.60) return 4;
    if (Isbt > 2.05) return 5;
    if (Isbt > 1.31) return 6;
    return 7;
}

// ---------------- Physics (unchanged; now reads current slice, not layer average) ----------------
double Q(double x)
{
    if (x <= 0.) x = 1.e-3;
    double Rf = cur_fs * m2k_fs / (cur_qc * m2k_qc); if (Rf <= 0.) Rf = 1e-3;
    double gt = (0.27 * log10(Rf * 100.) + 0.36 * log10(cur_qc * m2k_qc / pa) + 1.236) * gw; // Robertson
    double sv_t = cur_dp * gt;
    double sv_e = cur_dp * (gt - gw); if (sv_e < 0.) sv_e = 1.0e-8;
    double k0 = 1. - sin(x * d2r);
    double sh_t = k0 * sv_t;
    double sh_e = k0 * sv_e;
    double smean_t = (sv_t + 2. * sh_t) / 3.;
    double smean_e = (sv_e + 2. * sh_e) / 3.;
    double qt = cur_qc * m2k_qc;
    double shc_e = 7.89e-4 * pow(((qt - smean_t) / smean_e), 1.44) * sh_e; if (qt <= smean_t) shc_e = sh_e;
    double Nq = pow(tan(pi / 4. + x * d2r / 2.), 2.) * exp(pi * tan(x * d2r));
    double Nc = (Nq - 1.) / tan(x * d2r);
    double Ng = 2. * (Nq + 1.) * tan(x * d2r);
    double Nu = 6. * tan(x * d2r) * (1. + tan(x * d2r));
    double c = cur_fs * m2k_fs - shc_e * tan(2. * x * d2r / 3.); if (c < 0.) c = 0.;
    double q = (gt - gw) * cur_dp;
    double f = qt - Nu * cur_u2 * m2k_u2 - c * Nc - q * Nq - 0.5 * gt * 0.0357 * Ng;
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
        if (fr * fb > 0.) { b = r; fb = fr; if (side == -1) fa /= 2.; side = -1; }
        else if (fa * fr > 0.) { a = r; fa = fr; if (side == +1) fb /= 2.; side = +1; }
        else break;
    }
    return r;
}

double Find_Friction_Angle() { return Falsi_Method(&Q, 1.e-3, 60.); }

// ---------------- Determine unit conversion (from raw-data averages, as original) ----------------
static void DetermineUnits()
{
    double qc_total = 0., fs_total = 0., u2_total = 0.;
    for (const auto& p : raw) { qc_total += p.qc; fs_total += p.fs; u2_total += p.u2; }
    int n = (int)raw.size(); if (n == 0) n = 1;
    double qc_avg = qc_total / n; if (qc_avg <= 0.) qc_avg = 1e-8;
    double RF_avg = fs_total / n / qc_avg * 100.;
    double Bq_avg = u2_total / n / qc_avg * 100.;
    m2k_qc = 1000.;
    m2k_fs = (RF_avg > 8.0) ? 1. : 1000.;
    m2k_u2 = (Bq_avg > 1.4) ? 1. : 1000.;
}

// ---------------- Resample raw data to 0.1 m slices via linear interpolation ----------------
static void GenerateSlicedData()
{
    if (raw.size() < 2) return;
    sort(raw.begin(), raw.end(), [](const CPT& a, const CPT& b) { return a.depth < b.depth; });
    double z0 = raw.front().depth;
    double zN = raw.back().depth;
    if (zN - z0 < 1e-9) return;

    size_t idx = 0;
    for (double z = z0; z <= zN + 1e-9; z += 0.1) {
        while (idx + 1 < raw.size() && raw[idx + 1].depth <= z) idx++;
        size_t i0 = idx;
        size_t i1 = (idx + 1 < raw.size()) ? idx + 1 : idx;
        double d0 = raw[i0].depth, d1 = raw[i1].depth;
        double t = (d1 > d0) ? (z - d0) / (d1 - d0) : 0.0;
        if (t < 0.) t = 0.; if (t > 1.) t = 1.;
        s_depth.push_back(z);
        s_qc.push_back(raw[i0].qc + t * (raw[i1].qc - raw[i0].qc));
        s_fs.push_back(raw[i0].fs + t * (raw[i1].fs - raw[i0].fs));
        s_u2.push_back(raw[i0].u2 + t * (raw[i1].u2 - raw[i0].u2));
    }
}

// ---------------- Compute strength parameters at each slice ----------------
static void ComputeSliceStrength()
{
    for (size_t i = 0; i < s_depth.size(); ++i) {
        cur_qc = (s_qc[i] > 0.) ? s_qc[i] : 1e-8;
        cur_fs = (s_fs[i] > 0.) ? s_fs[i] : 1e-8;
        cur_u2 = s_u2[i];
        cur_dp = (s_depth[i] > 0.) ? s_depth[i] : 1e-8;

        double friction = Find_Friction_Angle();
        if (friction < 0.) friction = 0.001;

        double Rf = cur_fs * m2k_fs / (cur_qc * m2k_qc); if (Rf <= 0.) Rf = 1.e-3;
        double gt = (0.27 * log10(Rf * 100.) + 0.36 * log10(cur_qc * m2k_qc / pa) + 1.236) * gw;
        if (gt < 0.) gt = 1e-8;
        double sv_t = cur_dp * gt;
        double sv_e = cur_dp * (gt - gw); if (sv_e < 0.) sv_e = 1.0e-8;
        double k0 = 1. - sin(friction * d2r);
        double sh_t = k0 * sv_t;
        double sh_e = k0 * sv_e;
        double smean_t = (sv_t + 2. * sh_t) / 3.;
        double smean_e = (sv_e + 2. * sh_e) / 3.;
        double qt = cur_qc * m2k_qc;
        double shc_e = 7.89e-4 * pow(((qt - smean_t) / smean_e), 1.44) * sh_e; if (qt <= smean_t) shc_e = sh_e;
        double cohesion = cur_fs * m2k_fs - shc_e * tan(friction * d2r * 2. / 3.);
        if (cohesion < 0.) cohesion = 1e-8;

        double E = cur_qc * m2k_qc * qc2E;
        double Rf_pct = cur_fs * m2k_fs / (cur_qc * m2k_qc) * 100.;
        int sbt = classifySBT(cur_qc * m2k_qc, Rf_pct);

        s_c.push_back(cohesion);
        s_f.push_back(friction);
        s_g.push_back(gt);
        s_E.push_back(E);
        s_sbt.push_back(sbt);
    }
}

// ---------------- Uniform layer division over the slice range ----------------
static void DivideLayers(int n)
{
    if (n < 1) n = 1;
    num_of_layer = n;
    double z0 = s_depth.front();
    double zN = s_depth.back();
    Layer.clear();
    for (int i = 0; i <= n; ++i) Layer.push_back(z0 + (zN - z0) * i / n);
}

// ---------------- Build result: profile (slices) + layers (averaged) ----------------
static AnalysisResult BuildResult()
{
    AnalysisResult res;

    // profile = every slice
    for (size_t i = 0; i < s_depth.size(); ++i)
        res.profile.push_back(ProfilePoint{ s_depth[i], s_c[i], s_f[i], s_g[i], s_E[i], s_sbt[i] });

    // layers = average of slices within each uniform layer
    vector<double> Lc(num_of_layer, 0), Lf(num_of_layer, 0), Lg(num_of_layer, 0), LE(num_of_layer, 0);
    vector<double> Lqc(num_of_layer, 0), Lfs(num_of_layer, 0);
    vector<int> Lcnt(num_of_layer, 0);

    for (size_t i = 0; i < s_depth.size(); ++i) {
        double d = s_depth[i];
        int j = -1;
        for (int k = 0; k < num_of_layer; ++k) {
            double lo = (k == 0) ? -1e18 : Layer[k];
            if (d > lo && d <= Layer[k + 1] + 1e-9) { j = k; break; }
        }
        if (j < 0) j = num_of_layer - 1;
        Lc[j] += s_c[i]; Lf[j] += s_f[i]; Lg[j] += s_g[i]; LE[j] += s_E[i];
        Lqc[j] += s_qc[i]; Lfs[j] += s_fs[i]; Lcnt[j]++;
    }

    double pc = 0, pf = 0, pg = 0, pE = 0; int psbt = 0; // carry-forward for empty layers
    for (int j = 0; j < num_of_layer; ++j) {
        double c, f, g, E; int sbt;
        if (Lcnt[j] > 0) {
            c = Lc[j] / Lcnt[j]; f = Lf[j] / Lcnt[j]; g = Lg[j] / Lcnt[j]; E = LE[j] / Lcnt[j];
            double qc_avg = Lqc[j] / Lcnt[j], fs_avg = Lfs[j] / Lcnt[j];
            double Rf_pct = (qc_avg > 0.) ? (fs_avg * m2k_fs) / (qc_avg * m2k_qc) * 100. : 0.;
            sbt = classifySBT(qc_avg * m2k_qc, Rf_pct);
            pc = c; pf = f; pg = g; pE = E; psbt = sbt;
        } else { c = pc; f = pf; g = pg; E = pE; sbt = psbt; }
        res.layers.push_back(ProfilePoint{ Layer[j + 1], c, f, g, E, sbt });
    }
    return res;
}

// ---------------- Web boundary ----------------
static void ResetState() {
    m2k_qc = m2k_fs = m2k_u2 = 0.;
    cur_qc = cur_fs = cur_u2 = cur_dp = 0.;
    num_of_layer = 0;
    raw.clear();
    s_depth.clear(); s_qc.clear(); s_fs.clear(); s_u2.clear();
    s_c.clear(); s_f.clear(); s_g.clear(); s_E.clear(); s_sbt.clear();
    Layer.clear();
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
        raw.push_back(p);
    }
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

// Single entry point exposed to JS: CSV in, JSON out. num_layers now honoured.
string analyzeCSV(const string& csv, int num_layers) {
    ResetState();
    parseCSV(csv);
    if (raw.size() < 2) return "{\"layers\":[],\"profile\":[]}";
    DetermineUnits();
    GenerateSlicedData();
    if (s_depth.empty()) return "{\"layers\":[],\"profile\":[]}";
    ComputeSliceStrength();
    DivideLayers(num_layers);
    return resultToJSON(BuildResult());
}

#ifdef __EMSCRIPTEN__
EMSCRIPTEN_BINDINGS(cpt_module) {
    emscripten::function("analyzeCSV", &analyzeCSV);
}
#else
int main() {
    // dense-ish sample to mimic real CPT (0.5 m here for brevity)
    string csv = "depth,qc,fs,u2\n";
    // build a synthetic profile 0.5..10 m
    // (kept short; real data would be 0.01 m spacing)
    const char* rows =
        "0.5,1.8,0.020,0.08\n1.0,2.1,0.025,0.10\n1.5,2.7,0.030,0.12\n2.0,3.4,0.035,0.15\n"
        "2.5,4.3,0.043,0.18\n3.0,5.2,0.050,0.20\n3.5,6.0,0.055,0.21\n4.0,6.8,0.061,0.22\n"
        "4.5,7.5,0.066,0.24\n5.0,8.1,0.070,0.25\n6.0,10.5,0.082,0.28\n7.0,12.0,0.090,0.30\n"
        "8.0,14.3,0.101,0.33\n9.0,16.0,0.110,0.35\n10.0,18.5,0.120,0.38\n";
    csv += rows;
    string out = analyzeCSV(csv, 8);
    // print a compact summary
    cout << out.substr(0, 400) << " ...\n";
    return 0;
}
#endif
