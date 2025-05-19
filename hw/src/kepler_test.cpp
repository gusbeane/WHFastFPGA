#include "kepler.h"
#include "jump.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <cstring>
#include <cstdint>

// max number of bit errors in ULP
#define MAX_ULP_DIFF 0

// Helper to convert hex string to double
static double hex_to_double(const std::string& hexstr) {
    uint64_t u = 0;
    std::stringstream ss;
    ss << std::hex << hexstr;
    ss >> u;
    double d;
    std::memcpy(&d, &u, sizeof(double));
    return d;
}

// Enum for Stiefel mode
enum StiefelMode { GS03, GS13 };

// Struct for golden data
struct StiefelGolden {
    double beta, X;
    double Gs[4]; // Gs0, Gs1, Gs2, Gs3 (Gs0 unused for GS13)
};

// Read golden vectors for Gs03 or Gs13
std::vector<StiefelGolden> read_stiefel_golden(StiefelMode mode) {
    std::vector<StiefelGolden> data;
    std::string fname = (mode == GS03) ? "golden_stiefel_Gs03.csv" : "golden_stiefel_Gs13.csv";
    std::ifstream file(fname);
    if (!file.is_open()) {
        std::cerr << "Error opening " << fname << std::endl;
        std::exit(1);
    }
    std::string line;
    // Skip header
    if (!std::getline(file, line)) {
        std::cerr << "Error reading header from " << fname << std::endl;
        std::exit(1);
    }
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        StiefelGolden entry;
        std::string beta_hex, X_hex, Gs_hex[4];
        if (!std::getline(ss, beta_hex, ',')) break;
        if (!std::getline(ss, X_hex, ',')) break;
        entry.beta = hex_to_double(beta_hex);
        entry.X = hex_to_double(X_hex);
        int n = (mode == GS03) ? 4 : 3;
        for (int i = 0; i < n; ++i) {
            if (!std::getline(ss, Gs_hex[i], ',')) break;
            entry.Gs[i] = hex_to_double(Gs_hex[i]);
        }
        if (mode == GS13) entry.Gs[3] = 0.0; // unused
        data.push_back(entry);
    }
    file.close();
    return data;
}

// Unified test function (no function pointer, just mode)
void test_stiefel(StiefelMode mode) {
    const char *name = (mode == GS03) ? "test_stiefel_Gs03" : "test_stiefel_Gs13";

    auto golden = read_stiefel_golden(mode);
    int fail = 0, total = 0;
    unsigned long max_diff = 0;
    for (const auto& entry : golden) {
        double Gs[4] = {0,0,0,0};
        if (mode == GS03) {
            stiefel_Gs03(&Gs[0], &Gs[1], &Gs[2], &Gs[3], entry.beta, entry.X);
        } else {
            stiefel_Gs13(&Gs[0], &Gs[1], &Gs[2], entry.beta, entry.X); // Gs1, Gs2, Gs3
        }
        uint64_t bits[4], golden_bits[4];
        for (int i = 0; i < 4; ++i) {
            std::memcpy(&bits[i], &Gs[i], sizeof(double));
            std::memcpy(&golden_bits[i], &entry.Gs[i], sizeof(double));
        }
        auto ulp_diff = [](uint64_t a, uint64_t b) -> uint64_t {
            return (a > b) ? (a - b) : (b - a);
        };
        bool fail_any = false;
        for (int i = 0; i < ((mode == GS03) ? 4 : 3); ++i) {
            if (ulp_diff(bits[i], golden_bits[i]) > MAX_ULP_DIFF) fail_any = true;
            max_diff = std::max(max_diff, ulp_diff(bits[i], golden_bits[i]));
        }
        if (fail_any) {
            std::cout << std::setprecision(17);
            std::cout << "Mismatch at line " << (total+2) << ":\n";
            std::cout << "  beta=" << entry.beta << ", X=" << entry.X << std::endl;
            for (int i = 0; i < ((mode == GS03) ? 4 : 3); ++i) {
                std::cout << "  Gs" << i << ": computed=" << Gs[i] << " (0x" << std::hex << bits[i] << std::dec << ")"
                          << ", golden=" << entry.Gs[i] << " (0x" << std::hex << golden_bits[i] << std::dec << ")"
                          << " | ulp_diff=" << ulp_diff(bits[i], golden_bits[i]) << std::endl;
            }
            fail++;
        }
        total++;
    }
    if (fail) {
        std::cout << name << " failed with " << fail << " mismatches out of " << total << "." << std::endl;
        std::cout << "Maximum ULP difference: " << max_diff << std::endl;
        std::exit(fail);
    } else {
        std::cout << name << " passed! (" << total << " cases, max_diff=" << max_diff << ")" << std::endl;
    }
}

// Enum for Newton/Halley mode
enum NewtonHalleyMode { NEWTON, HALLEY };

// Struct for golden data
struct NewtonHalleyGolden {
    double X, beta, r0, eta0, zeta0, Xout;
};

// Read golden vectors for Newton or Halley
struct NewtonHalleyGoldenSet {
    double dt;
    std::vector<NewtonHalleyGolden> data;
};

NewtonHalleyGoldenSet read_newton_halley_golden(NewtonHalleyMode mode) {
    NewtonHalleyGoldenSet result;
    std::string fname = (mode == NEWTON) ? "golden_newton_step.csv" : "golden_halley_step.csv";
    std::ifstream file(fname);
    if (!file.is_open()) {
        std::cerr << "Error opening " << fname << std::endl;
        std::exit(1);
    }
    std::string line;
    // Skip header
    if (!std::getline(file, line)) {
        std::cerr << "Error reading header from " << fname << std::endl;
        std::exit(1);
    }
    // Read dt line (format: dt=hexnumber)
    if (!std::getline(file, line)) {
        std::cerr << "Error reading dt from " << fname << std::endl;
        std::exit(1);
    }
    size_t eqpos = line.find('=');
    if (eqpos == std::string::npos) {
        std::cerr << "Malformed dt line in " << fname << std::endl;
        std::exit(1);
    }
    std::string dt_hex = line.substr(eqpos + 1);
    result.dt = hex_to_double(dt_hex);
    // Read data lines
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string X_hex, beta_hex, r0_hex, eta0_hex, zeta0_hex, Xout_hex;
        if (!std::getline(ss, X_hex, ',')) break;
        if (!std::getline(ss, beta_hex, ',')) break;
        if (!std::getline(ss, r0_hex, ',')) break;
        if (!std::getline(ss, eta0_hex, ',')) break;
        if (!std::getline(ss, zeta0_hex, ',')) break;
        if (!std::getline(ss, Xout_hex, ',')) break;
        NewtonHalleyGolden entry;
        entry.X = hex_to_double(X_hex);
        entry.beta = hex_to_double(beta_hex);
        entry.r0 = hex_to_double(r0_hex);
        entry.eta0 = hex_to_double(eta0_hex);
        entry.zeta0 = hex_to_double(zeta0_hex);
        entry.Xout = hex_to_double(Xout_hex);
        result.data.push_back(entry);
    }
    file.close();
    return result;
}

// Unified test function for Newton/Halley
void test_newton_halley(NewtonHalleyMode mode) {
    const char *name = (mode == NEWTON) ? "test_newton_step" : "test_halley_step";
    auto golden_set = read_newton_halley_golden(mode);
    double dt = golden_set.dt;
    int fail = 0, total = 0;
    unsigned long max_diff = 0;
    for (const auto& entry : golden_set.data) {
        double Xout;
        if (mode == NEWTON) {
            Xout = newton_step(entry.X, entry.beta, entry.r0, entry.eta0, entry.zeta0, dt);
        } else {
            Xout = halley_step(entry.X, entry.beta, entry.r0, entry.eta0, entry.zeta0, dt);
        }
        uint64_t bits, golden_bits;
        std::memcpy(&bits, &Xout, sizeof(double));
        std::memcpy(&golden_bits, &entry.Xout, sizeof(double));
        auto ulp_diff = [](uint64_t a, uint64_t b) -> uint64_t {
            return (a > b) ? (a - b) : (b - a);
        };
        unsigned long diff = ulp_diff(bits, golden_bits);
        max_diff = std::max(max_diff, diff);
        if (diff > MAX_ULP_DIFF) {
            std::cout << std::setprecision(17);
            std::cout << "Mismatch at line " << (total+3) << ":\n";
            std::cout << "  X=" << entry.X << ", beta=" << entry.beta << ", r0=" << entry.r0 << ", eta0=" << entry.eta0 << ", zeta0=" << entry.zeta0 << ", dt=" << dt << std::endl;
            std::cout << "  Xout: computed=" << Xout << " (0x" << std::hex << bits << std::dec << ")"
                      << ", golden=" << entry.Xout << " (0x" << std::hex << golden_bits << std::dec << ")"
                      << " | ulp_diff=" << diff << std::endl;
            fail++;
        }
        total++;
    }
    if (fail) {
        std::cout << name << " failed with " << fail << " mismatches out of " << total << "." << std::endl;
        std::cout << "Maximum ULP difference: " << max_diff << std::endl;
        std::exit(fail);
    } else {
        std::cout << name << " passed! (" << total << " cases, max_diff=" << max_diff << ")" << std::endl;
    }
}

// Struct for golden data
struct KeplerGoldenSet {
    double M0, dt;
    bodies_t ss_in;
    bodies_t ss_golden;
};

KeplerGoldenSet read_kepler_golden() {
    KeplerGoldenSet result;
    std::string fname = "golden_kepler_step.csv";
    std::ifstream file(fname);
    if (!file.is_open()) {
        std::cerr << "Error opening " << fname << std::endl;
        std::exit(1);
    }
    std::string line;
    // Skip header
    if (!std::getline(file, line)) {
        std::cerr << "Error reading header from " << fname << std::endl;
        std::exit(1);
    }
    
    // Read dt line (format: dt=hexnumber)
    if (!std::getline(file, line)) {
        std::cerr << "Error reading dt from " << fname << std::endl;
        std::exit(1);
    }
    size_t eqpos = line.find('=');
    if (eqpos == std::string::npos) {
        std::cerr << "Malformed dt line in " << fname << std::endl;
        std::exit(1);
    }
    std::string dt_hex = line.substr(eqpos + 1);
    result.dt = hex_to_double(dt_hex);

    // Read M0 line (format: M0=hexnumber)
    if (!std::getline(file, line)) {
        std::cerr << "Error reading M0 from " << fname << std::endl;
        std::exit(1);
    }
    eqpos = line.find('=');
    if (eqpos == std::string::npos) {
        std::cerr << "Malformed M0 line in " << fname << std::endl;
        std::exit(1);
    }
    std::string M0_hex = line.substr(eqpos + 1);
    result.M0 = hex_to_double(M0_hex);

    // Read input data for N_PLANETS
    for (int i = 0; i < N_PLANETS; ++i) {
        if (!std::getline(file, line)) {
            std::cerr << "Error: not enough planet lines in " << fname << std::endl;
            std::exit(1);
        }
        std::stringstream ss(line);
        std::string fields[6];
        for (int j = 0; j < 6; ++j) {
            if (!std::getline(ss, fields[j], ',')) {
                std::cerr << "Malformed planet line in " << fname << " at planet " << i << std::endl;
                std::exit(1);
            }
        }
        // Input state
        result.ss_in.x_vec[i]  = hex_to_double(fields[0]);
        result.ss_in.y_vec[i]  = hex_to_double(fields[1]);
        result.ss_in.z_vec[i]  = hex_to_double(fields[2]);
        result.ss_in.vx_vec[i] = hex_to_double(fields[3]);
        result.ss_in.vy_vec[i] = hex_to_double(fields[4]);
        result.ss_in.vz_vec[i] = hex_to_double(fields[5]);
    }

    // Read output data for N_PLANETS
    for (int i = 0; i < N_PLANETS; ++i) {
        if (!std::getline(file, line)) {
            std::cerr << "Error: not enough planet lines in " << fname << std::endl;
            std::exit(1);
        }
        std::stringstream ss(line);
        std::string fields[6];
        for (int j = 0; j < 6; ++j) {
            if (!std::getline(ss, fields[j], ',')) {
                std::cerr << "Malformed planet line in " << fname << " at planet " << i << std::endl;
                std::exit(1);
            }
        }
        // Input state
        result.ss_golden.x_vec[i]  = hex_to_double(fields[0]);
        result.ss_golden.y_vec[i]  = hex_to_double(fields[1]);
        result.ss_golden.z_vec[i]  = hex_to_double(fields[2]);
        result.ss_golden.vx_vec[i] = hex_to_double(fields[3]);
        result.ss_golden.vy_vec[i] = hex_to_double(fields[4]);
        result.ss_golden.vz_vec[i] = hex_to_double(fields[5]);
    }

    file.close();
    return result;
}

void test_kepler_step() {
    struct KeplerGoldenSet kepler_gold = read_kepler_golden();
    
    // Run computation
    bodies_t ss_in = kepler_gold.ss_in;
    bodies_t ss_out = kepler_step(ss_in, kepler_gold.M0, kepler_gold.dt);
    bodies_t ss_golden = kepler_gold.ss_golden;

    int fail = 0, total = 0;
    unsigned long max_diff = 0;
    const char* vec_names[6] = {"x_vec", "y_vec", "z_vec", "vx_vec", "vy_vec", "vz_vec"};
    for (int i = 0; i < N_PLANETS; ++i) {
        double* out_vecs[6] = {ss_out.x_vec, ss_out.y_vec, ss_out.z_vec, ss_out.vx_vec, ss_out.vy_vec, ss_out.vz_vec};
        double* gold_vecs[6] = {ss_golden.x_vec, ss_golden.y_vec, ss_golden.z_vec, ss_golden.vx_vec, ss_golden.vy_vec, ss_golden.vz_vec};
        for (int v = 0; v < 6; ++v) {
            uint64_t bits, golden_bits;
            std::memcpy(&bits, &out_vecs[v][i], sizeof(double));
            std::memcpy(&golden_bits, &gold_vecs[v][i], sizeof(double));
            auto ulp_diff = [](uint64_t a, uint64_t b) -> uint64_t {
                return (a > b) ? (a - b) : (b - a);
            };
            unsigned long diff = ulp_diff(bits, golden_bits);
            max_diff = std::max(max_diff, diff);
            if (diff > MAX_ULP_DIFF) {
                std::cout << std::setprecision(17);
                std::cout << "Mismatch for planet " << i << ", " << vec_names[v] << ":\n";
                std::cout << "  computed=" << out_vecs[v][i] << " (0x" << std::hex << bits << std::dec << ")"
                          << ", golden=" << gold_vecs[v][i] << " (0x" << std::hex << golden_bits << std::dec << ")"
                          << " | ulp_diff=" << diff << std::endl;
                fail++;
            }
            total++;
        }
    }
    if (fail) {
        std::cout << "test_kepler_step failed with " << fail << " mismatches out of " << total << "." << std::endl;
        std::cout << "Maximum ULP difference: " << max_diff << std::endl;
        std::exit(fail);
    } else {
        std::cout << "test_kepler_step passed! (" << total << " cases, max_diff=" << max_diff << ")" << std::endl;
    }
}

// Struct for golden data
struct JumpGoldenSet {
    double M0, dt;
    bodies_t ss_in;
    bodies_t ss_golden;
};

JumpGoldenSet read_jump_golden() {
    JumpGoldenSet result;
    std::string fname = "golden_jump_step.csv";
    std::ifstream file(fname);
    if (!file.is_open()) {
        std::cerr << "Error opening " << fname << std::endl;
        std::exit(1);
    }
    std::string line;
    // Skip header
    if (!std::getline(file, line)) {
        std::cerr << "Error reading header from " << fname << std::endl;
        std::exit(1);
    }
    
    // Read dt line (format: dt=hexnumber)
    if (!std::getline(file, line)) {
        std::cerr << "Error reading dt from " << fname << std::endl;
        std::exit(1);
    }
    size_t eqpos = line.find('=');
    if (eqpos == std::string::npos) {
        std::cerr << "Malformed dt line in " << fname << std::endl;
        std::exit(1);
    }
    std::string dt_hex = line.substr(eqpos + 1);
    result.dt = hex_to_double(dt_hex);

    // Read M0 line (format: M0=hexnumber)
    if (!std::getline(file, line)) {
        std::cerr << "Error reading M0 from " << fname << std::endl;
        std::exit(1);
    }
    eqpos = line.find('=');
    if (eqpos == std::string::npos) {
        std::cerr << "Malformed M0 line in " << fname << std::endl;
        std::exit(1);
    }
    std::string M0_hex = line.substr(eqpos + 1);
    result.M0 = hex_to_double(M0_hex);

    // Read input data for N_PLANETS
    for (int i = 0; i < N_PLANETS; ++i) {
        if (!std::getline(file, line)) {
            std::cerr << "Error: not enough planet lines in " << fname << std::endl;
            std::exit(1);
        }
        std::stringstream ss(line);
        std::string fields[7];
        for (int j = 0; j < 7; ++j) {
            if (!std::getline(ss, fields[j], ',')) {
                std::cerr << "Malformed planet line in " << fname << " at planet " << i << std::endl;
                std::exit(1);
            }
        }
        // Input state
        result.ss_in.x_vec[i]  = hex_to_double(fields[0]);
        result.ss_in.y_vec[i]  = hex_to_double(fields[1]);
        result.ss_in.z_vec[i]  = hex_to_double(fields[2]);
        result.ss_in.vx_vec[i] = hex_to_double(fields[3]);
        result.ss_in.vy_vec[i] = hex_to_double(fields[4]);
        result.ss_in.vz_vec[i] = hex_to_double(fields[5]);
        result.ss_in.m_vec[i] = hex_to_double(fields[6]);
    }

    // Read output data for N_PLANETS
    for (int i = 0; i < N_PLANETS; ++i) {
        if (!std::getline(file, line)) {
            std::cerr << "Error: not enough planet lines in " << fname << std::endl;
            std::exit(1);
        }
        std::stringstream ss(line);
        std::string fields[7];
        for (int j = 0; j < 7; ++j) {
            if (!std::getline(ss, fields[j], ',')) {
                std::cerr << "Malformed planet line in " << fname << " at planet " << i << std::endl;
                std::exit(1);
            }
        }
        // Input state
        result.ss_golden.x_vec[i]  = hex_to_double(fields[0]);
        result.ss_golden.y_vec[i]  = hex_to_double(fields[1]);
        result.ss_golden.z_vec[i]  = hex_to_double(fields[2]);
        result.ss_golden.vx_vec[i] = hex_to_double(fields[3]);
        result.ss_golden.vy_vec[i] = hex_to_double(fields[4]);
        result.ss_golden.vz_vec[i] = hex_to_double(fields[5]);
        result.ss_golden.m_vec[i] = hex_to_double(fields[6]);
    }

    file.close();
    return result;
}

void test_jump_step() {
    struct JumpGoldenSet jump_gold = read_jump_golden();
    
    // Run computation
    bodies_t ss_in = jump_gold.ss_in;
    bodies_t ss_out = jump_step(ss_in, jump_gold.M0, jump_gold.dt);
    bodies_t ss_golden = jump_gold.ss_golden;

    int fail = 0, total = 0;
    unsigned long max_diff = 0;
    const char* vec_names[6] = {"x_vec", "y_vec", "z_vec", "vx_vec", "vy_vec", "vz_vec"};
    for (int i = 0; i < N_PLANETS; ++i) {
        double* out_vecs[6] = {ss_out.x_vec, ss_out.y_vec, ss_out.z_vec, ss_out.vx_vec, ss_out.vy_vec, ss_out.vz_vec};
        double* gold_vecs[6] = {ss_golden.x_vec, ss_golden.y_vec, ss_golden.z_vec, ss_golden.vx_vec, ss_golden.vy_vec, ss_golden.vz_vec};
        for (int v = 0; v < 6; ++v) {
            uint64_t bits, golden_bits;
            std::memcpy(&bits, &out_vecs[v][i], sizeof(double));
            std::memcpy(&golden_bits, &gold_vecs[v][i], sizeof(double));
            auto ulp_diff = [](uint64_t a, uint64_t b) -> uint64_t {
                return (a > b) ? (a - b) : (b - a);
            };
            unsigned long diff = ulp_diff(bits, golden_bits);
            max_diff = std::max(max_diff, diff);
            if (diff > MAX_ULP_DIFF) {
                std::cout << std::setprecision(17);
                std::cout << "Mismatch for planet " << i << ", " << vec_names[v] << ":\n";
                std::cout << "  computed=" << out_vecs[v][i] << " (0x" << std::hex << bits << std::dec << ")"
                          << ", golden=" << gold_vecs[v][i] << " (0x" << std::hex << golden_bits << std::dec << ")"
                          << " | ulp_diff=" << diff << std::endl;
                fail++;
            }
            total++;
        }
    }
    if (fail) {
        std::cout << "test_jump_step failed with " << fail << " mismatches out of " << total << "." << std::endl;
        std::cout << "Maximum ULP difference: " << max_diff << std::endl;
        std::exit(fail);
    } else {
        std::cout << "test_jump_step passed! (" << total << " cases, max_diff=" << max_diff << ")" << std::endl;
    }
}

int main()
{
    // ---------------------------
    // Kepler step tests
    // ---------------------------

    // 1. Stiefel function tests
    test_stiefel(GS13);
    test_stiefel(GS03);

    // 2. Newton and Halley step tests
    test_newton_halley(NEWTON);
    test_newton_halley(HALLEY);

    // 3. Full Kepler step tests
    test_kepler_step();

     // ---------------------------
    // Jump step tests
    // ---------------------------
    test_jump_step();

}