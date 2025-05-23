#include "kepler.h"
#include "jump.h"
#include "interaction.h"
#include "whfast.h"
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

// Helper to convert hex string to long (matches long_to_hex)
static long hex_to_long(const std::string& hexstr) {
    unsigned long ul = 0;
    std::stringstream ss;
    ss << std::hex << hexstr;
    ss >> ul;
    return static_cast<long>(ul);
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
    std::ifstream file("golden/" + fname);
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
    std::ifstream file("golden/" + fname);
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
struct IntegrateGoldenSet {
    double M0, dt, dt_half; // Added dt_half
    bodies_t ss_ics;
    bodies_t ss_kepler;
    bodies_t ss_jump;
    bodies_t ss_interact;
};

real_t read_header_line(std::ifstream& file, const std::string& fname) {
    std::string line;
    if (!std::getline(file, line)) {
        std::cerr << "Error reading dt from " << fname << std::endl;
        std::exit(1);
    }
    size_t eqpos = line.find('=');
    if (eqpos == std::string::npos) {
        std::cerr << "Malformed line in " << line << std::endl;
        std::exit(1);
    }
    std::string hex_value = line.substr(eqpos + 1);
    return hex_to_double(hex_value);
}

long read_header_line_long(std::ifstream& file, const std::string& fname) {
    std::string line;
    if (!std::getline(file, line)) {
        std::cerr << "Error reading dt from " << fname << std::endl;
        std::exit(1);
    }
    size_t eqpos = line.find('=');
    if (eqpos == std::string::npos) {
        std::cerr << "Malformed line in " << line << std::endl;
        std::exit(1);
    }
    std::string hex_value = line.substr(eqpos + 1);
    return hex_to_long(hex_value);
}

struct bodies_t read_particle_block(std::ifstream& file, const std::string& fname) {
    bodies_t result;
    std::string line;
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
        result.x_vec[i]  = hex_to_double(fields[0]);
        result.y_vec[i]  = hex_to_double(fields[1]);
        result.z_vec[i]  = hex_to_double(fields[2]);
        result.vx_vec[i] = hex_to_double(fields[3]);
        result.vy_vec[i] = hex_to_double(fields[4]);
        result.vz_vec[i] = hex_to_double(fields[5]);
        result.m_vec[i]  = hex_to_double(fields[6]);
    }
    return result;
}

IntegrateGoldenSet read_integrate_golden() {
    IntegrateGoldenSet result;
    std::string fname = "golden_integrate_step.csv";
    
    std::ifstream file("golden/" + fname);
    if (!file.is_open()) {
        std::cerr << "Error opening " << fname << std::endl;
        std::exit(1);
    }
    
    std::string line;

    // Skip x,y,z,vx,vy,vz,m header
    if (!std::getline(file, line)) {
        std::cerr << "Error reading header from " << fname << std::endl;
        std::exit(1);
    }
    
    // Read dt line (format: dt=hexnumber)
    result.dt = read_header_line(file, fname);

    // Read dt_half line (format: dt_half=hexnumber)
    result.dt_half = read_header_line(file, fname);
    
    // Read M0 line (format: M0=hexnumber)
    result.M0 = read_header_line(file, fname);

    // Read ics
    result.ss_ics = read_particle_block(file, fname);
    // Read kepler step
    result.ss_kepler = read_particle_block(file, fname);
    // Read jump step
    result.ss_jump = read_particle_block(file, fname);
    // Read interaction step
    result.ss_interact = read_particle_block(file, fname);

    // Compute gr prefactors
    result.ss_ics = set_gr_prefac(result.ss_ics, result.M0);
    result.ss_kepler = set_gr_prefac(result.ss_kepler, result.M0);
    result.ss_jump = set_gr_prefac(result.ss_jump, result.M0);
    result.ss_interact = set_gr_prefac(result.ss_interact, result.M0);

    file.close();
    return result;
}

void test_single_step(struct bodies_t ss, struct bodies_t ss_gold, std::string name)
{
    int fail = 0, total = 0;
    unsigned long max_diff = 0;
    const char *vec_names[7] = {"x_vec", "y_vec", "z_vec", "vx_vec", "vy_vec", "vz_vec", "m_vec"};
    for (int i = 0; i < N_PLANETS; ++i)
    {
        double *out_vecs[7] = {ss.x_vec, ss.y_vec, ss.z_vec, ss.vx_vec, ss.vy_vec, ss.vz_vec, ss.m_vec};
        double *gold_vecs[7] = {ss_gold.x_vec, ss_gold.y_vec, ss_gold.z_vec, ss_gold.vx_vec, ss_gold.vy_vec, ss_gold.vz_vec, ss_gold.m_vec};
        for (int v = 0; v < 7; ++v)
        {
            uint64_t bits, golden_bits;
            std::memcpy(&bits, &out_vecs[v][i], sizeof(double));
            std::memcpy(&golden_bits, &gold_vecs[v][i], sizeof(double));
            auto ulp_diff = [](uint64_t a, uint64_t b) -> uint64_t
            {
                return (a > b) ? (a - b) : (b - a);
            };
            unsigned long diff = ulp_diff(bits, golden_bits);
            max_diff = std::max(max_diff, diff);
            if (diff > MAX_ULP_DIFF)
            {
                std::cout << std::setprecision(20);
                std::cout << "Mismatch for planet " << i << ", " << vec_names[v] << ":\n";
                std::cout << "  computed=" << out_vecs[v][i] << " (0x" << std::hex << bits << std::dec << ")"
                          << ", golden=" << gold_vecs[v][i] << " (0x" << std::hex << golden_bits << std::dec << ")"
                          << " | ulp_diff=" << diff << std::endl;
                fail++;
            }
            total++;
        }
    }
    if (fail)
    {
        std::cout << name << " failed with " << fail << " mismatches out of " << total << "." << std::endl;
        std::cout << "Maximum ULP difference: " << max_diff << std::endl;
        std::exit(fail);
    }
    else
    {
        std::cout << name << " passed! (" << total << " cases, max_diff=" << max_diff << ")" << std::endl;
    }
}

void test_integrate_step() {
    struct IntegrateGoldenSet int_gold = read_integrate_golden();
    
    // Run computation
    bodies_t ss_ics = int_gold.ss_ics;
    bodies_t ss_out = kepler_step(ss_ics, int_gold.M0, int_gold.dt_half);
    bodies_t ss_gold = int_gold.ss_kepler;
    test_single_step(ss_out, ss_gold, "test_integrate_step_kepler");

    ss_out = jump_step(ss_out, int_gold.M0, int_gold.dt_half);
    ss_gold = int_gold.ss_jump;
    test_single_step(ss_out, ss_gold, "test_integrate_step_jump");

    ss_out = interaction_step(ss_out, int_gold.M0, int_gold.dt);
    ss_gold = int_gold.ss_interact;
    test_single_step(ss_out, ss_gold, "test_integrate_step_interact");
}

// Struct for golden data
struct KernelGoldenSet {
    double M0, dt, dt_half; // Added dt_half
    long Nint;
    bodies_t ss_ics;
    bodies_t ss_kernel;
};

KernelGoldenSet read_kernel_golden() {
    KernelGoldenSet result;
    std::string fname = "golden_kernel_step.csv";
    
    std::ifstream file("golden/" + fname);
    if (!file.is_open()) {
        std::cerr << "Error opening " << fname << std::endl;
        std::exit(1);
    }
    
    std::string line;

    // Skip x,y,z,vx,vy,vz,m header
    if (!std::getline(file, line)) {
        std::cerr << "Error reading header from " << fname << std::endl;
        std::exit(1);
    }
    
    // Read dt line (format: dt=hexnumber)
    result.dt = read_header_line(file, fname);

    // Read dt_half line (format: dt_half=hexnumber)
    // result.dt_half = read_header_line(file, fname);
    
    // Read M0 line (format: M0=hexnumber)
    result.M0 = read_header_line(file, fname);

    result.Nint = read_header_line_long(file, fname);
    std::cout << "Nint: " << result.Nint << std::endl;

    // Read ics
    result.ss_ics = read_particle_block(file, fname);
    // Read golden step
    result.ss_kernel = read_particle_block(file, fname);

    file.close();
    return result;
}

void test_kernel_step() {
    struct KernelGoldenSet krnl_gold = read_kernel_golden();
    
    // Run computation
    bodies_t ss_ics = krnl_gold.ss_ics;
    bodies_t ss_out = whfast_kernel(ss_ics, krnl_gold.M0, krnl_gold.dt, krnl_gold.Nint);
    bodies_t ss_gold = krnl_gold.ss_kernel;
    test_single_step(ss_out, ss_gold, "test_kernel");
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

    // ---------------------------
    // Integrate tests
    // ---------------------------
    test_integrate_step();

    // ---------------------------
    // Kernel tests
    // ---------------------------
    test_kernel_step();

}