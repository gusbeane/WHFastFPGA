#include "kepler.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <cstring>
#include <cstdint>

// max number of bit errors in ULP
#define MAX_GS13 1

void read_input(double x_vec[N_PLANETS], double y_vec[N_PLANETS], double z_vec[N_PLANETS],
                double vx_vec[N_PLANETS], double vy_vec[N_PLANETS], double vz_vec[N_PLANETS],
                double m_vec[N_PLANETS], double *M0, double *dt)
{
    std::ifstream file("golden_drift_inputvectors.csv");

    if (!file.is_open())
    {
        std::cerr << "Error opening file." << std::endl;
    }

    // Read the first line for M0
    std::string line;
    if (!std::getline(file, line))
    {
        std::cerr << "Error reading M0." << std::endl;
    }
    std::istringstream m0_stream(line);
    char comma;
    m0_stream >> *M0 >> comma >> *dt;

    // Read the rest of the lines for x, y, z, vx, vy, vz, m
    for (int i = 0; i < N_PLANETS; ++i)
    {
        if (!std::getline(file, line))
        {
            std::cerr << "Error reading planet data." << std::endl;
        }
        std::istringstream ss(line);
        ss >> x_vec[i] >> comma >> y_vec[i] >> comma >> z_vec[i] >> comma >> vx_vec[i] >> comma >> vy_vec[i] >> comma >> vz_vec[i] >> comma >> m_vec[i];
    }
    file.close();
}

void verify_output(double x_vec[N_PLANETS], double y_vec[N_PLANETS], double z_vec[N_PLANETS],
                double vx_vec[N_PLANETS], double vy_vec[N_PLANETS], double vz_vec[N_PLANETS],
                double m_vec[N_PLANETS])
{
    std::ifstream file("golden_drift_outputvectors.csv");

    if (!file.is_open())
    {
        std::cerr << "Error opening file." << std::endl;
    }

    double x_vec_golden[N_PLANETS], y_vec_golden[N_PLANETS], z_vec_golden[N_PLANETS];
    double vx_vec_golden[N_PLANETS], vy_vec_golden[N_PLANETS], vz_vec_golden[N_PLANETS];
    double m_vec_golden[N_PLANETS];

    // Ignore first line
    std::string line;
    if (!std::getline(file, line))
    {
        std::cerr << "Error reading M0." << std::endl;
    }
    // Read the rest of the lines for x, y, z, vx, vy, vz, m
    for (int i = 0; i < N_PLANETS; ++i)
    {
        char comma;
        if (!std::getline(file, line))
        {
            std::cerr << "Error reading planet data." << std::endl;
        }
        std::istringstream ss(line);
        ss >> x_vec_golden[i] >> comma >> y_vec_golden[i] >> comma >> z_vec_golden[i] >> comma >> vx_vec_golden[i] >> comma >> vy_vec_golden[i] >> comma >> vz_vec_golden[i] >> comma >> m_vec_golden[i];
    }

    // Compare the values
    int fail = 0;
    std::cout << std::fixed << std::setprecision(20);
    for (int i = 0; i < N_PLANETS; ++i)
    {
        if (fabs(x_vec[i] - x_vec_golden[i]) > 1e-10)
        {
            std::cout << "x_vec[" << i << "] mismatch: " << x_vec[i] << " vs " << x_vec_golden[i] << std::endl;
            fail++;
        }
        if (fabs(y_vec[i] - y_vec_golden[i]) > 1e-10)
        {
            std::cout << "y_vec[" << i << "] mismatch: " << y_vec[i] << " vs " << y_vec_golden[i] << std::endl;
            fail++;
        }
        if (fabs(z_vec[i] - z_vec_golden[i]) > 1e-10)
        {
            std::cout << "z_vec[" << i << "] mismatch: " << z_vec[i] << " vs " << z_vec_golden[i] << std::endl;
            fail++;
        }
        if (fabs(vx_vec[i] - vx_vec_golden[i]) > 1e-10)
        {
            std::cout << "vx_vec[" << i << "] mismatch: " << vx_vec[i] << " vs " << vx_vec_golden[i] << std::endl;
            fail++;
        }
        if (fabs(vy_vec[i] - vy_vec_golden[i]) > 1e-10)
        {
            std::cout << "vy_vec[" << i << "] mismatch: " << vy_vec[i] << " vs " << vy_vec_golden[i] << std::endl;
            fail++;
        }
        if (fabs(vz_vec[i] - vz_vec_golden[i]) > 1e-10)
        {
            std::cout << "vz_vec[" << i << "] mismatch: " << vz_vec[i] << " vs " << vz_vec_golden[i] << std::endl;
            fail++;
        }
    }
    file.close();
    if(fail)
    {
        std::cout << "Test failed with " << fail << " mismatches." << std::endl;
        std::exit(fail);
    }
    else
    {
        std::cout << "Test passed!" << std::endl;
    }
}

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
            if (ulp_diff(bits[i], golden_bits[i]) > MAX_GS13) fail_any = true;
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
        if (diff > MAX_GS13) {
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

int main()
{
    // First read in the initial conditions
    double x_vec[N_PLANETS], y_vec[N_PLANETS], z_vec[N_PLANETS];
    double vx_vec[N_PLANETS], vy_vec[N_PLANETS], vz_vec[N_PLANETS];
    double m_vec[N_PLANETS];
    double M0, dt;

    // read_input(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, &M0, &dt);

    // printf("M0: %.20f, dt: %.20f\n", M0, dt);
    // printf("x: %.20f, %.20f, %.20f, %.20f, %.20f, %.20f, %.20f, %.20f\n", x_vec[0], x_vec[1], x_vec[2], x_vec[3], x_vec[4], x_vec[5], x_vec[6], x_vec[7]);

    // kepler_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, M0, dt);

    // verify_output(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec);

    test_stiefel(GS13);
    test_stiefel(GS03);

    test_newton_halley(NEWTON);
    test_newton_halley(HALLEY);

    // Now run kepler step
}