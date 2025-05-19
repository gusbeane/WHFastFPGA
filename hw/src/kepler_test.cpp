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

void test_stiefel_Gs13() {
    std::ifstream file("golden_stiefel_Gs13.csv");
    if (!file.is_open()) {
        std::cerr << "Error opening golden_stiefel_Gs13.csv" << std::endl;
        std::exit(1);
    }
    std::string line;
    // Skip header
    if (!std::getline(file, line)) {
        std::cerr << "Error reading header from golden_stiefel_Gs13.csv" << std::endl;
        std::exit(1);
    }
    int fail = 0, total = 0;
    unsigned long max_diff = 0;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string beta_hex, X_hex, Gs1_hex, Gs2_hex, Gs3_hex;
        if (!std::getline(ss, beta_hex, ',')) break;
        if (!std::getline(ss, X_hex, ',')) break;
        if (!std::getline(ss, Gs1_hex, ',')) break;
        if (!std::getline(ss, Gs2_hex, ',')) break;
        if (!std::getline(ss, Gs3_hex, ',')) break;
        double beta = hex_to_double(beta_hex);
        double X = hex_to_double(X_hex);
        double Gs1_golden = hex_to_double(Gs1_hex);
        double Gs2_golden = hex_to_double(Gs2_hex);
        double Gs3_golden = hex_to_double(Gs3_hex);
        double Gs1, Gs2, Gs3;
        stiefel_Gs13(&Gs1, &Gs2, &Gs3, beta, X);
        uint64_t Gs1_bits, Gs2_bits, Gs3_bits;
        uint64_t Gs1_golden_bits, Gs2_golden_bits, Gs3_golden_bits;
        std::memcpy(&Gs1_bits, &Gs1, sizeof(double));
        std::memcpy(&Gs2_bits, &Gs2, sizeof(double));
        std::memcpy(&Gs3_bits, &Gs3, sizeof(double));
        std::memcpy(&Gs1_golden_bits, &Gs1_golden, sizeof(double));
        std::memcpy(&Gs2_golden_bits, &Gs2_golden, sizeof(double));
        std::memcpy(&Gs3_golden_bits, &Gs3_golden, sizeof(double));
        // Compute bitwise difference and allow 1 LSB difference
        auto ulp_diff = [](uint64_t a, uint64_t b) -> uint64_t {
            // Handle sign bit: treat as unsigned distance in bit pattern
            return (a > b) ? (a - b) : (b - a);
        };
        bool fail_Gs1 = ulp_diff(Gs1_bits, Gs1_golden_bits) > MAX_GS13;
        bool fail_Gs2 = ulp_diff(Gs2_bits, Gs2_golden_bits) > MAX_GS13;
        bool fail_Gs3 = ulp_diff(Gs3_bits, Gs3_golden_bits) > MAX_GS13;
        max_diff = std::max(max_diff,
            std::max(ulp_diff(Gs1_bits, Gs1_golden_bits),
            std::max(ulp_diff(Gs2_bits, Gs2_golden_bits),
                     ulp_diff(Gs3_bits, Gs3_golden_bits))));
        if (fail_Gs1 || fail_Gs2 || fail_Gs3) {
            std::cout << std::setprecision(17);
            std::cout << "Mismatch at line " << (total+2) << ":\n";
            std::cout << "  beta=" << beta << ", X=" << X << std::endl;
            std::cout << "  Gs1: computed=" << Gs1 << " (0x" << std::hex << Gs1_bits << std::dec << "), golden=" << Gs1_golden << " (0x" << std::hex << Gs1_golden_bits << std::dec << ") | ulp_diff=" << ulp_diff(Gs1_bits, Gs1_golden_bits) << std::endl;
            std::cout << "  Gs2: computed=" << Gs2 << " (0x" << std::hex << Gs2_bits << std::dec << "), golden=" << Gs2_golden << " (0x" << std::hex << Gs2_golden_bits << std::dec << ") | ulp_diff=" << ulp_diff(Gs2_bits, Gs2_golden_bits) << std::endl;
            std::cout << "  Gs3: computed=" << Gs3 << " (0x" << std::hex << Gs3_bits << std::dec << "), golden=" << Gs3_golden << " (0x" << std::hex << Gs3_golden_bits << std::dec << ") | ulp_diff=" << ulp_diff(Gs3_bits, Gs3_golden_bits) << std::endl;
            fail++;
        }
        total++;
    }
    file.close();
    if (fail) {
        std::cout << "test_stiefel_Gs13 failed with " << fail << " mismatches out of " << total << "." << std::endl;
        std::cout << "Maximum ULP difference: " << max_diff << std::endl;
        std::exit(fail);
    } else {
        std::cout << "test_stiefel_Gs13 passed! (" << total << " cases, max_diff=" << max_diff << ")" << std::endl;
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

    test_stiefel_Gs13();

    // Now run kepler step
}