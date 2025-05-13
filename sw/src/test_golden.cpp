#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <chrono> // For timing
#include "whfastfpga.h"
#include "whfast.h"
#include "whfast_constants.h"
#include "whfast_kernel.h"
#include "util.h"

// Initial conditions copied from main.cpp
constexpr std::array<Body, N_BODIES> solarsystem_ics = {
    Body{{-0.008816286905115728, -0.0010954664916791675, 0.0002143249385447027},
         {0.00014315746073017681, -0.0004912441820893999, 8.127678560998346e-07},
         0.9999999999950272},
    Body{{-0.05942272929227954, -0.46308699693348293, -0.032897989948949075},
         {1.2978664284760637, -0.09524541469911743, -0.12677574364801253},
         1.6601208254808336e-07},
    Body{{-0.7276101005375593, 0.006575003332463933, 0.041795901908847084},
         {-0.019239782390457125, -1.1813975672919448, -0.01509392594251431},
         2.447838287784771e-06},
    Body{{-0.5789530452882667, -0.8361530119313055, 0.0002611520181901174},
         {0.8098712561282222, -0.5682496529341624, 2.6169897281383047e-05},
         3.0404326489511185e-06},
    Body{{-1.45202492400084, 0.827519404194876, 0.052981833432457694},
         {-0.37436417754222295, -0.6365841544564991, -0.004143932260467942},
         3.2271560828978514e-07},
    Body{{4.492983939852296, 2.0661626247490354, -0.10909246996001629},
         {-0.18818907783656452, 0.41919544951404614, 0.0024710497024424977},
         0.0009547919099366768},
    Body{{8.4974210980544, -4.8620394993693585, -0.2537835862373596},
         {0.14292308496870448, 0.2808676923735748, -0.010574288572728728},
         0.0002858856700231729},
    Body{{12.959111916929283, 14.760785302864473, -0.1130656917933948},
         {-0.1734971049470612, 0.14019515029516152, 0.0027683484887051457},
         4.366249613200406e-05},
    Body{{29.787987348666505, -2.51460654509393, -0.6347108842010732},
         {0.014142947617173336, 0.18292110872737416, -0.004092845767710294},
         5.151383772628957e-05}
};

bool compare_double(double a, double b, double epsilon = 1e-10) {
    return std::fabs(a - b) < epsilon;
}

void validate_csv(const std::string& filename, const std::array<Body, N_BODIES>& expected) {
    std::ifstream file("golden/" + filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    std::string line;
    size_t index = 0;
    while (std::getline(file, line)) {
        if (index >= expected.size()) {
            throw std::runtime_error("CSV file has more rows than expected.");
        }

        std::istringstream ss(line);
        std::string value;
        std::array<double, 7> values;
        size_t value_index = 0;

        while (std::getline(ss, value, ',')) {
            if (value_index >= values.size()) {
                throw std::runtime_error("CSV row has more columns than expected.");
            }
            values[value_index++] = std::stod(value);
        }

        if (value_index != values.size()) {
            throw std::runtime_error("CSV row has fewer columns than expected.");
        }

        const Body& expected_body = expected[index];
        if (!compare_double(values[0], expected_body.pos[0], 1.0E-7) ||
            !compare_double(values[1], expected_body.pos[1], 1.0E-7) ||
            !compare_double(values[2], expected_body.pos[2], 1.0E-7) ||
            !compare_double(values[3], expected_body.vel[0], 1.0E-7) ||
            !compare_double(values[4], expected_body.vel[1], 1.0E-7) ||
            !compare_double(values[5], expected_body.vel[2], 1.0E-7) ||
            !compare_double(values[6], expected_body.mass)) {
            std::cerr << "Mismatch found in file: " << filename << " at row " << index << ":\n"
                      << "Expected: " << expected_body.pos[0] << ", " << expected_body.pos[1] << ", " << expected_body.pos[2] << ", "
                      << expected_body.vel[0] << ", " << expected_body.vel[1] << ", " << expected_body.vel[2] << ", "
                      << expected_body.mass << "\n"
                      << "Found: " << values[0] << ", " << values[1] << ", " << values[2] << ", "
                      << values[3] << ", " << values[4] << ", " << values[5] << ", " << values[6] << std::endl;
            throw std::runtime_error("Mismatch found in file: " + filename);
        }

        ++index;
    }

    if (index != expected.size()) {
        throw std::runtime_error("CSV file has fewer rows than expected.");
    }

    file.close();
}

void test_golden(double tmax_inyr, double dt, const std::string& filename) {
    auto start = std::chrono::high_resolution_clock::now(); // Start timing

    std::cout << "Testing: " << filename << "... ";

    // Generate integration
    std::array<Body, N_BODIES> solarsystem = solarsystem_ics;
    move_to_center_of_mass(solarsystem);

    Body com;
    double tmax = 2.0 * M_PI * tmax_inyr;
    long Nint = static_cast<long>(tmax / dt);
    whfast_integrate(solarsystem, &com, dt, Nint);

    // Validate the CSV file against the computed values
    validate_csv(filename, solarsystem);

    auto end = std::chrono::high_resolution_clock::now(); // End timing
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Done. Time taken: " << elapsed.count() << " seconds." << std::endl;
}

// Single drift step
void test_golden_drift(double dt, const std::string &filename)
{
     auto start = std::chrono::high_resolution_clock::now(); // Start timing

     std::cout << "Generating:" << filename << "... ";
     // Generate 100 yr integration
     std::array<Body, N_BODIES> solarsystem = solarsystem_ics;
     move_to_center_of_mass(solarsystem);

     Body com;
     // Prepare AVX-512 vectors for bodies 1-8 (ignore solarsystem[0])
     __m512d x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec;

     inertial_to_democraticheliocentric_posvel(solarsystem, &com, &x_vec, &y_vec, &z_vec,
                                               &vx_vec, &vy_vec, &vz_vec, &m_vec);

     // Calculate necessary constants
     // Constructs a struct called kConsts with the constants
     initialize_constants(solarsystem[0].mass, m_vec);

     // Do drift step
     whfast_drift_step(&x_vec, &y_vec, &z_vec, &vx_vec, &vy_vec, &vz_vec, m_vec, &com, dt);

     // Convert back to inertial coordinates and store into solarsystem
     democraticheliocentric_to_inertial_posvel(solarsystem, &com, x_vec, y_vec, z_vec,
                                               vx_vec, vy_vec, vz_vec, m_vec);

     // Output mass, position, and velocities of the particles
     validate_csv(filename, solarsystem);

     auto end = std::chrono::high_resolution_clock::now(); // End timing
     std::chrono::duration<double> elapsed = end - start;
     std::cout << "Done. Time taken: " << elapsed.count() << " seconds." << std::endl;
}

int main() {
    try {
        double dt = 5.0 / 365.25 * 2 * M_PI; // 5 days
        test_golden(1e2, dt, "solarsystem_100yr.csv");
        test_golden(1e3, dt, "solarsystem_1kyr.csv");
        test_golden(1e4, dt, "solarsystem_10kyr.csv");
        test_golden_drift(dt, "solarsystem_1drift.csv");
        std::cout << "All golden vectors validated successfully." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Validation failed: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
