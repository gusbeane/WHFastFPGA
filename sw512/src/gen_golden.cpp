#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <array>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <iomanip>
#include <chrono> // For timing
#include "whfastfpga.h"
#include "whfast512.h"
#include "whfast512_constants.h"
#include "whfast512_kernel.h" // Ensure this header contains the declaration of whfast512_drift_step
#include "util.h"
#include <immintrin.h> // For AVX-512

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

void output_to_csv(const std::array<Body, N_BODIES>& solarsystem, const std::string& filename) {
    // Create directory "golden/" if it does not exist
    struct stat info;
    if (stat("golden", &info) != 0) {
        if (mkdir("golden", 0777) != 0) {
            throw std::runtime_error("Failed to create directory 'golden'");
        }
    }

    // Open the file for writing
    std::ofstream file("golden/" + filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    // Set precision for floating-point numbers
    file << std::fixed << std::setprecision(12);

    // Write data to the CSV file
    for (const auto& body : solarsystem) {
        file << body.pos[0] << "," << body.pos[1] << "," << body.pos[2] << ","
             << body.vel[0] << "," << body.vel[1] << "," << body.vel[2] << ","
             << body.mass << "\n";
    }

    file.close();
}

// Rewrite to accept AVX-512 vectors and extract elements for CSV output
void output_vec_to_csv(double M0, double dt, __m512d x_vec, __m512d y_vec, __m512d z_vec,
                       __m512d vx_vec, __m512d vy_vec, __m512d vz_vec,
                       __m512d m_vec, const std::string &filename)
{
     // Create directory "golden/" if it does not exist
     struct stat info;
     if (stat("golden", &info) != 0)
     {
          if (mkdir("golden", 0777) != 0)
          {
               throw std::runtime_error("Failed to create directory 'golden'");
          }
     }

     // Open the file for writing
     std::ofstream file("golden/" + filename);
     if (!file.is_open())
     {
          throw std::runtime_error("Failed to open file: " + filename);
     }

     // Set precision for floating-point numbers
     file << std::fixed << std::setprecision(20);

     // Extract vector elements into arrays
     alignas(64) double x[N_BODIES-1], y[N_BODIES-1], z[N_BODIES-1];
     alignas(64) double vx[N_BODIES-1], vy[N_BODIES-1], vz[N_BODIES-1], m[N_BODIES-1];
     _mm512_store_pd(x, x_vec);
     _mm512_store_pd(y, y_vec);
     _mm512_store_pd(z, z_vec);
     _mm512_store_pd(vx, vx_vec);
     _mm512_store_pd(vy, vy_vec);
     _mm512_store_pd(vz, vz_vec);
     _mm512_store_pd(m, m_vec);

     // Write data to the CSV file
     file << M0 << "," << dt << "\n"; // Write M0 as the first line
     for (int i = 0; i < N_BODIES-1; ++i) {
         file << x[i] << "," << y[i] << "," << z[i] << ","
              << vx[i] << "," << vy[i] << "," << vz[i] << "," << m[i] << "\n";
     }

     file.close();
}

void gen_golden(double tmax_inyr, double dt, const std::string& filename)
{
    auto start = std::chrono::high_resolution_clock::now(); // Start timing

    std::cout << "Generating:" << filename << "... ";
    // Generate 100 yr integration
    std::array<Body, N_BODIES> solarsystem = solarsystem_ics;
    move_to_center_of_mass(solarsystem);
    
    Body com;
    double tmax = 2.0 * M_PI * tmax_inyr; // 100 yr
    long Nint = static_cast<long>(tmax / dt);
    whfast512_integrate(solarsystem, &com, dt, Nint);

    // Output mass, position, and velocities of the particles
    output_to_csv(solarsystem, filename);

    auto end = std::chrono::high_resolution_clock::now(); // End timing
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Done. Time taken: " << elapsed.count() << " seconds." << std::endl;
}

// Single drift step
void gen_golden_drift(double dt, const std::string &filename)
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

     // Do drift step, first outputting the input vectors
     output_vec_to_csv(kConsts->M0, dt, x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, "golden_drift_inputvectors.csv");

     whfast512_drift_step(&x_vec, &y_vec, &z_vec, &vx_vec, &vy_vec, &vz_vec, m_vec, &com, dt);
     
     output_vec_to_csv(kConsts->M0, dt, x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, "golden_drift_outputvectors.csv");

     // Convert back to inertial coordinates and store into solarsystem
     democraticheliocentric_to_inertial_posvel(solarsystem, &com, x_vec, y_vec, z_vec,
                                               vx_vec, vy_vec, vz_vec, m_vec);

     // Output mass, position, and velocities of the particles
     output_to_csv(solarsystem, filename);

     auto end = std::chrono::high_resolution_clock::now(); // End timing
     std::chrono::duration<double> elapsed = end - start;
     std::cout << "Done. Time taken: " << elapsed.count() << " seconds." << std::endl;
}

// Single drift step
void print_kepler_timing(double dt)
{

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

     auto start = std::chrono::high_resolution_clock::now(); // Start timing
     for(int i=0; i<10000000; i++)
          whfast512_kepler_step(&x_vec, &y_vec, &z_vec, &vx_vec, &vy_vec, &vz_vec, m_vec, dt);
     auto end = std::chrono::high_resolution_clock::now(); // End timing
     std::chrono::duration<double> elapsed = end - start;

     std::cout << "Average time per kepler step taken: " << 1e9*elapsed.count()/10000000.0 << " nanoseconds." << std::endl;

     // Now get timing for Stiefel function
     double X = 0.0925188;
     double beta = 2.58286;

     __m512d X_ = _mm512_set1_pd(X);
     __m512d beta_ = _mm512_set1_pd(beta);

     __m512d Gs1;
     __m512d Gs2;
     __m512d Gs3;

     start = std::chrono::high_resolution_clock::now(); // Start timing
     for(int i=0; i<10000000; i++)
          mm_stiefel_Gs13_avx512(&Gs1, &Gs2, &Gs3, beta_, X_);     
     end = std::chrono::high_resolution_clock::now(); // End timing
     elapsed = end - start;

     std::cout << "Average time per stiefel Gs13 call taken: " << 1e9*elapsed.count()/10000000.0 << " nanoseconds." << std::endl;
}

int main() {
    double dt = 5.0 / 365.25 * 2 * M_PI; // 5 days
    gen_golden(1e2, dt, "solarsystem_100yr.csv");

    gen_golden(1e3, dt, "solarsystem_1kyr.csv");

    gen_golden(1e4, dt, "solarsystem_10kyr.csv");
    
    gen_golden_drift(dt, "solarsystem_1drift.csv");

    print_kepler_timing(dt);

    std::cout << "All golden vectors made." << std::endl;
    return 0;
}
