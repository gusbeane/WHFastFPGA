#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <array>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <iomanip>
#include <chrono> // For timing
#include <sstream>
#include "whfast.h"
#include "whfast_kernel.h"
#include "whfast_constants.h"
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

double compute_energy(const std::array<Body, N_BODIES>& solarsystem, const Body& com) {
    double total_energy = 0.0;
    
    double x, y, z;
    double vx, vy, vz;

    // compute kinetic energy
    for (int i = 0; i < N_BODIES; ++i) {
        x = solarsystem[i].pos[0] - com.pos[0];
        y = solarsystem[i].pos[1] - com.pos[1];
        z = solarsystem[i].pos[2] - com.pos[2];
        vx = solarsystem[i].vel[0] - com.vel[0];
        vy = solarsystem[i].vel[1] - com.vel[1];
        vz = solarsystem[i].vel[2] - com.vel[2];

        double kinetic_energy = 0.5 * solarsystem[i].mass * (vx * vx + vy * vy + vz * vz);
        total_energy += kinetic_energy;
    }

    // compute potential energy
    for (int i = 0; i < N_BODIES; ++i) {
        for (int j = i + 1; j < N_BODIES; ++j) {
            double dx = solarsystem[i].pos[0] - solarsystem[j].pos[0];
            double dy = solarsystem[i].pos[1] - solarsystem[j].pos[1];
            double dz = solarsystem[i].pos[2] - solarsystem[j].pos[2];
            double r = sqrt(dx * dx + dy * dy + dz * dz);
            total_energy -= (solarsystem[i].mass * solarsystem[j].mass) / r;
        }
    }
    // std::cout << "Total energy: " << total_energy << " (" << double_to_hex(total_energy) << ")" << std::endl;
    return total_energy;
}

double gen_golden(double tmax_inyr, double dt, const std::string& filename)
{
    auto start = std::chrono::high_resolution_clock::now(); // Start timing

    std::cout << "Generating:" << filename << "... ";
    // Generate 100 yr integration
    std::array<Body, N_BODIES> solarsystem = solarsystem_ics;
    move_to_center_of_mass(solarsystem);
    
    Body com;
    double tmax = 2.0 * M_PI * tmax_inyr; // 100 yr
    long Nint = static_cast<long>(tmax / dt);
    std::cout << "Nint=" << Nint << "  ";
    if(Nint > 0)
        whfast_integrate(solarsystem, &com, dt, Nint);

    // Output mass, position, and velocities of the particles
    output_to_csv(solarsystem, filename);

    auto end = std::chrono::high_resolution_clock::now(); // End timing
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Done. Time taken: " << elapsed.count() << " seconds." << std::endl;

    return compute_energy(solarsystem, com);
}

void gen_integrate_step_csv(const std::string &filename)
{
    std::ofstream file("golden/" + filename);
    std::cout << "Generating:" << filename << "... ";
    // Generate 100 yr integration
    std::array<Body, N_BODIES> solarsystem = solarsystem_ics;
    move_to_center_of_mass(solarsystem);

    double x_vec[N_PLANETS], y_vec[N_PLANETS], z_vec[N_PLANETS];
    double vx_vec[N_PLANETS], vy_vec[N_PLANETS], vz_vec[N_PLANETS];
    double m_vec[N_PLANETS];
    Body com;
    double dt = 5.0 / 365.25 * 2 * M_PI; // 5 days in radians
    double dt_half = dt / 2.0;

    inertial_to_democraticheliocentric_posvel(solarsystem, &com, x_vec, y_vec, z_vec,
                                              vx_vec, vy_vec, vz_vec, m_vec);

    initialize_constants(solarsystem[0].mass, m_vec);

    // Print initial conditions to file
    file << "x,y,z,vx,vy,vz,m\n";
    file << "dt=" << double_to_hex(dt) << "\n";
    file << "dt_half=" << double_to_hex(dt_half) << "\n";
    file << "M0=" << double_to_hex(kConsts->M0) << "\n";
    for (int i = 0; i < N_PLANETS; ++i)
    {
        file << double_to_hex(x_vec[i]) << ","
             << double_to_hex(y_vec[i]) << ","
             << double_to_hex(z_vec[i]) << ","
             << double_to_hex(vx_vec[i]) << ","
             << double_to_hex(vy_vec[i]) << ","
             << double_to_hex(vz_vec[i]) << ","
             << double_to_hex(m_vec[i]) << "\n";
    }

    whfast_kepler_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, dt_half);

    // Print results to file_kepler
    for (int i = 0; i < N_PLANETS; ++i)
    {
        file << double_to_hex(x_vec[i]) << ","
             << double_to_hex(y_vec[i]) << ","
             << double_to_hex(z_vec[i]) << ","
             << double_to_hex(vx_vec[i]) << ","
             << double_to_hex(vy_vec[i]) << ","
             << double_to_hex(vz_vec[i]) << ","
             << double_to_hex(m_vec[i]) << "\n";
    }

    // Do jump step
    whfast_jump_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt_half);

    // Print results to file_jump
    for (int i = 0; i < N_PLANETS; ++i)
    {
        file << double_to_hex(x_vec[i]) << ","
             << double_to_hex(y_vec[i]) << ","
             << double_to_hex(z_vec[i]) << ","
             << double_to_hex(vx_vec[i]) << ","
             << double_to_hex(vy_vec[i]) << ","
             << double_to_hex(vz_vec[i]) << ","
             << double_to_hex(m_vec[i]) << "\n";
    }

    // Do interaction step
    whfast_interaction_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt);

    // Print results to file
    for (int i = 0; i < N_PLANETS; ++i)
    {
        file << double_to_hex(x_vec[i]) << ","
             << double_to_hex(y_vec[i]) << ","
             << double_to_hex(z_vec[i]) << ","
             << double_to_hex(vx_vec[i]) << ","
             << double_to_hex(vy_vec[i]) << ","
             << double_to_hex(vz_vec[i]) << ","
             << double_to_hex(m_vec[i]) << "\n";
    }

    file.close();

    std::cout << "Done." << std::endl;
}

// Remove extern "C" and use regular C++ function declarations
void stiefel_Gs03(double* Gs0, double* Gs1, double* Gs2, double* Gs3, double beta, double X);
void stiefel_Gs13(double* Gs1, double* Gs2, double* Gs3, double beta, double X);

#define BETA_MIN 0.01
#define BETA_MAX 3.0
#define X_MIN 0.001
#define X_MAX 0.3
#define N_BETA 100
#define N_X 100

void gen_stiefel_Gs03_csv(const std::string& filename) {
    std::ofstream file("golden/" + filename);
    file << "beta,X,Gs0,Gs1,Gs2,Gs3\n";
    for (int i = 0; i < N_BETA; ++i) {
        double beta = BETA_MIN + (BETA_MAX - BETA_MIN) * i / (N_BETA - 1);
        for (int j = 0; j < N_X; ++j) {
            double X = X_MIN + (X_MAX - X_MIN) * j / (N_X - 1);
            double Gs0, Gs1, Gs2, Gs3;
            stiefel_Gs03(&Gs0, &Gs1, &Gs2, &Gs3, beta, X);
            file << double_to_hex(beta) << "," << double_to_hex(X) << ","
                 << double_to_hex(Gs0) << "," << double_to_hex(Gs1) << ","
                 << double_to_hex(Gs2) << "," << double_to_hex(Gs3) << "\n";
        }
    }
    file.close();
}

void gen_stiefel_Gs13_csv(const std::string& filename) {
    std::ofstream file("golden/" + filename);
    file << "beta,X,Gs1,Gs2,Gs3\n";
    for (int i = 0; i < N_BETA; ++i) {
        double beta = BETA_MIN + (BETA_MAX - BETA_MIN) * i / (N_BETA - 1);
        for (int j = 0; j < N_X; ++j) {
            double X = X_MIN + (X_MAX - X_MIN) * j / (N_X - 1);
            double Gs1, Gs2, Gs3;
            stiefel_Gs13(&Gs1, &Gs2, &Gs3, beta, X);
            file << double_to_hex(beta) << "," << double_to_hex(X) << ","
                 << double_to_hex(Gs1) << "," << double_to_hex(Gs2) << ","
                 << double_to_hex(Gs3) << "\n";
        }
    }
    file.close();
}

#define N_NEWTON_HALLEY 8
#define R0_MIN 0.25
#define R0_MAX 50.0
#define ETA0_MIN -0.25
#define ETA0_MAX 0.25
#define ZETA0_MIN -0.25
#define ZETA0_MAX 0.25

// Generates a CSV of input vectors and outputs for newton_step and halley_step
void gen_newton_halley_csv(const std::string &filename_newton, const std::string &filename_halley)
{
    std::ofstream file_newton("golden/" + filename_newton);
    std::ofstream file_halley("golden/" + filename_halley);
    
    file_newton << "X,beta,r0,eta0,zeta0,Xout\n";
    file_halley << "X,beta,r0,eta0,zeta0,Xout\n";

    const double dt = 5.0 / 365.25 * 2 * M_PI;
    file_newton << "dt=" << double_to_hex(dt) << "\n";
    file_halley << "dt=" << double_to_hex(dt) << "\n";

    for (int i = 0; i < N_NEWTON_HALLEY; ++i) {
        double r0 = R0_MIN + (R0_MAX - R0_MIN) * i / (N_NEWTON_HALLEY - 1);
        for (int j = 0; j < N_NEWTON_HALLEY; ++j) {
            double eta0 = ETA0_MIN + (ETA0_MAX - ETA0_MIN) * j / (N_NEWTON_HALLEY - 1);
            for (int k = 0; k < N_NEWTON_HALLEY; ++k) {
                double zeta0 = ZETA0_MIN + (ZETA0_MAX - ZETA0_MIN) * k / (N_NEWTON_HALLEY - 1);
                for (int l = 0; l < N_NEWTON_HALLEY; ++l) {
                    double beta = BETA_MIN + (BETA_MAX - BETA_MIN) * l / (N_NEWTON_HALLEY - 1);
                    for (int m = 0; m < N_NEWTON_HALLEY; ++m) {
                        double X = X_MIN + (X_MAX - X_MIN) * m / (N_NEWTON_HALLEY - 1);
                        double X_newton = X;
                        double X_halley = X;
                        // Call the step functions
                        extern void newton_step(double*, double, double, double, double, double);
                        extern void halley_step(double*, double, double, double, double, double);
                        
                        newton_step(&X_newton, beta, r0, eta0, zeta0, dt);
                        file_newton << double_to_hex(X) << ","
                             << double_to_hex(beta) << ","
                             << double_to_hex(r0) << ","
                             << double_to_hex(eta0) << ","
                             << double_to_hex(zeta0) << ","
                             << double_to_hex(X_newton) << "\n";
                        
                        halley_step(&X_halley, beta, r0, eta0, zeta0, dt);
                        file_halley << double_to_hex(X) << ","
                             << double_to_hex(beta) << ","
                             << double_to_hex(r0) << ","
                             << double_to_hex(eta0) << ","
                             << double_to_hex(zeta0) << ","
                             << double_to_hex(X_halley) << "\n";
                    }
                }
            }
        }
    }
    file_newton.close();
    file_halley.close();
}

void gen_kernel_step_csv(const std::string &filename)
{
    std::ofstream file("golden/" + filename);
    std::cout << "Generating:" << filename << "... ";
    // Generate 100 yr integration
    std::array<Body, N_BODIES> solarsystem = solarsystem_ics;
    move_to_center_of_mass(solarsystem);

    double x_vec[N_PLANETS], y_vec[N_PLANETS], z_vec[N_PLANETS];
    double vx_vec[N_PLANETS], vy_vec[N_PLANETS], vz_vec[N_PLANETS];
    double m_vec[N_PLANETS];
    Body com;
    double dt = 5.0 / 365.25 * 2 * M_PI; // 5 days in radians
    double dt_half = dt / 2.0;
    long Nint = 100;

    inertial_to_democraticheliocentric_posvel(solarsystem, &com, x_vec, y_vec, z_vec,
                                              vx_vec, vy_vec, vz_vec, m_vec);

    initialize_constants(solarsystem[0].mass, m_vec);

    // First do half-drift step
    whfast_drift_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, &com, dt_half);

    // Print initial conditions to file
    file << "x,y,z,vx,vy,vz,m\n";
    file << "dt=" << double_to_hex(dt) << "\n";
    file << "M0=" << double_to_hex(kConsts->M0) << "\n";
    file << "Nint=" << long_to_hex(Nint) << "\n";
    for (int i = 0; i < N_PLANETS; ++i)
    {
        file << double_to_hex(x_vec[i]) << ","
             << double_to_hex(y_vec[i]) << ","
             << double_to_hex(z_vec[i]) << ","
             << double_to_hex(vx_vec[i]) << ","
             << double_to_hex(vy_vec[i]) << ","
             << double_to_hex(vz_vec[i]) << ","
             << double_to_hex(m_vec[i]) << "\n";
    }

    // Now call the kernel Nint times
    // In whfast_integrate, the kernel is called Nint-1 times
    whfast_kernel(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, &com, dt, Nint);

    // Print results
    for (int i = 0; i < N_PLANETS; ++i)
    {
        file << double_to_hex(x_vec[i]) << ","
             << double_to_hex(y_vec[i]) << ","
             << double_to_hex(z_vec[i]) << ","
             << double_to_hex(vx_vec[i]) << ","
             << double_to_hex(vy_vec[i]) << ","
             << double_to_hex(vz_vec[i]) << ","
             << double_to_hex(m_vec[i]) << "\n";
    }

    file.close();

    std::cout << "Done." << std::endl;
}

int main() {
    double dt = 5.0 / 365.25 * 2 * M_PI; // 5 days
    double energy_ics  = gen_golden(0, dt, "solarsystem_ic.csv");

    double energy_1e2 = gen_golden(1e2, dt, "solarsystem_100yr.csv");

    double energy_1e3 = gen_golden(1e3, dt, "solarsystem_1kyr.csv");

    double energy_1e4 = gen_golden(1e4, dt, "solarsystem_10kyr.csv");

    // double energy_2e4 = gen_golden(2e4, dt, "solarsystem_20kyr.csv");
    
    // double energy_4e4 = gen_golden(4e4, dt, "solarsystem_40kyr.csv");
    
    // double energy_1e5 = gen_golden(1e5, dt, "solarsystem_100kyr.csv");

    std::cout << "Energy difference (100 yr - ics): " << (energy_1e2 - energy_ics)/energy_ics << std::endl;
    std::cout << "Energy difference (  1 kyr - ics): " << (energy_1e3 - energy_ics)/energy_ics << std::endl;
    std::cout << "Energy difference ( 10 kyr - ics): " << (energy_1e4 - energy_ics)/energy_ics << std::endl;
    // std::cout << "Energy difference ( 20 kyr - ics): " << (energy_2e4 - energy_ics)/energy_ics << std::endl;
    // std::cout << "Energy difference ( 40 kyr - ics): " << (energy_4e4 - energy_ics)/energy_ics << std::endl;
    // std::cout << "Energy difference (100 kyr - ics): " << (energy_1e5 - energy_ics)/energy_ics << std::endl;

    // Generate Stiefel function golden CSVs
    gen_stiefel_Gs03_csv("golden_stiefel_Gs03.csv");
    gen_stiefel_Gs13_csv("golden_stiefel_Gs13.csv");

    // Generate Newton-Halley input CSV
    gen_newton_halley_csv("golden_newton_step.csv", "golden_halley_step.csv");

    // Generate integration step CSV
    gen_integrate_step_csv("golden_integrate_step.csv");

    // Generate kernel step CSV
    gen_kernel_step_csv("golden_kernel_step.csv");

    std::cout << "All golden vectors made." << std::endl;

    std::cout << "All Stiefel Gs03 and Gs13 samples made." << std::endl;
    return 0;
}
