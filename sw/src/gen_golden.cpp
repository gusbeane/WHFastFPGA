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
#include "whfastfpga.h"
#include "whfast.h"
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
    whfast_integrate(solarsystem, &com, dt, Nint);

    // Output mass, position, and velocities of the particles
    output_to_csv(solarsystem, filename);

    auto end = std::chrono::high_resolution_clock::now(); // End timing
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Done. Time taken: " << elapsed.count() << " seconds." << std::endl;
}

// Remove extern "C" and use regular C++ function declarations
void stiefel_Gs03(double* Gs0, double* Gs1, double* Gs2, double* Gs3, double beta, double X);
void stiefel_Gs13(double* Gs1, double* Gs2, double* Gs3, double beta, double X);

// Helper to get hex string of a double
std::string double_to_hex(double d) {
    union { double d; uint64_t u; } u;
    u.d = d;
    std::ostringstream oss;
    oss << std::hex << std::setw(16) << std::setfill('0') << u.u;
    return oss.str();
}

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


int main() {
    double dt = 5.0 / 365.25 * 2 * M_PI; // 5 days
    gen_golden(1e2, dt, "solarsystem_100yr.csv");

    gen_golden(1e3, dt, "solarsystem_1kyr.csv");

    gen_golden(1e4, dt, "solarsystem_10kyr.csv");

    // Generate Stiefel function golden CSVs
    gen_stiefel_Gs03_csv("golden_stiefel_Gs03.csv");
    gen_stiefel_Gs13_csv("golden_stiefel_Gs13.csv");

    std::cout << "All golden vectors made." << std::endl;

    std::cout << "All Stiefel Gs03 and Gs13 samples made." << std::endl;
    return 0;
}
