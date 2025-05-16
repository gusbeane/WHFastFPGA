#include "kepler.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

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

int main()
{
    // First read in the initial conditions
    double x_vec[N_PLANETS], y_vec[N_PLANETS], z_vec[N_PLANETS];
    double vx_vec[N_PLANETS], vy_vec[N_PLANETS], vz_vec[N_PLANETS];
    double m_vec[N_PLANETS];
    double M0, dt;

    read_input(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, &M0, &dt);

    // printf("M0: %.20f, dt: %.20f\n", M0, dt);
    // printf("x: %.20f, %.20f, %.20f, %.20f, %.20f, %.20f, %.20f, %.20f\n", x_vec[0], x_vec[1], x_vec[2], x_vec[3], x_vec[4], x_vec[5], x_vec[6], x_vec[7]);

    kepler_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, M0, dt);

    verify_output(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec);

    // Now run kepler step
}