#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <array>
#include <iostream>
#include "whfastfpga.h"
#include "whfast512.h"
#include "whfast512_constants.h"
#include "util.h"
#include <chrono> // For timing

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
         5.151383772628957e-05}};

int main(int argc, char **argv)
{
      std::array<Body, N_BODIES> solarsystem = solarsystem_ics;
      move_to_center_of_mass(solarsystem);

      Body com;

      double tmax = 2.0 * M_PI * 1e5;      // 1 kyr
      double dt = 5.0 / 365.25 * 2 * M_PI; // 5 days
      long Nint = static_cast<long>(tmax / dt);

      auto start = std::chrono::high_resolution_clock::now(); // Start timing
      whfast512_integrate(solarsystem, &com, dt, Nint);
      auto end = std::chrono::high_resolution_clock::now(); // End timing
      std::chrono::duration<double> elapsed = end - start;
      
      printf("Time taken: %f seconds.\n", elapsed.count());
      printf("Average time per step taken: %f nanoseconds.\n", 1e9 * elapsed.count() / Nint);
      printf("\n");
}
