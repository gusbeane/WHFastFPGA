# WHFastFPGA

> Sandbox project for learning how to write FPGA kernels using a basic symplectic solar system integrator.

## Overview

This project implements the WHFast integrator as described in [Rein and Tamayo, 2015](https://arxiv.org/abs/1506.01084), starting from the WHFast512 implementation that uses AVX512 instructions ([Javaheri et al., 2023](https://arxiv.org/abs/2307.05683)). I began by forking REBOUND ([hannorein/rebound](https://github.com/hannorein/rebound)) and extracting the WHFast512 integrator, removing the AVX512-specific code to create a cleaner baseline for experimentation.

My initial goal is to learn how to write an FPGA kernel using Vitis HLS and the corresponding host-side program. Once that foundation is working, I plan to explore several performance-related improvements, including:

- Custom bit-widths to efficiently integrate systems that require more than double precision
- Pipelining to support concurrent integration of large numbers of planetary systems
- Algorithm rearchitecting to minimize latency

## Repository Structure

```
/
├── sw/         # C++ version of WHFast (WHFast512 with AVX512 instructions removed)
├── hw/         # FPGA implementation
├── sw512/      # C++ version of WHFast512 (archival; not actively developed)
└── README.md   # Project overview
```

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).