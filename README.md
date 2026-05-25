# OpenDSP
OpenDSP is an optimised for ARM toolbox library for audio digital signal processing. It features filtering and dynamic 
processing implementations.

# Features
Full floating-point operation, optimised with the NEON instruction set whenever possible.  
Every coefficients are computed dynamically, no LUTs are used.

### Filtering
IIR Filter with analog based filters approximated with the bilinear transform. All filters types (HPF, LPF, BPF) are
supported, 1st and 2nd order. There also is an implementation of an IIR Peak filter.

### Dynamics
Support for custom arbitrary transfer-function dynamics processing. Based on an ASR envelope, support down to the ms tuning
of attack, hold and release.
