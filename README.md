# BdNMC

This code requires the BdNMC_LANL submodule. All submodules can be grabbed with "git submodule update --init --recursive".

Full documentation can be found in the appendix of http://arxiv.org/abs/1609.01770.

To run the code, execute ./bin/BDNMC \<filepath\> in the base BdNMC directory, where the filepath points to a valid parameter file. This will compile and run the code with the supplied parameter file. If no \<filepath\> is specified, it will use parameter.dat found in the base BdNMC directory.

A sample parameter file is included: paramater.dat.
