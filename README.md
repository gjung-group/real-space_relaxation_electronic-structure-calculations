# Content

Also available in lammps/potential folder:
- CBN_RPA.drip: DRIP parameters based on EXX-RPA DFT data for use with LAMMPS.  
- CBN_LDA.drip: DRIP parameters based on LDA DFT data for use with LAMMPS. 
- CH.rebo-LB: Reparametrization of the CH.rebo file, for use with LAMMPS.

Input file:
- lammps.in.gap2020: input file for LAMMPS that leverages the gap2020 potential for intralayer interactions.

Output scripts:
- `getOutput.py` and `getOutputOnlyG.py`: extract data from LAMMPS `dump.minimization` file
- getInterlayerDistace.py: extract interlayer distances from files in previous step
- getDisplacements.py: extract displacements from files in first step

# Future content
Code:
- Gendata.in: input file to perform bandstructure calculation
- grabnes: binary to perform the bandstructure calculation
 
