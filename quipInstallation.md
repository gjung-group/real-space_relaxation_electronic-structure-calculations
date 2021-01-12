To install the quip machine learning potential with Lammps, a couple of notes:
- Install quip first, do not forget to activate the GAP library while doing so. Details here. https://libatoms.github.io/GAP/installation.html. Important, use gfortran compiler. Check the lammps specific section to turn quip into a library.
- Install LAMMPS from source as usual, but add the path to the quip library you have just created. 
```
cmake -D PKG_MANYBODY=yes -D PKG_USER-MISC=yes -D PKG_MOLECULE=yes -D PKG_USER-QUIP=yes -D QUIP_LIBRARY=/home/nleconte/github/QUIP/build/linux_x86_64_gfortran/libquip.a  ../cmake
```
Do not mix yes, on, etc within cmake. gfortran should be your standard compiler chosen by cmake, do not change it, otherwise you’ll have trouble interfacing with quip. In case you really need it, you need to activate it with the following 
```
-D -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpiicpc -D CMAKE_Fortran_COMPILER=mpiifort
```
but you might have to recompile quip with a different compiler. Unfortunately, several attemps at quip compilation using intel compiler failed for me (both icc and gcc).
- To run the calculation, e.g. using drip for interlayer interactions, and tersoff for B N interactions:
```
 pair_style hybrid/overlay drip quip extep
 pair_coeff * * drip  C.drip   B N C
 pair_coeff * * quip  Carbon_GAP_20_potential/Carbon_GAP_20.xml "" NULL NULL 6
 pair_coeff  * * extep BN.extep  B  N NULL # chemical
 ```
- The Carbon_GAP_20_potential folder, I have a local copy. The usual open database is not available these days, but you can download it from https://doi.org/10.17863/CAM.54529.
