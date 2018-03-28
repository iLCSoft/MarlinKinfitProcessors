# v00-04

* 2017-12-15 Graham Wilson ([PR#4](https://github.com/iLCSoft/MarlinKinfitProcessors/pull/4))
  - Add steering files for J/psi -> mu+ mu-, Higgs -> mu+ mu-, Higgs -> mu+mu-mu+mu-, K+->pi+pi-pi+
  - eta-> pi+ pi- gamma, eta-> pi+pi-pi0
  - Fix most - but not yet all - of the compiler warnings in the MassConstraintFitter related processors.

* 2017-11-29 Graham Wilson ([PR#3](https://github.com/iLCSoft/MarlinKinfitProcessors/pull/3))
  - First commit of code of MassConstraintFitter from Justin Anguiano for 
    generalized mass-constrained fitting -- Graham. 
    Example .xml steering file is for simple use case of pi0-> gamma gamma.

* 2017-11-23 Graham Wilson ([PR#2](https://github.com/iLCSoft/MarlinKinfitProcessors/pull/2))
  - add author

* 2017-12-20 Graham Wilson ([PR#5](https://github.com/iLCSoft/MarlinKinfitProcessors/pull/5))
  - All compiler warnings now fixed.
  - Pi0Fitter.xml now configured and tested for v01-19-05-p01 test production.
  - a) angular smearing now turned off (was work-around for photon position issues in v01-19-04)
  - b) energy scale set to 1.0
  - c) additional processor parameters in MCParticleFilter to deal with generator status differences

# v00-03

# v00-02

J. List
   - added ZH5CFit processor

F. Gaede
   - made compatible with c++11
   - removed -ansi -pedantic -Wno-long-long
  
# v00-01

 moved example processors to this new package in order to
 reduce cross-dependencies with MarlinReco (eg via GammaGammaFinder).
