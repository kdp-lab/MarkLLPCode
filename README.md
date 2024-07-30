# MarkLLPCode

For running sim, digi, reco clone MuC-Tutorial and a fork of mucoll-benchmarks (to incorporate LCRelations and Displaced Tracking options) into main directory:

```bash

git clone git@github.com:MuonColliderSoft/MuC-Tutorial.git
git clone git@github.com:mlarson02/mucoll-benchmarks-LLPs.git
mv mucoll-benchmarks-LLPs/ mucoll-benchmarks/
```

To incorporate HitSlimmer code from https://github.com/madbaron/MyBIBUtils/tree/master to use displaced tracking config
Note: need to be sourced to cvmfs to compile, and need to export MARLIN_DLL each time before running reco (will add to multi_digi_reco.py script):

```bash

git clone git@github.com:madbaron/MyBIBUtils.git
cd MyBIBUtils
mkdir build
cd build
cmake ../
make
export MARLIN_DLL=$(realpath lib/libMyBIBUtils.so):${MARLIN_DLL}
```
To incorporate TrackSelector code from https://github.com/madbaron/MarlinTrkProcessors/tree/selector_updates to use displaced tracking config (note: need to be sourced to cvmfs to compile): 

```bash

git clone git@github.com:madbaron/MarlinTrkProcessors.git
git checkout selector_updates
cd MarlinTrkProcessors
mkdir build
cd build
cmake ../
make
export MARLIN_DLL=$(realpath lib/libMarlinTrkProcessors.so):${MARLIN_DLL}
```
