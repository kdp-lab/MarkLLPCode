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

```

After cloning, to fix a small issue with MyBibUtils/src/HitSlimmer.cc, change 'itHit' to 'jitHit' on line 169 

```bash
cd MyBIBUtils
mkdir build
cd build
cmake ../
make
export MARLIN_DLL=$(realpath lib/libMyBIBUtils.so):${MARLIN_DLL}
```

For running on osg cluster, to setup environment for either sim, digi, or reco use: 
```bash
apptainer run \
    -B /scratch/${USER} \
    -B /ospool/uc-shared/project/muoncollider/ \
    /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/muon-collider/mucoll-deploy/mucoll:2.9-alma9

source /opt/setup_mucoll.sh
```
