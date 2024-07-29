# MarkLLPCode

For running sim, digi, reco clone MuC-Tutorial and a fork of mucoll-benchmarks (to incorporate LCRelations and Displaced Tracking options) into main directory:

```bash

git clone git@github.com:MuonColliderSoft/MuC-Tutorial.git
git clone git@github.com:mlarson02/mucoll-benchmarks-LLPs.git
mv mucoll-benchmarks-LLPs/ mucoll-benchmarks/
```

To incorporate HitSlimmer code from https://github.com/madbaron/MyBIBUtils/tree/master to use displaced tracking config:

```bash

git clone git@github.com:madbaron/MyBIBUtils.git
cd MyBIBUtils
mkdir build
cd build
cmake ../
make
export MARLIN_DLL=$(realpath lib/libMyBIBUtils.so):${MARLIN_DLL}
```
