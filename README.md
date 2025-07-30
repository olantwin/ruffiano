# ruffiano

## Introduction

Simulation framework for the SHiP charm cross-section measurement. This package is pure python, but depends on the `ruffiano` branch of [`FairShip`](https://github.com/ShipSoft/FairShip).

For now based on the the upstream `run_simScript.py` with irrelevant parts
removed and the legacy SHiP-Charm geometry (the muon flux geometry actually). 

The way the geometry is specified will undergo further streamlining.

## Getting started

To set up the environment, source the SHiP development software stack:

```
source /cvmfs/ship.cern.ch/dev-root-master/setUp.sh
```

Clone `FairShip`:

```
git clone https://github.com/ShipSoft/FairShip
cd FairShip; git checkout ruffiano; cd ..
```

Build a development version of `FairShip`:

```
bits build FairShip --default release --always-prefer-system --config-dir $SHIPDIST
```

Enter the environment:

```
bits enter FairShip/latest-ruffiano-release
```

Finally, run the simulation from this repo:

```
python run_simScript.py -h
usage: run_simScript.py [-h] [--Pythia8] [--PG] [--pID PID] [--Estart ESTART] [--Eend EEND] [--FixedTarget] [--Ntuple] [--MuonBack] [--FollowMuon] [--FastMuon] [--phiRandom]
                        [--Cosmics COSMICS] [-n NEVENTS] [-i FIRSTEVENT] [-s THESEED] [-S SAMESEED] [-f INPUTFILE] [-g GEOFILE] [-o OUTPUTDIR] [-F] [-t] [--dry-run] [-D] [--debug {0,1,2}]

optional arguments:
  -h, --help            show this help message and exit
  --Pythia8             Use Pythia8
  --PG                  Use Particle Gun
  --pID PID             id of particle used by the gun (default=22)
  --Estart ESTART       start of energy range of particle gun (default=10 GeV)
  --Eend EEND           end of energy range of particle gun (default=10 GeV)
  --FixedTarget         Use FixedTarget generator
  --Ntuple              Use ntuple as input
  --MuonBack            Generate events from muon background file, --Cosmics=0 for cosmic generator data
  --FollowMuon          Make muonshield active to follow muons
  --FastMuon            Only transport muons for a fast muon only background estimate
  --phiRandom           only relevant for muon background generator, random phi
  --Cosmics COSMICS     Use cosmic generator, argument switch for cosmic generator 0 or 1
  -n NEVENTS, --nEvents NEVENTS
                        Number of events to generate
  -i FIRSTEVENT, --firstEvent FIRSTEVENT
                        First event of input file to use
  -s THESEED, --seed THESEED
                        Seed for random number. Only for experts, see TRrandom::SetSeed documentation
  -S SAMESEED, --sameSeed SAMESEED
                        can be set to an integer for the muonBackground simulation with specific seed for each muon, only for experts!
  -f INPUTFILE          Input file if not default file
  -g GEOFILE            geofile for muon shield geometry, for experts only
  -o OUTPUTDIR, --output OUTPUTDIR
                        Output directory
  -F                    default = False: copy only stable particles to stack, except for HNL events
  -t, --test            quick test
  --dry-run             stop after initialize
  -D, --display         store trajectories
  --debug {0,1,2}       1: print weights and field 2: make overlap check
  ```

Where the Pythia8 generator uses the charm cascade from the previous target production, and FixedTarget produces minimum bias Pythia8 events. The Particle Gun should also work.

The location and size of the tracking stations can be defined in lines 308ff of
run_simScript.py. For now, the target and other parts of the geometry are
defined in C++, partially configurable through `charmDet_conf.py`. This will
require further cleanup.
