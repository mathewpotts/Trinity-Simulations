# Trinity-Simulations
#### Authors: 
Andrew Wang, Otte Nepomuk, Mathew Potts

## Processed CORSIKA ROOT Files (Processed_CORSIKA_ROOT_Files)
Authors: Andrew Wang, Mathew Potts

## Performance Simulations (TrinityPerformaceSims)
#### Author: 
Otte Nepomuk

#### Description:

#### Functions:


## Photo-Electron Simulations (Trinity_PE_Sims)
#### Author: 
Andrew Wang, Mathew Potts

#### Description: 
This directory consists of a single C++ script that uses parameterization to simulate the Trinity camera's response to Extensive Air Showers (EAS's) caused by upward going Taus.  

#### Useage:
The current usage of this script requires you to load the script into ROOT and executing functions via the terminal. If you aren't using an interactive desktop, X11 forwarding, or other such method of graphic deplay you will need to modify the script to save the histograms. Any changes that are made requries you to re-make the file. 
```bash
make clean
make
```

#### Functions:
##### Double_t myPEfunction2(Double_t azi, Double_t elv, Double_t l)
##### Double_t DistanceThroughEarth(Double_t azimuth, Double_t elevation, Double_t y)
##### Double_t GH(Double_t *x, Double_t *p) 
##### Double_t Gauss2D(Double_t *x, Double_t *p)
##### Double_t DistanceFovBelow(Double_t fovbelow, Double_t hdist)
##### Double_t *GetTauEarthCoords(Double_t azi, Double_t elv, Double_t l, Double_t hdist)
##### Double_t *TauTrajectoryFoV(Double_t *coords, Double_t azi, Double_t elv, Double_t fov, Double_t hdist, bool lower)
##### bool TauEmergeBelowVFoV(Double_t *coords, Double_t vfovbelow, Double_t hdist)
##### bool AirShowerAppearsInVFoV(Double_t *coords, Double_t dprime, Double_t d, Double_t azi, Double_t elv, Double_t l, Double_t vfovabove, Double_t vfovbelow, Double_t hdist, Double_t s, short tau_config)
##### bool *ASCutoff(Double_t *base, Double_t *tip, Double_t vfovbelow, Double_t vfovabove, Double_t hdist)
##### void SimPEinCam(Double_t azi, Double_t elv, Double_t l, Double_t E)
##### void SimNAirShowers(Double_t azi, Double_t elv, Double_t l, Double_t E, int trials)
