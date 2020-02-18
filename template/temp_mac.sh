#!/bin/sh

#define parameters which are passed in.
particle=$1
energy=$2
number=$3

#define the template.
cat  << EOF
#/tracking/verbose 1

# Macro file for the
# initialization of G4GeneralParticleSource

# cosmics
/gps/energy $energy GeV
#/gps/particle mu-
#/gps/particle pi-

#/gps/particle gamma
/gps/particle $particle

/gps/pos/type Plane
/gps/pos/shape Square


/gps/pos/halfx 0.1 mm
/gps/pos/halfy 0.1 mm
/gps/pos/centre 0 0 -1900 mm
/gps/direction 0 0 1



#beam divergence
#/gps/ang/type beam2d
#/gps/pos/centre 0 0 -1.9 mm    # start beam close to crystal to have 0 time close to 0

#/gps/ang/rot1 -1. 0. 0.    # parallel beam --> eta = 0
#/gps/ang/rot1 -1. 0. 0.521    # 27.52 deg angle beam --> eta = 0.5
#/gps/ang/rot1 -1. 0. 1.175    # 49.6  deg angle beam --> eta = 1
#/gps/ang/rot1 -1. 0. 2.083    # 64.35 deg angle beam --> eta = 1.48

# parallel beam in opposite direction (SiPM will be in the front)
#/gps/pos/centre 0 0 35 mm
#/gps/ang/rot1 1. 0. 0.    

#/gps/ang/rot1 -1. 0. 0.018    # 1 deg angle beam
#/gps/ang/rot1 -1. 0. 0.036    # 2 deg angle beam
#/gps/ang/rot2 0. 0. -1.
#/gps/ang/sigma_x 2 deg
#/gps/ang/sigma_y 2 deg


# radioactive source
# 511 for Sodium, 662 for Cesium, 160 for Co57, 1250 for Co60
#/gps/energy 160 keV    
#/gps/particle gamma
#/gps/ang/type iso


/run/beamOn $number
EOF
