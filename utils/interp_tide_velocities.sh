#!/bin/bash

ulimit -s unlimited

for Wave in M2 S2 N2 K1 O1 ; do
FES=FES2014b
TIDDIR=$DDIR/$FES/Elevation_Currents/INTERP

SOSIE=$DEVGIT/sosie-3.0-base/bin/sosie3.x
CORVECT=$DEVGIT/sosie-3.0-base/bin/corr_vect.x

for zone in BAFBDY  GINBDY  HBDY  SBDY   ; do
  cd $zone

# eastward
  for RI in real imag ; do
    case $RI in
    (real)  NEMORI=u1 ;;
    (imag)  NEMORI=u2 ;;
    esac
    cat ../namelist_eastward.tmp  | sed -e "s/<ZONE>/$zone/" -e "s/<WAVE>/$Wave/g"  -e "s/<RI>/$RI/g" -e "s/<NEMORI>/$NEMORI/g" > ${Wave}_${RI}.namelist_x
    ln -sf $TIDDIR/${Wave}.${FES}_eastward_RI.nc  ./
    $SOSIE -f ${Wave}_${RI}.namelist_x
  done

# northward
  for RI in real imag ; do
    case $RI in
    (real)  NEMORI=v1 ;;
    (imag)  NEMORI=v2 ;;
    esac
    cat ../namelist_northward.tmp  | sed -e "s/<ZONE>/$zone/" -e "s/<WAVE>/$Wave/g"  -e "s/<RI>/$RI/g" -e "s/<NEMORI>/$NEMORI/g" >  ${Wave}_${RI}.namelist_y
    ln -sf $TIDDIR/${Wave}.${FES}_northward_RI.nc  ./
    $SOSIE -f  ${Wave}_${RI}.namelist_y
  done

# correct orientation
  for RI in real imag ; do
  $CORVECT  -f ${Wave}_${RI}.namelist -G U -m ${zone}36_mesh_mask.nc 
  done

# clean in zone
  rm -f *raw*nc
  rm -f *namelist* 
  rm -f ${Wave}.${FES}_*.nc

  cd ../

done

done






exit
BAFBDY36_domain_cfg_v3.3.2.nc   M2.FES2014b_elevation_RI.nc  M2_z1_FES2014b-BAFBDY36_tid.nc  O1_z1_FES2014b-BAFBDY36_tid.nc  uraw_FES2014b-BAFBDY36_M2_u1.nc
BAFBDY36_mesh_mask.nc           M2.FES2014b_northward_RI.nc  M2_z2_FES2014b-BAFBDY36_tid.nc  O1_z2_FES2014b-BAFBDY36_tid.nc  uraw_FES2014b-BAFBDY36_M2_u2.nc
K1.FES2014b_elevation_RI.nc     M2_imag.namelist_x           N2.FES2014b_elevation_RI.nc     S2.FES2014b_elevation_RI.nc     vraw_FES2014b-BAFBDY36_M2_v1.nc
K1_z1_FES2014b-BAFBDY36_tid.nc  M2_imag.namelist_y           N2_z1_FES2014b-BAFBDY36_tid.nc  S2_z1_FES2014b-BAFBDY36_tid.nc  vraw_FES2014b-BAFBDY36_M2_v2.nc
K1_z2_FES2014b-BAFBDY36_tid.nc  M2_real.namelist_x           N2_z2_FES2014b-BAFBDY36_tid.nc  S2_z2_FES2014b-BAFBDY36_tid.nc
M2.FES2014b_eastward_RI.nc      M2_real.namelist_y           O1.FES2014b_elevation_RI.nc     namelist_sosie
[rcli002@jean-zay1: BAFBDY]$ ls
BAFBDY36_domain_cfg_v3.3.2.nc   M2.FES2014b_elevation_RI.nc  M2_z1_FES2014b-BAFBDY36_tid.nc  O1_z1_FES2014b-BAFBDY36_tid.nc  uraw_FES2014b-BAFBDY36_M2_u1.nc
BAFBDY36_mesh_mask.nc           M2.FES2014b_northward_RI.nc  M2_z2_FES2014b-BAFBDY36_tid.nc  O1_z2_FES2014b-BAFBDY36_tid.nc  uraw_FES2014b-BAFBDY36_M2_u2.nc
K1.FES2014b_elevation_RI.nc     M2_imag.namelist_x           N2.FES2014b_elevation_RI.nc     S2.FES2014b_elevation_RI.nc     vraw_FES2014b-BAFBDY36_M2_v1.nc
K1_z1_FES2014b-BAFBDY36_tid.nc  M2_imag.namelist_y           N2_z1_FES2014b-BAFBDY36_tid.nc  S2_z1_FES2014b-BAFBDY36_tid.nc  vraw_FES2014b-BAFBDY36_M2_v2.nc
K1_z2_FES2014b-BAFBDY36_tid.nc  M2_real.namelist_x           N2_z2_FES2014b-BAFBDY36_tid.nc  S2_z2_FES2014b-BAFBDY36_tid.nc
M2.FES2014b_eastward_RI.nc      M2_real.namelist_y           O1.FES2014b_elevation_RI.nc     namelist_sosie


M2.FES2014b_eastward_RI.nc  M2.FES2014b_elevation_RI.nc  M2.FES2014b_northward_RI.nc
