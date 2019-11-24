#!/bin/bash

ulimit -s unlimited

for Wave in M2 S2 N2 K1 O1 ; do
FES=FES2014b
TIDDIR=$DDIR/$FES/Elevation_Currents/INTERP

SOSIE=$DEVGIT/sosie-3.0-base/bin/sosie3.x
CORVECT=$DEVGIT/sosie-3.0-base/bin/corr_vect.x

for zone in BAFBDY  GINBDY  HBDY  SBDY   ; do
  cd $zone
  for RI in real imag ; do
    case $RI in
    (real)  NEMORI=z1 ;;
    (imag)  NEMORI=z2 ;;
    esac
    cat ../namelist_elev.tmp  | sed -e "s/<ZONE>/$zone/" -e "s/<WAVE>/$Wave/g"  -e "s/<RI>/$RI/g" -e "s/<NEMORI>/$NEMORI/g" > namelist_sosie
    ln -sf $TIDDIR/${Wave}.${FES}_elevation_RI.nc  ./
    $SOSIE -f namelist_sosie
  done
  cd ../
done

done






exit

   clname(1)  =  'M2'
   clname(2)  =  'S2'
   clname(3)  =  'N2'
   clname(4)  =  'K1'
   clname(5)  =  'O1

M2.FES2014b_eastward_RI.nc  M2.FES2014b_elevation_RI.nc  M2.FES2014b_northward_RI.nc
BAFBDY  GINBDY  HBDY  SBDY

!!               Ex: you want all points where your field is <= 0 to become land mask,
cf_src     = '<WAVE>.FES2014b_elevation_RI.nc'
cv_src     = 'elevation_<RI>'
cf_x_src   = '<WAVE>.FES2014b_elevation_RI.nc'
ctarget    = '<TGT>'
cf_x_trg   = '<ZONE>36_domain_cfg_v3.3.2.nc'
cf_lsm_trg = '<ZONE>36_mesh_mask.nc'
cv_out    = '<WAVE>_<NEMORI>'

