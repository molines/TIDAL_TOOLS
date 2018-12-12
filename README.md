# TIDAL_TOOLS
A set of tools  linked with the use of tidal forcing in NEMO

## Description :
 This package is a set of fortran90 programs, useful in the processing of NEMO simulation when tides are involved. There
are 3 kinds of tools in the package, so far :
 1. Harmonic analysis for tidal constituents : **tid_harm_ana**
 1. Tidal prediction : **tid_predict**
 1. Management of tidal constituents :
  * convert from Amplitude/Phase to real/imaginary : **tid_conv_ri**
  * convert from real/imaginary to Amplitude/Phase : **tid_conv_ag**

 Although developped in the frame of NEMO, those tools can be easily used with any other kind of data, provided they use NetCdf format.

 The harmonic analysis and tidal prediction programs are just a rewriting of an old fortran package whose origin in uncertain. We improved the
readibility of the code, and used the NEMO coding rule to make it a NEMO companion. The user interface has also been improved via an extented use
of the namelist.

## Usage
### **tid_harm_ana**
     usage :  tid_harm_ana -l LST-files [-o HARM-file] [-n NAMLIST-file] [-nc4]
         ... [-zeromean]
       
      PURPOSE :
        Perform the harmonic analysis on the corresponding time series 
        represented by the list of files given as arguments.
       
      ARGUMENTS :
        -l LST-files : a blank separated list of input files
       
      OPTIONS :
        -n NAMLIST-file : Input namelist file. Default is 'namelist' 
        -o HARM-file : Name of the output file with analysed harmonic 
            constituents. Default is  res_harm.nc, or set in the namelist.
        -nc4 : output file is in Netcdf4/Hdf5 with chunking and deflation
        -zeromean : subtract spatial mean of the field before analysis.
       
      REQUIRED FILES :
         If -zeromean option, mesh_hgr and mask files are required.
      OUTPUT : 
        netcdf file : res_harm.nc

### **tid_predict**
     usage :  tid_predict -n NAMLIST-file -f INPUT-file
       
      PURPOSE :
        Perform tide prediction using a set of harmonic constituents.
        The prediction will correspond to the time period covered by
        the input file, on the same grid. 
        The output file will have the same netcdf format than the input file.
       
      ARGUMENTS :
        -n NAMLIST-file : give the name of the namelist to be used.
        -f INPUT-file   : give the name of the input file.
       
      OPTIONS :
         none so far ... 
        -nc4 : use netcdf4 with chinking and deflation 
       
      REQUIRED FILES :
          none 
       
      OUTPUT : 
        netcdf file : name specified in the namelist
          Variable : same as in the input file, specified in the namelist.

### **tid_conv_ag**

     usage : conv_ag <Tidal_File_RI> [VAR-rootname]
  
     PURPOSE:
       Compute amplitude and phase from the real/imaginary part
       tidal constituent given as input.
  
     ARGUMENTS:
       Tidal file with real/imaginary part (m) 
  
     OUTPUT:
       Netcdf file names <INPUT_FILE%_RI.nc>_AG.nc
       Variables : elevation_a, elevation_G

### **tid_conv_ri**

     usage : tid_conv_ri  <Tidal_File_AG> [VAR-rootname]

     PURPOSE:
       Compute real and imaginary part of the tides, in order to use SOSIE
       on continuous fields. An unlimited time axis is added in the output
       file for sosie.

     ARGUMENTS:
       Tidal file with amplitude and phase (degrees)
  
     OPTIONS:
       VAR-rootname :Root name of the variable to work with (default is 
       elevation )
       The program will look for <VAR-rootname>_a and <VAR-rootname>_G
  
     OUTPUT:
       Netcdf file names <INPUT_FILE%.nc>_<VAR-rootname>_RI.nc
       Variables : <VAR-rootname>_real, <VAR-rootname>_imag

### sample namelist used in **tid_harm_ana** and **tid_predict**
 A [sample namelist](../master/namelists/namelist_tideharm)  is provided in the namelist sub-directory.
 

