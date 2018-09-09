PROGRAM tid_conv_ri
  !!======================================================================
  !!                     ***  PROGRAM  tid_conv_ri  ***
  !!======================================================================
  !!  ** Purpose :  convert a tidal constituent given as Amplitude, Phase into
  !!               its equivalent Real and Imaginary part
  !!  ** Method  : Real= Ampli*cos(pi/180*phase)
  !!           Imag= Ampli*sin(pi/180*phase)
  !!
  !! History :  1.0  : 03/2018  : J.M. Molines : 
  !!----------------------------------------------------------------------
  USE netcdf

  IMPLICIT NONE
  INTEGER(KIND=4)                           :: npiglo, npjglo
  INTEGER(KIND=4)                           :: narg, ipos
  INTEGER(KIND=4)                           :: ncid, ierr, id, idx, idy, idt
  INTEGER(KIND=4)                           :: idlon, idlat, idreal, idimag, idtim

  REAL(KIND=4)                              :: spval
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ampli, phase

  REAL(KIND=8)                              :: dpi
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dreal, dimag
  REAL(KIND=8), DIMENSION(:  ), ALLOCATABLE :: dlon, dlat

  CHARACTER(LEN=255)  :: cf_in
  CHARACTER(LEN=255)  :: cf_out
! Default names 
  CHARACTER(LEN=80 )  :: cld_x='nx'  ! i dimension name
  CHARACTER(LEN=80 )  :: cld_y='ny'  ! j dimension name
  CHARACTER(LEN=80 )  :: cld_t='nt'  ! j dimension name

  CHARACTER(LEN=80 )  :: cv_lon='longitude'   ! longitude name
  CHARACTER(LEN=80 )  :: cv_lat='latitude'    ! longitude name
  CHARACTER(LEN=80 )  :: cv_tim='time'        ! longitude name
  CHARACTER(LEN=80 )  :: cv_root='elevation'  ! variable rootname

  !!----------------------------------------------------------------------
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : tid_conv_ri  <Tidal_File_AG> [VAR-rootname]'
     PRINT *, ' '
     PRINT *, ' PURPOSE:'
     PRINT *, '    Compute real and imaginary part of the tides, in order to use SOSIE'
     PRINT *, '    on continuous fields. An unlimited time axis is added in the output'
     PRINT *,'     file for sosie.'
     PRINT *, ' '
     PRINT *, ' ARGUMENTS:'
     PRINT *, '    Tidal file with amplitude and phase (degrees)'
     PRINT *, ' '
     PRINT *, ' OPTIONS:'
     PRINT *, '    VAR-rootname :Root name of the variable to work with (default is '
     PRINT *, '  ',TRIM(cv_root)//' )'
     PRINT *, '      The program will look for <VAR-rootname>_a and <VAR-rootname>_G'
     PRINT *, ' '
     PRINT *, ' OUTPUT:'
     PRINT *, '    Netcdf file names <INPUT_FILE%.nc>_<VAR-rootname>_RI.nc'
     PRINT *, '    Variables : <VAR-rootname>_real, <VAR-rootname>_imag'
     PRINT *, ' '
     STOP 0
  ENDIF

  CALL getarg(1,cf_in)
  IF ( narg == 2 ) THEN
    CALL getarg( 2,cv_root )
  ENDIF
  ipos=INDEX(cf_in,".nc") 
  cf_out=TRIM(cf_in(1:ipos-1))//'_'//TRIM(cv_root)//'_RI.nc'
  PRINT *, ' Output file : ',TRIM(cf_out)

  ! read input file 
  ierr = NF90_OPEN(cf_in,NF90_NOWRITE, ncid)
  !  read dimension
  ierr = NF90_INQ_DIMID(ncid,cld_x,id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npiglo)
  ierr = NF90_INQ_DIMID(ncid,cld_y,id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npjglo)
  !  Allocate arrays
  ALLOCATE ( ampli(npiglo, npjglo), phase(npiglo,npjglo) )
  ALLOCATE ( dreal(npiglo, npjglo), dimag(npiglo,npjglo) )
  ALLOCATE ( dlon(npiglo), dlat(npjglo) )
  ! read lon/lat
  ierr = NF90_INQ_VARID(ncid,cv_lon, id )
    ierr = NF90_GET_VAR(ncid,id,dlon)
  ierr = NF90_INQ_VARID(ncid,cv_lat, id )
    ierr = NF90_GET_VAR(ncid,id,dlat)

  ! read arrays
  ierr = NF90_INQ_VARID(ncid,TRIM(cv_root)//'_a', id )
    ierr = NF90_GET_VAR(ncid,id,ampli)
  ierr = NF90_INQ_VARID(ncid,TRIM(cv_root)//'_G', id )
    ierr = NF90_GET_VAR(ncid,id,phase)
  ! get spval ( identical for both)
  ierr = NF90_GET_ATT(ncid,id,'_FillValue',spval)
  ! close input file
  ierr = NF90_CLOSE(ncid)

  ! COMPUTE real and imaginary part
  dpi=acos(-1.d0)
  WHERE(phase /= spval)
     dreal(:,:)=ampli(:,:)*cos(dpi/180.d0*phase(:,:))
  ELSEWHERE
     dreal(:,:)= spval*1.d0
  ENDWHERE

  WHERE(phase /= spval)
     dimag(:,:)=ampli(:,:)*sin(dpi/180.d0*phase(:,:))
  ELSEWHERE
     dimag(:,:)= spval*1.d0
  ENDWHERE

  ! Create output file
  ierr=NF90_CREATE(cf_out,OR(NF90_NOCLOBBER,NF90_64BIT_OFFSET),ncid)
  ! Add dimension
  ierr=NF90_DEF_DIM(ncid,cld_x,npiglo,idx)
  ierr=NF90_DEF_DIM(ncid,cld_y,npjglo,idy)
  ierr=NF90_DEF_DIM(ncid,cld_t,NF90_UNLIMITED,idt)
  ! add variables 
  ierr=NF90_DEF_VAR(ncid,cv_lon,NF90_DOUBLE,(/idx/), idlon)
    ierr=NF90_PUT_ATT(ncid,idlon,'standard_name','longitude')
    ierr=NF90_PUT_ATT(ncid,idlon,'units','degrees')
    ierr=NF90_PUT_ATT(ncid,idlon,'axis','X')

  ierr=NF90_DEF_VAR(ncid,cv_lat,NF90_DOUBLE,(/idy/), idlat)
    ierr=NF90_PUT_ATT(ncid,idlat,'standard_name','latitude')
    ierr=NF90_PUT_ATT(ncid,idlat,'units','degrees')
    ierr=NF90_PUT_ATT(ncid,idlat,'axis','Y')

  ierr=NF90_DEF_VAR(ncid,cv_tim,NF90_DOUBLE,(/idt/), idtim)

  ierr=NF90_DEF_VAR(ncid,TRIM(cv_root)//'_real',NF90_DOUBLE,(/idx,idy,idt/), idreal)
    ierr=NF90_PUT_ATT(ncid,idreal,'coordinates','latitude longitude')
    ierr=NF90_PUT_ATT(ncid,idreal,'standard_name','tidal_real_part')
    ierr=NF90_PUT_ATT(ncid,idreal,'units','m')
    ierr=NF90_PUT_ATT(ncid,idreal,'_FillValue',spval*1.d0)

  ierr=NF90_DEF_VAR(ncid,TRIM(cv_root)//'_imag',NF90_DOUBLE,(/idx,idy,idt/), idimag)
    ierr=NF90_PUT_ATT(ncid,idimag,'coordinates','latitude longitude')
    ierr=NF90_PUT_ATT(ncid,idimag,'standard_name','tidal_imag_part')
    ierr=NF90_PUT_ATT(ncid,idimag,'units','m')
    ierr=NF90_PUT_ATT(ncid,idimag,'_FillValue',spval*1.d0)

  ! end definition
  ierr = NF90_ENDDEF(ncid)

  ! Put variables
  ierr = NF90_PUT_VAR(ncid,idlon,dlon)
  ierr = NF90_PUT_VAR(ncid,idlat,dlat)
  ierr = NF90_PUT_VAR(ncid,idtim,(/0.d0/) )
  ierr = NF90_PUT_VAR(ncid,idreal,dreal)
  ierr = NF90_PUT_VAR(ncid,idimag,dimag)

  ! close file
  ierr = NF90_CLOSE(ncid)
  
  
END PROGRAM tid_conv_ri
