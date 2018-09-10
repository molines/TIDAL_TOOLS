PROGRAM tid_conv_ag
  !!======================================================================
  !!                     ***  PROGRAM  tid_conv_ag ***
  !!======================================================================
  !!  ** Purpose : convert a tidal constituent given as real/imaginary into
  !!           its equivalent amplitude and phase. Assume that the input file
  !!           is a NEMO-like file with lon and lat 2D and time axis
  !!
  !!  ** Method  : ampli = sqrt( real*real + imag*imag)
  !!               phase = atan2(imag, real) * 180./pi
  !!
  !! History :  1.0  : 03/2018  : J.M. Molines : 
  !!----------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE
  INTEGER(KIND=4)                           :: npiglo, npjglo, npt
  INTEGER(KIND=4)                           :: narg, ipos
  INTEGER(KIND=4)                           :: ncid, ierr, id, idx, idy, idt
  INTEGER(KIND=4)                           :: idlon, idlat, idampl, idphas, idtime

  REAL(KIND=4)                              :: spval

  REAL(KIND=8)                              :: dpi
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dreal, dimag
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dlon2d, dlat2d
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dampli, dphase
  REAL(KIND=8), DIMENSION(:)  , ALLOCATABLE :: dtime

  CHARACTER(LEN=255)  :: cf_in
  CHARACTER(LEN=255)  :: cf_out
! Default names 
  CHARACTER(LEN=80 )  :: cld_x='x'  ! i dimension name
  CHARACTER(LEN=80 )  :: cld_y='y'  ! j dimension name
  CHARACTER(LEN=80 )  :: cld_t='time_counter'  ! t dimension name (if any)

  CHARACTER(LEN=80 )  :: cv_lon='nav_lon'   ! longitude name
  CHARACTER(LEN=80 )  :: cv_lat='nav_lat'   ! longitude name
  CHARACTER(LEN=80 )  :: cv_tim='time_counter'        ! longitude name
  CHARACTER(LEN=80 )  :: cv_root='elevation'  ! variable rootname

  !!----------------------------------------------------------------------

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : conv_ag <Tidal_File_RI> [VAR-rootname]'
     PRINT *, ' '
     PRINT *, ' PURPOSE:'
     PRINT *, '    Compute amplitude and phase from the real/imaginary part'
     PRINT *, '    tidal constituent given as input.'
     PRINT *, ' '
     PRINT *, ' ARGUMENTS:'
     PRINT *, '    Tidal file with real/imaginary part (m) '
     PRINT *, ' '
     PRINT *, ' OUTPUT:'
     PRINT *, '    Netcdf file names <INPUT_FILE%_RI.nc>_AG.nc'
     PRINT *, '    Variables : elevation_a, elevation_G'
     PRINT *, ' '
     STOP 0
  ENDIF

  CALL getarg(1,cf_in)
  ipos=INDEX(cf_in,"_RI.nc") 
  cf_out=TRIM(cf_in(1:ipos-1))//'_AG.nc'
  PRINT *, ' Output file : ',TRIM(cf_out)

  ! read input file 
  ierr = NF90_OPEN(cf_in,NF90_NOWRITE, ncid)
  !  read dimension
  ierr = NF90_INQ_DIMID(ncid,cld_x,id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npiglo)
  ierr = NF90_INQ_DIMID(ncid,cld_y,id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npjglo)
  ierr = NF90_INQ_DIMID(ncid,cld_t,id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npt)
  !  Allocate arrays
  ALLOCATE ( dampli(npiglo, npjglo), dphase(npiglo,npjglo) )
  ALLOCATE ( dreal(npiglo, npjglo) , dimag(npiglo,npjglo)  )
  ALLOCATE ( dlon2d(npiglo, npjglo) , dlat2d(npiglo,npjglo)  )
  ALLOCATE ( dtime(1) )
  ! read lon/lat
  ierr = NF90_INQ_VARID(ncid,cv_lon, id )
    ierr = NF90_GET_VAR(ncid,id,dlon2d)
  ierr = NF90_INQ_VARID(ncid,cv_lat, id )
    ierr = NF90_GET_VAR(ncid,id,dlat2d)

  ! read arrays
  ierr = NF90_INQ_VARID(ncid,TRIM(cv_root)//'_real', id )
    ierr = NF90_GET_VAR(ncid,id,dreal)
  ierr = NF90_INQ_VARID(ncid,TRIM(cv_root)//'_imag', id )
    ierr = NF90_GET_VAR(ncid,id,dimag)
  ! get spval ( identical for both)
  ierr = NF90_GET_ATT(ncid,id,'_FillValue',spval)
  ! close input file
  ierr = NF90_CLOSE(ncid)

  ! COMPUTE real and imaginary part
  dpi=acos(-1.d0)
  WHERE(dreal /= spval)
     dampli(:,:)=sqrt( dreal*dreal + dimag*dimag )
  ELSEWHERE
     dampli(:,:)= spval*1.d0
  ENDWHERE

  WHERE(dreal /= spval)
     dphase(:,:)=ATAN2(dimag,dreal)*180.d0/dpi
  ELSEWHERE
     dphase(:,:)= spval*1.d0
  ENDWHERE

  ! Create output file
  ierr=NF90_CREATE(cf_out,OR(NF90_CLOBBER,NF90_64BIT_OFFSET),ncid)
  ! Add dimension
  ierr=NF90_DEF_DIM(ncid,cld_x,npiglo,idx)
  ierr=NF90_DEF_DIM(ncid,cld_y,npjglo,idy)
  ierr=NF90_DEF_DIM(ncid,cld_t,NF90_UNLIMITED,idt)
  ! add variables 
  ierr=NF90_DEF_VAR(ncid,cv_lon,NF90_DOUBLE,(/idx,idy/), idlon)
    ierr=NF90_PUT_ATT(ncid,idlon,'standard_name','longitude')
    ierr=NF90_PUT_ATT(ncid,idlon,'units','degrees')
    ierr=NF90_PUT_ATT(ncid,idlon,'axis','X')

  ierr=NF90_DEF_VAR(ncid,cv_lat,NF90_DOUBLE,(/idx,idy/), idlat)
    ierr=NF90_PUT_ATT(ncid,idlat,'standard_name','latitude')
    ierr=NF90_PUT_ATT(ncid,idlat,'units','degrees')
    ierr=NF90_PUT_ATT(ncid,idlat,'axis','Y')

  ierr=NF90_DEF_VAR(ncid,cv_tim,NF90_DOUBLE,(/idt/), idtime)

  ierr=NF90_DEF_VAR(ncid,TRIM(cv_root)//'_a',NF90_FLOAT,(/idx,idy,idt/), idampl)
    ierr=NF90_PUT_ATT(ncid,idampl,'coordinates','nav_lon nav_lat')
    ierr=NF90_PUT_ATT(ncid,idampl,'standard_name','Amplitude')
    ierr=NF90_PUT_ATT(ncid,idampl,'units','m')
    ierr=NF90_PUT_ATT(ncid,idampl,'_FillValue',spval)

  ierr=NF90_DEF_VAR(ncid,TRIM(cv_root)//'_G',NF90_FLOAT,(/idx,idy,idt/), idphas)
    ierr=NF90_PUT_ATT(ncid,idphas,'coordinates','nav_lon nav_lat')
    ierr=NF90_PUT_ATT(ncid,idphas,'standard_name','Phase')
    ierr=NF90_PUT_ATT(ncid,idphas,'units','degrees')
    ierr=NF90_PUT_ATT(ncid,idphas,'_FillValue',spval)

  ! end definition
  ierr = NF90_ENDDEF(ncid)

  ! Put variables
   dtime=0.d0
  ierr = NF90_PUT_VAR(ncid,idlon,dlon2d)
  ierr = NF90_PUT_VAR(ncid,idlat,dlat2d)
  ierr = NF90_PUT_VAR(ncid,idtime,dtime)
  ierr = NF90_PUT_VAR(ncid,idampl,dampli, start=(/1,1,1/), count=(/npiglo,npjglo,1/) )
  ierr = NF90_PUT_VAR(ncid,idphas,dphase, start=(/1,1,1/), count=(/npiglo,npjglo,1/) )

  ! close file
  ierr = NF90_CLOSE(ncid)
  
  
END PROGRAM tid_conv_ag
