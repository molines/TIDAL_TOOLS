PROGRAM tid_conv_ag
  !!======================================================================
  !!                     ***  PROGRAM  tid_conv_ag ***
  !!======================================================================
  !!  ** Purpose : convert a tidal constituent given as real/imaginary into
  !!           its equivalent amplitude and phase.
  !!
  !!  ** Method  : ampli = sqrt( real*real + imag*imag)
  !!               phase = atan2(imag, real) * 180./pi
  !!
  !! History :  1.0  : 03/2018  : J.M. Molines : 
  !!----------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE
  INTEGER(KIND=4)                           :: npiglo, npjglo
  INTEGER(KIND=4)                           :: narg, ipos
  INTEGER(KIND=4)                           :: ncid, ierr, id, idx, idy
  INTEGER(KIND=4)                           :: idlon, idlat, idampl, idphas

  REAL(KIND=4)                              :: spval

  REAL(KIND=8)                              :: dpi
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dreal, dimag
  REAL(KIND=8), DIMENSION(:  ), ALLOCATABLE :: dlon, dlat
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dampli, dphase

  CHARACTER(LEN=255)  :: cf_in
  CHARACTER(LEN=255)  :: cf_out

  !!----------------------------------------------------------------------

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : conv_ag <Tidal_File_RI>'
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
  ierr = NF90_INQ_DIMID(ncid,'nx',id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npiglo)
  ierr = NF90_INQ_DIMID(ncid,'ny',id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npjglo)
  !  Allocate arrays
  ALLOCATE ( dampli(npiglo, npjglo), dphase(npiglo,npjglo) )
  ALLOCATE ( dreal(npiglo, npjglo) , dimag(npiglo,npjglo)  )
  ALLOCATE ( dlon(npiglo), dlat(npjglo) )
  ! read lon/lat
  ierr = NF90_INQ_VARID(ncid,'longitude', id )
    ierr = NF90_GET_VAR(ncid,id,dlon)
  ierr = NF90_INQ_VARID(ncid,'latitude', id )
    ierr = NF90_GET_VAR(ncid,id,dlat)

  ! read arrays
  ierr = NF90_INQ_VARID(ncid,'tid_real', id )
    ierr = NF90_GET_VAR(ncid,id,dreal)
  ierr = NF90_INQ_VARID(ncid,'tid_imag', id )
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
     dphase(:,:)=ATAN2(dimag,dreal)
  ELSEWHERE
     dphase(:,:)= spval*1.d0
  ENDWHERE

  ! Create output file
  ierr=NF90_CREATE(cf_out,OR(NF90_CLOBBER,NF90_64BIT_OFFSET),ncid)
  ! Add dimension
  ierr=NF90_DEF_DIM(ncid,'nx',npiglo,idx)
  ierr=NF90_DEF_DIM(ncid,'ny',npjglo,idy)
  ! add variables 
  ierr=NF90_DEF_VAR(ncid,'longitude',NF90_DOUBLE,(/idx/), idlon)
    ierr=NF90_PUT_ATT(ncid,idlon,'standard_name','longitude')
    ierr=NF90_PUT_ATT(ncid,idlon,'units','degrees')
    ierr=NF90_PUT_ATT(ncid,idlon,'axis','X')

  ierr=NF90_DEF_VAR(ncid,'latitude',NF90_DOUBLE,(/idy/), idlat)
    ierr=NF90_PUT_ATT(ncid,idlat,'standard_name','latitude')
    ierr=NF90_PUT_ATT(ncid,idlat,'units','degrees')
    ierr=NF90_PUT_ATT(ncid,idlat,'axis','Y')

  ierr=NF90_DEF_VAR(ncid,'elevation_a',NF90_FLOAT,(/idx,idy/), idampl)
    ierr=NF90_PUT_ATT(ncid,idampl,'coordinates','latitude longitude')
    ierr=NF90_PUT_ATT(ncid,idampl,'standard_name','Amplitude')
    ierr=NF90_PUT_ATT(ncid,idampl,'units','m')
    ierr=NF90_PUT_ATT(ncid,idampl,'_FillValue',spval)

  ierr=NF90_DEF_VAR(ncid,'elevation_G',NF90_FLOAT,(/idx,idy/), idphas)
    ierr=NF90_PUT_ATT(ncid,idphas,'coordinates','latitude longitude')
    ierr=NF90_PUT_ATT(ncid,idphas,'standard_name','Phase')
    ierr=NF90_PUT_ATT(ncid,idphas,'units','degrees')
    ierr=NF90_PUT_ATT(ncid,idphas,'_FillValue',spval)

  ! end definition
  ierr = NF90_ENDDEF(ncid)

  ! Put variables
  ierr = NF90_PUT_VAR(ncid,idlon,dlon)
  ierr = NF90_PUT_VAR(ncid,idlat,dlat)
  ierr = NF90_PUT_VAR(ncid,idampl,dampli)
  ierr = NF90_PUT_VAR(ncid,idphas,dphase)

  ! close file
  ierr = NF90_CLOSE(ncid)
  
  
END PROGRAM tid_conv_ag
