PROGRAM tid_conv_ri
  !! -------------------------------------------------------
  !!      PROGRAM tid_conv_ri
  !!      
  !!  Purpose : convert a tidal consituant giben as Amplitude, Phase into
  !!           its equivalent Real and Imaginary part
  !!
  !!  Method : Real= Ampli*cos(pi/180*phase)
  !!           Imag= Amplu*sin(pi/180*phase)
  !! ---------------------------------------------------------
  USE netcdf
  IMPLICIT NONE
  INTEGER(KIND=4)                           :: npiglo, npjglo
  INTEGER(KIND=4)                           :: narg, ipos
  INTEGER(KIND=4)                           :: ncid, ierr, id, idx, idy
  INTEGER(KIND=4)                           :: idlon, idlat, idreal, idimag

  REAL(KIND=4)                              :: spval
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: ampli, phase

  REAL(KIND=8)                              :: dpi
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dreal, dimag
  REAL(KIND=8), DIMENSION(:  ), ALLOCATABLE :: dlon, dlat

  CHARACTER(LEN=255)  :: cf_in
  CHARACTER(LEN=255)  :: cf_out

  !! ---------------------------------------------------------------
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage : conv_RI <Tidal_File_AG>'
     PRINT *, ' '
     PRINT *, ' PURPOSE:'
     PRINT *, '    Compute real and imaginary part of the tides, in order to use SOSIE'
     PRINT *, '    on continuous fields.'
     PRINT *, ' '
     PRINT *, ' ARGUMENTS:'
     PRINT *, '    Tidal file with amplitude and phase (degrees)'
     PRINT *, ' '
     PRINT *, ' OUTPUT:'
     PRINT *, '    Netcdf file names <INPUT_FILE%.nc>_RI.nc'
     PRINT *, '    Variables : tid_real, tid_imag'
     PRINT *, ' '
     STOP 0
  ENDIF

  CALL getarg(1,cf_in)
  ipos=INDEX(cf_in,".nc") 
  cf_out=TRIM(cf_in(1:ipos-1))//'_RI.nc'
  PRINT *, ' Output file : ',TRIM(cf_out)

  ! read input file 
  ierr = NF90_OPEN(cf_in,NF90_NOWRITE, ncid)
  !  read dimension
  ierr = NF90_INQ_DIMID(ncid,'nx',id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npiglo)
  ierr = NF90_INQ_DIMID(ncid,'ny',id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npjglo)
  !  Allocate arrays
  ALLOCATE ( ampli(npiglo, npjglo), phase(npiglo,npjglo) )
  ALLOCATE ( dreal(npiglo, npjglo), dimag(npiglo,npjglo) )
  ALLOCATE ( dlon(npiglo), dlat(npjglo) )
  ! read lon/lat
  ierr = NF90_INQ_VARID(ncid,'longitude', id )
    ierr = NF90_GET_VAR(ncid,id,dlon)
  ierr = NF90_INQ_VARID(ncid,'latitude', id )
    ierr = NF90_GET_VAR(ncid,id,dlat)

  ! read arrays
  ierr = NF90_INQ_VARID(ncid,'elevation_a', id )
    ierr = NF90_GET_VAR(ncid,id,ampli)
  ierr = NF90_INQ_VARID(ncid,'elevation_G', id )
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

  ierr=NF90_DEF_VAR(ncid,'tid_real',NF90_DOUBLE,(/idx,idy/), idreal)
    ierr=NF90_PUT_ATT(ncid,idreal,'coordinates','latitude longitude')
    ierr=NF90_PUT_ATT(ncid,idreal,'standard_name','tidal_real_part')
    ierr=NF90_PUT_ATT(ncid,idreal,'units','m')
    ierr=NF90_PUT_ATT(ncid,idreal,'_FillValue',spval*1.d0)

  ierr=NF90_DEF_VAR(ncid,'tid_imag',NF90_DOUBLE,(/idx,idy/), idimag)
    ierr=NF90_PUT_ATT(ncid,idimag,'coordinates','latitude longitude')
    ierr=NF90_PUT_ATT(ncid,idimag,'standard_name','tidal_imag_part')
    ierr=NF90_PUT_ATT(ncid,idimag,'units','m')
    ierr=NF90_PUT_ATT(ncid,idimag,'_FillValue',spval*1.d0)

  ! end definition
  ierr = NF90_ENDDEF(ncid)

  ! Put variables
  ierr = NF90_PUT_VAR(ncid,idlon,dlon)
  ierr = NF90_PUT_VAR(ncid,idlat,dlat)
  ierr = NF90_PUT_VAR(ncid,idreal,dreal)
  ierr = NF90_PUT_VAR(ncid,idimag,dimag)

  ! close file
  ierr = NF90_CLOSE(ncid)
  
  
END PROGRAM tid_conv_ri
