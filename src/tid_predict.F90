PROGRAM tid_predict
  !!======================================================================
  !!                     ***  PROGRAM  tid_predict  ***
  !!=====================================================================
  !!  ** Purpose : Perform a tidal prediction from cotidal maps
  !!
  !!  ** Method  : Just apply tidal prediction formulae ..
  !!
  !! History :  1.0  : 12/2018  : J.M. Molines  from Mercator .. from ...
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  !!
  !!----------------------------------------------------------------------
  USE tide
  USE utils
  USE netcdf

  IMPLICIT NONE

  INTEGER(KIND=4)               :: npconst
  INTEGER(KIND=4)               :: jt, jn, ji, jj
  INTEGER(KIND=4)               :: narg, iargc, ijarg
  INTEGER(KIND=4)               :: npi   , npj, npt
  INTEGER(KIND=4)               :: npiout, npjout
  INTEGER(KIND=4)               :: npih  , npjh
  INTEGER(KIND=4)               :: ncid, istatus, id_var, id_vart, id_varx, id_vary
  INTEGER(KIND=4)               :: ii0, ij0, inum = 20
  INTEGER(KIND=4)               :: iyy, imm, idd, ihh, imn, isec
  INTEGER(KIND=4)               :: ijul0, ijuli, it0
  INTEGER(KIND=4), DIMENSION(3) :: id_dim

  REAL(KIND=8)                                :: dl_amp, dl_ph, dl_hfrac, dcoef
  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: domega, dft, dut, dvt, dvt0
  REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: dg_time
  REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dlon, dlat
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dres
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dtempx, dtempy

  CHARACTER(LEN=100), DIMENSION(jpmax_harmo) :: cname
  CHARACTER(LEN=120)                         :: cldum
  CHARACTER(LEN=120)                         :: cf_in
  CHARACTER(LEN=120)                         :: cf_namli
  CHARACTER(LEN=100)                         :: ca_units     
  CHARACTER(LEN=120)                         :: cn_fharm, cn_fout  
  CHARACTER(LEN=100)                         :: cn_var_nm
  CHARACTER(LEN=255)                         :: cn_var_ln, cn_var_unit
  CHARACTER(LEN=120)                         :: cn_dim_x, cn_dim_y, cn_dim_t
  CHARACTER(LEN=120)                         :: cn_var_lon,cn_var_lat, cn_var_time
  CHARACTER(LEN=120)                         :: cn_att_miss
  CHARACTER(LEN=120)                         :: cn_dim_x_harm, cn_dim_y_harm
  CHARACTER(LEN=120)                         :: cn_var_lon_harm,cn_var_lat_harm
  CHARACTER(LEN=120)                         :: cn_att_miss_harm
  CHARACTER(LEN=120)                         :: cn_real_ext, cn_imag_ext

  LOGICAL :: ln_lonlat_2d_in=.TRUE.
  LOGICAL :: ln_lonlat_2d_harm=.TRUE.
  LOGICAL :: l_moor =.TRUE.
  LOGICAL :: lnc4   =.FALSE.
  LOGICAL :: l_exist, lchk=.FALSE.

  NAMELIST /prediction_param/ cn_var_nm, cn_var_ln, cn_var_unit,  cn_fharm, cn_fout
  NAMELIST /predic_constituents/ cname
  NAMELIST /input_file_format/ cn_dim_x, cn_dim_y, cn_dim_t,        &
       &                       cn_var_lon, cn_var_lat, cn_var_time, &
       &                       cn_att_miss, ln_lonlat_2d_in

   NAMELIST /harm_file_format/ cn_dim_x_harm,   cn_dim_y_harm,      &
        &                      cn_var_lon_harm, cn_var_lat_harm,    &
        &                      cn_att_miss_harm, ln_lonlat_2d_harm, &
        &                      cn_real_ext, cn_imag_ext

  !!----------------------------------------------------------------------
  !! TIDAL_TOOLS , MEOM 2018
  !! $Id$
  !! Copyright (c) 2018, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/TIDAL_TOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  ! constituents must be chosen in the namelist
  cname(:)='none'

  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  tid_predict -n NAMLIST-file -f INPUT-file'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Perform tide prediction using a set of harmonic constituents.' 
     PRINT *,'       The prediction will correspond to the time period covered by'
     PRINT *,'       the input file, on the same grid.  A residual file (de-tided'
     PRINT *,'       signal) will be produced.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -n NAMLIST-file : give the name of the namelist to be used.' 
     PRINT *,'       -f INPUT-file : give the name of the input file.'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        none so far ... '
     PRINT *,'       -nc4 : use netcdf4 with chinking and deflation '
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'         none '
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file '
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       tid_harm_ana'
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-n'   ) ; CALL getarg(ijarg, cf_namli ) ; ijarg=ijarg+1
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_in    ) ; ijarg=ijarg+1
        ! option
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  ! sanity check
  lchk = lchk .OR.  chkfile(cf_namli)
  lchk = lchk .OR.  chkfile(cf_in   )
  IF ( lchk ) STOP   ! missing input file

  ! 0. Read the namelist
  ! -------------------------------

  OPEN( inum, file = cf_namli,status = 'old' ,  form   = 'formatted' )
  READ( inum, NML = prediction_param )

  REWIND(inum)
  READ( inum, NML = predic_constituents )

  ! REWIND(inum)
  ! READ( inum, NML = input_file_format )

  ! REWIND(inum)
  ! READ( inum, NML = output_file_format )
  CLOSE( inum )
  lchk = lchk .OR. chkfile(cn_fharm)
  IF ( lchk ) STOP   ! missing input file


  ! look for the number of required constituents for analysis
  npconst=COUNT ( (cname(:) /= 'none') )
  ALLOCATE ( domega(npconst), dft(npconst), dut(npconst), dvt(npconst), dvt0(npconst) )

  PRINT *, ' '
  PRINT *, 'Variable to be predicted = ', TRIM(cn_var_nm)
  PRINT *, 'Variable long name       = ', TRIM(cn_var_ln)
  PRINT *, 'Variable units           = ', TRIM(cn_var_unit)
  PRINT *, 'Harmonic data file       = ', TRIM(cn_fharm)
  PRINT *, 'Prediction output file   = ', TRIM(cn_fout)
  PRINT *, ' '

  ! Read input file
  !----------------- 
  CALL ReadFile

  WRITE(6,*)'INPUT FILE START DATE (ddmmyy) ', idd, imm, iyy
  WRITE(6,*) 'Number OF time dumps ', npt
  WRITE(6,*) 'Number OF harmonics  ', npconst
  IF (l_moor) THEN
     WRITE(6,*)'Mooring file indices ', ii0, ij0
  ENDIF
  CALL FLUSH(6)


  ! Set tidal waves pulsations and nodal corrections (at initial time)
  !-------------------------------------------------------------------
  CALL tide_pulse( domega(1:npconst), cname(1:npconst) ,npconst)
  CALL tide_vuf( dvt(1:npconst), dut(1:npconst) , dft(1:npconst), cname(1:npconst) ,npconst, iyy, imm, idd, dl_hfrac)
  dvt0 =dvt

  ! Read harmonic data file
  !-------------------------      
  CALL ReadHarmonics

  ! Build residual
  !---------------

  it0 = FLOOR(dg_time(1)/dcoef)*86400.
  CALL CreateOutput (ncid) 

  DO jt=1,npt 

     ijul0 = ijuli + FLOOR(dg_time(jt)/dcoef) ! Julian day of first time dump in file
     CALL caldat(ijul0,imm,idd,iyy)
     dl_hfrac = ijuli+dg_time(jt)/dcoef-ijul0
     !
     !          write(6,*)
     !         write(6,*)' time ',dg_time(jt)
     !         write(6,*)' ijuli ijul0 dl_hfrac ',ijuli,ijul0,dl_hfrac
     !         write(6,*) imm,idd,iyy
     !         write(6,*)
     !
     CALL tide_vuf( dvt(1:npconst), dut(1:npconst) , dft(1:npconst), cname(1:npconst) ,npconst, iyy, imm, idd, dl_hfrac)
     dres(:,:,:) = 0.d0
     DO jn=1, npconst
        DO ji=1,npiout
           DO jj=1,npjout
              dl_amp = SQRT(dtempx(ji,jj,jn)**2 + dtempy(ji,jj,jn)**2)
              dl_ph = ATAN2(dtempy(ji,jj,jn),dtempx(ji,jj,jn)) + dvt0(jn) + dut(jn)
              dres(ji,jj,1) = dres(ji,jj,1) + dft(jn) * dl_amp * COS(dl_ph) * COS(domega(jn)*(dg_time(jt)-it0)) &
                   &                      - dft(jn) * dl_amp * SIN(dl_ph) * SIN(domega(jn)*(dg_time(jt)-it0))
           END DO
        END DO
        WRITE(6,'(a,1x,4(f18.8,2x))') TRIM(cname(jn)),dft(jn),dl_amp,SIN(dl_ph),dg_time(jt)-it0
     END DO
     istatus = NF90_PUT_VAR(ncid, id_var,dres,start = (/ 1, 1, jt /),count = (/ npiout, npjout, 1 /))
  END DO

  istatus = NF90_CLOSE(ncid)
CONTAINS

  SUBROUTINE ReadFile
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ReadFile  ***
    !!
    !! ** Purpose :  read basic information from the input file (time, map ) 
    !!
    !! ** Method  :  read file assuming some format defined in the namelist 
    !!
    !!----------------------------------------------------------------------
    ! local Variables :
    INTEGER(KIND=4) :: istatus
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dl_lon, dl_lat
    !!----------------------------------------------------------------------
    istatus = NF90_OPEN (cf_in, NF90_NOWRITE, ncid)
    ! check for moring file 
    istatus = NF90_INQ_DIMID (ncid,cn_dim_x,id_var)
    IF (istatus == NF90_NOERR) THEN 
       l_moor=.FALSE.
    ELSE
       l_moor=.TRUE.
    ENDIF

    ! Get mooring indice on global grid
    IF (l_moor) THEN
       istatus = NF90_GET_ATT(ncid, NF90_GLOBAL, "i-indice", ii0)
       IF (istatus /= NF90_NOERR) THEN
          PRINT *, 'Can not find mooring i-indice global attribute in: ', cf_in
          !              STOP
       ENDIF
       istatus = NF90_GET_ATT(ncid, NF90_GLOBAL, "j-indice", ij0)
       IF (istatus /= NF90_NOERR) THEN
          PRINT *, 'Can not find mooring j-indice global attribute in: ', cf_in
          !              STOP
       ENDIF
    ENDIF

    istatus=NF90_INQ_DIMID (ncid,cn_dim_t,id_var)
    IF (istatus /= NF90_NOERR) THEN
       PRINT *, 'Pb with time-dimension in file: ', cf_in
       STOP
    ENDIF
    istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=npt)

    istatus = NF90_INQ_VARID(ncid,cn_var_time,id_var)
    IF (istatus /= NF90_NOERR) THEN
       PRINT *, 'Pb with time-dimension in file: ', cf_in
       STOP
    ENDIF

    ALLOCATE(dg_time(npt))
    istatus = NF90_GET_VAR(ncid,id_var, dg_time)
    istatus = NF90_GET_ATT(ncid, id_var, "units", ca_units)

    print *, ca_units
    IF ( index(ca_units,'hours') /= 0 ) THEN
       READ(ca_units,7000)  iyy, imm, idd, ihh, imn, isec
       dcoef=24.d0
    ELSEIF ( index(ca_units,'seconds') /= 0 ) THEN
       READ(ca_units,7001)  iyy, imm, idd, ihh, imn, isec
       dcoef=86400.d0
    ELSE
       PRINT *, 'time unit not understood'
       STOP
    ENDIF

7000 FORMAT('hours since ', I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)
7001 FORMAT('seconds since ', I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)

    ijuli = julday(imm,idd,iyy)             ! Julian day of time origin in data file
    ijul0 = ijuli + FLOOR(dg_time(1)/dcoef) ! Julian day of first time dump in file
    dl_hfrac=ihh+imn/60.+isec/3600.
    CALL caldat(ijul0,imm,idd,iyy)
    !
    WRITE(6,*)
    WRITE(6,*)' ijuli ijul0 dl_hfrac ',ijuli,ijul0,dl_hfrac
    WRITE(6,*) imm,idd,iyy
    WRITE(6,*)
    !
    IF (l_moor) THEN
       ALLOCATE(dlon(1,1), dlat(1,1))
       istatus = NF90_INQ_VARID(ncid, cn_var_lon, id_var) ; istatus = NF90_GET_VAR(ncid, id_var, dlon)
       istatus = NF90_INQ_VARID(ncid, cn_var_lat, id_var) ; istatus = NF90_GET_VAR(ncid, id_var, dlat)
    ELSE
       istatus=NF90_INQ_DIMID (ncid, cn_dim_x, id_var)  ; istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=npi)
       istatus=NF90_INQ_DIMID (ncid, cn_dim_y, id_var)  ; istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=npj)
       IF ( ln_lonlat_2d_in ) THEN
          ALLOCATE( dlon(npi,npj), dlat(npi,npj) )
          istatus = NF90_INQ_VARID(ncid, cn_var_lon, id_var) ; istatus = NF90_GET_VAR(ncid, id_var, dlon)
          istatus = NF90_INQ_VARID(ncid, cn_var_lat, id_var) ; istatus = NF90_GET_VAR(ncid, id_var, dlat)
       ELSE
          ALLOCATE( dlon(npi,1), dlat(1,npj) ,dl_lon(npi), dl_lat(npj))
          istatus = NF90_INQ_VARID(ncid, cn_var_lon, id_var) ; istatus = NF90_GET_VAR(ncid, id_var, dl_lon)
          istatus = NF90_INQ_VARID(ncid, cn_var_lat, id_var) ; istatus = NF90_GET_VAR(ncid, id_var, dl_lat)
          dlon(:,1)=dl_lon
          dlat(1,:)=dl_lat
          DEALLOCATE (dl_lon, dl_lat )
       ENDIF
    ENDIF
    istatus = NF90_CLOSE(ncid)

  END SUBROUTINE ReadFile


  SUBROUTINE ReadHarmonics
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ReadHarmonics ***
    !!
    !! ** Purpose :  Read the harmonic constituent file produced by a previous
    !!               analysis. 
    !!
    !! ** Method  :   
    !!
    !!----------------------------------------------------------------------
    ! local variables
    INTEGER(KIND=4) :: istatus
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dtemp2d
    !!----------------------------------------------------------------------

  istatus = NF90_OPEN (cn_fharm, NF90_NOWRITE, ncid)

  istatus=NF90_INQ_DIMID (ncid,cn_dim_x_harm,id_var)
  IF (istatus /= NF90_NOERR) THEN
     PRINT *, 'Pb with x-dimension in file: ', cn_fharm
     STOP
  ENDIF
  istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=npih)

  istatus = NF90_INQ_DIMID (ncid,cn_dim_y_harm,id_var)
  IF (istatus /= NF90_NOERR) THEN
     PRINT *, 'Pb with y-dimension in file: ', cn_fharm
     STOP
  ENDIF
  istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=npjh)

  ! Hopefully at this level, both cf_in and cf_harm grid corresponds
  ! JMM ==> sanity check ?

  ALLOCATE(dtemp2d(npih,npjh))
  IF (l_moor) THEN
     ALLOCATE(dtempx(1,1,npconst), dtempy(1,1,npconst))
  ELSE
     npi=npih
     npj=npjh
     ALLOCATE(dtempx(npi,npj,npconst), dtempy(npi,npj,npconst))
  ENDIF

  DO jt=1, npconst
     istatus = NF90_INQ_VARID( ncid,TRIM(cname(jt))//TRIM(cn_real_ext),id_var)
     IF (istatus /= NF90_NOERR) THEN
        PRINT *, 'Can not find harmonic data ', TRIM(cname(jt))//TRIM(cn_real_ext)//' in file: ', cn_fharm
        STOP
     ENDIF
     istatus = NF90_GET_VAR(ncid,id_var, dtemp2d)

     IF (l_moor) THEN 
        IF ((npih==1).AND.(npjh==1)) THEN   
           dtempx(1,1,jt)=dtemp2d(1,1)
        ELSE       
           dtempx(1,1,jt)=dtemp2d(ii0,ij0)
        ENDIF
     ELSE
        dtempx(:,:,jt)=dtemp2d(:,:)
     ENDIF

     istatus = NF90_INQ_VARID( ncid,TRIM(cname(jt))//TRIM(cn_imag_ext),id_var)
     IF (istatus /= NF90_NOERR) THEN
        PRINT *, 'Can not find harmonic data ', TRIM(cname(jt))//TRIM(cn_imag_ext)//' in file: ', cn_fharm
        STOP
     ENDIF
     istatus = NF90_GET_VAR(ncid,id_var, dtemp2d)
     IF (l_moor) THEN     
        IF ((npih==1).AND.(npjh==1)) THEN  
           dtempy(1,1,jt)=dtemp2d(1,1)
        ELSE 
           dtempy(1,1,jt)=dtemp2d(ii0,ij0)
        ENDIF
     ELSE
        dtempy(:,:,jt)=dtemp2d(:,:)
     ENDIF
  END DO
  istatus = NF90_CLOSE(ncid)
  DEALLOCATE(dtemp2d)

  END SUBROUTINE ReadHarmonics
  
  SUBROUTINE CreateOutput (kcid)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateOutput  ***
    !!
    !! ** Purpose :  Create output file and build header, then
    !!               return its ncid to the main program 
    !!               Assume the file format is identical to cf_in
    !!
    !! ** Method  :  Netcdf. Take care of lnc4 
    !!
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(out) :: kcid   ! netcdf handle for the output file
    ! local variables
    INTEGER(KIND=4) :: istatus
    !!----------------------------------------------------------------------

  istatus = NF90_CREATE(cn_fout, NF90_CLOBBER, kcid)
  IF (l_moor) THEN
     npiout=1
     npjout=1
     istatus = nf90_DEF_DIM(kcid, cn_dim_t, NF90_UNLIMITED, id_dim(3))
     istatus = NF90_DEF_VAR(kcid, cn_var_time, NF90_FLOAT, id_dim(3), id_vart)
        ; istatus = NF90_PUT_ATT(kcid, id_vart, 'long_name','Time')
        ; istatus = NF90_PUT_ATT(kcid, id_vart, 'units', ca_units)
     istatus = nf90_DEF_DIM(kcid, cn_dim_x,npiout, id_dim(1))
     istatus = NF90_DEF_VAR(kcid, cn_var_lon, NF90_FLOAT, id_dim(1), id_varx)
        ; istatus = NF90_PUT_ATT(kcid, id_varx, 'long_name','Longitude')
        ; istatus = NF90_PUT_ATT(kcid, id_varx, 'units', 'degrees')
     istatus = nf90_DEF_DIM(kcid, cn_dim_y, npjout, id_dim(2))
     istatus = NF90_DEF_VAR(kcid, cn_var_lat, NF90_FLOAT, id_dim(2), id_vary)
        ; istatus = NF90_PUT_ATT(kcid, id_vary, 'long_name','latitude')
        ; istatus = NF90_PUT_ATT(kcid, id_vary, 'units', 'degrees')
  ELSE
     npiout=npi
     npjout=npj
     istatus = nf90_DEF_DIM(kcid, cn_dim_t, NF90_UNLIMITED, id_dim(3))
     istatus = NF90_DEF_VAR(kcid,'time_counter', NF90_FLOAT, id_dim(3), id_vart)
     istatus = NF90_PUT_ATT(kcid, id_vart, 'long_name','Time')
     istatus = NF90_PUT_ATT(kcid, id_vart, 'units', ca_units)
     istatus = nf90_DEF_DIM(kcid,'x',npiout, id_dim(1))
     istatus = nf90_DEF_DIM(kcid,'y',npjout, id_dim(2))
     istatus = NF90_DEF_VAR(kcid,'nav_lon', NF90_FLOAT, id_dim(1:2), id_varx)
     istatus = NF90_PUT_ATT(kcid, id_varx, 'long_name','Longitude')
     istatus = NF90_PUT_ATT(kcid, id_varx, 'units', 'degrees')
     istatus = NF90_DEF_VAR(kcid,'nav_lat', NF90_FLOAT, id_dim(1:2), id_vary)
     istatus = NF90_PUT_ATT(kcid, id_vary, 'long_name','latitude')
     istatus = NF90_PUT_ATT(kcid, id_vary, 'units', 'degrees')

  ENDIF
  istatus = nf90_DEF_VAR(kcid,TRIM(cn_var_nm), NF90_FLOAT, id_dim(1:3), id_var) 
  istatus = NF90_PUT_ATT(kcid, id_var, 'long_name', TRIM(cn_var_ln))
  istatus = NF90_PUT_ATT(kcid, id_var, 'units', TRIM(cn_var_unit))
  istatus = NF90_PUT_ATT(kcid, id_var, 'missing_value', 0.)
  istatus = NF90_ENDDEF(kcid)

  ALLOCATE(dres(npiout,npjout,1))

  istatus = NF90_PUT_VAR(kcid, id_vart,dg_time(1:npt))
  istatus = NF90_PUT_VAR(kcid, id_varx,dlon)
  istatus = NF90_PUT_VAR(kcid, id_vary,dlat)
  END SUBROUTINE CreateOutput

END PROGRAM tid_predict
