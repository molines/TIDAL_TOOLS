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
  USE netcdf

  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER :: nconst=11
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

  REAL(KIND=8), DIMENSION(nconst)             :: domega, dft, dut, dvt, dvt0
  REAL(KIND=8)                                :: dl_amp, dl_ph, dl_hfrac
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE     :: dg_time
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE   :: dtemp2d, dlon, dlat
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dres
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dtempx, dtempy

  CHARACTER(LEN=100), DIMENSION(jpmax_harmo) :: cname
  CHARACTER(LEN=120)                         :: cldum
  CHARACTER(LEN=120)                         :: cf_in
  CHARACTER(LEN=120)                         :: cf_namli
  CHARACTER(LEN=100)                         :: cl_units     
  CHARACTER(LEN=120)                         :: cn_fharm, cn_fout  
  CHARACTER(LEN=100)                         :: cn_var_nm
  CHARACTER(LEN=255)                         :: cn_var_ln, cn_var_unit

  LOGICAL :: l_moor =.TRUE.
  LOGICAL :: lnc4   =.FALSE.
  LOGICAL :: l_exist

  NAMELIST / prediction_param / cn_var_nm, cn_var_ln, cn_var_unit,  cn_fharm, cn_fout
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017
  !! $Id$
  !! Copyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  tid_predict -n NAMLIST-file -h HARM-file'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Perform tide prediction using harmonic constituents.' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -n NAMLIST-file : give the name of the namelist to be used.' 
     PRINT *,'       -h HARM-file : give the name of the constituents file.'
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
     CASE ( '-h'   ) ; CALL getarg(ijarg, cf_in    ) ; ijarg=ijarg+1
        ! option
     CASE ( '-nc4' ) ; lnc4 = .TRUE.
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  ! sanity check (tbd) (namelist cf_in)

  ! 0. Read the namelist
  ! -------------------------------

  OPEN( inum, file = 'namelist',status = 'old' ,  form   = 'formatted' )
  READ( inum, NML = prediction_param )
  CLOSE( inum )

  PRINT *, ' '
  PRINT *, 'Variable to be predicted = ', TRIM(cn_var_nm)
  PRINT *, 'Variable long name       = ', TRIM(cn_var_ln)
  PRINT *, 'Variable units           = ', TRIM(cn_var_unit)
  PRINT *, 'Harmonic data file       = ', TRIM(cn_fharm)
  PRINT *, 'Prediction output file   = ', TRIM(cn_fout)
  PRINT *, ' '

  cname(1 )='M2'
  cname(2 )='N2'
  cname(3 )='S2'
  cname(4 )='K2'
  cname(5 )='Q1'
  cname(6 )='K1'
  cname(7 )='O1'
  cname(8 )='P1'
  cname(9 )='Mf'
  cname(10)='Mm'
  cname(11)='M4'
  !         cname(12)='M6'

  ! Read input file
  !----------------- 
  istatus = NF90_OPEN (cf_in, NF90_NOWRITE, ncid)
  istatus = NF90_INQ_DIMID (ncid,'x',id_var)
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

  istatus=NF90_INQ_DIMID (ncid,'time',id_var)
  IF (istatus /= NF90_NOERR) THEN
     istatus=NF90_INQ_DIMID (ncid,'time_counter',id_var)
     IF (istatus /= NF90_NOERR) THEN
        PRINT *, 'Pb with time-dimension in file: ', cf_in
        STOP
     ENDIF
  ENDIF
  istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=npt)

  istatus = NF90_INQ_VARID(ncid,'time',id_var)
  IF (istatus /= NF90_NOERR) THEN
     istatus = NF90_INQ_VARID(ncid,'time_counter',id_var)
     IF (istatus /= NF90_NOERR) THEN
        PRINT *, 'Pb with time-dimension in file: ', cf_in
        STOP
     ENDIF
  ENDIF
  ALLOCATE(dg_time(npt))
  istatus = NF90_GET_VAR(ncid,id_var, dg_time)
  istatus = NF90_GET_ATT(ncid, id_var, "units", cl_units)
  READ(cl_units,7000) iyy, imm, idd, ihh, imn, isec
7000 FORMAT('seconds since ', I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)

  ijuli = julday(imm,idd,iyy)             ! Julian day of time origin in data file
  ijul0 = ijuli + FLOOR(dg_time(1)/86400.) ! Julian day of first time dump in file
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
     istatus = NF90_INQ_VARID(ncid,'longitude',id_var)
     istatus = NF90_GET_VAR(ncid,id_var, dlon)
     istatus = NF90_INQ_VARID(ncid,'latitude',id_var)
     istatus = NF90_GET_VAR(ncid,id_var, dlat)
  ELSE
     istatus=NF90_INQ_DIMID (ncid,'x',id_var)
     istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=npi)
     istatus=NF90_INQ_DIMID (ncid,'y',id_var)
     istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=npj)
     ALLOCATE(dlon(npi,npj), dlat(npi,npj))
     istatus = NF90_INQ_VARID(ncid,'nav_lon',id_var)
     istatus = NF90_GET_VAR(ncid,id_var, dlon)
     istatus = NF90_INQ_VARID(ncid,'nav_lat',id_var)
     istatus = NF90_GET_VAR(ncid,id_var, dlat)
  ENDIF

  WRITE(6,*)'INPUT FILE START DATE (ddmmyy) ', idd, imm, iyy
  WRITE(6,*) 'Number OF time dumps ', npt
  WRITE(6,*) 'Number OF harmonics  ', nconst
  IF (l_moor) THEN
     WRITE(6,*)'Mooring file indices ', ii0, ij0
  ENDIF
  CALL FLUSH(6)

  istatus = NF90_CLOSE(ncid)

  ! Set tidal waves pulsations and nodal corrections (at initial time)
  !-------------------------------------------------------------------
  CALL tide_pulse( domega(1:nconst), cname(1:nconst) ,nconst)
  CALL tide_vuf( dvt(1:nconst), dut(1:nconst) , dft(1:nconst), cname(1:nconst) ,nconst, iyy, imm, idd, dl_hfrac)
  dvt0 =dvt

  ! Read harmonic data file
  !-------------------------      
  istatus = NF90_OPEN (cn_fharm, NF90_NOWRITE, ncid)
  IF (istatus /= NF90_NOERR) THEN
     PRINT *, 'Can not open file: ', cn_fharm
     STOP
  ENDIF

  istatus=NF90_INQ_DIMID (ncid,'x',id_var)
  IF (istatus /= NF90_NOERR) THEN
     istatus=NF90_INQ_DIMID (ncid,'longitude',id_var)
     IF (istatus /= NF90_NOERR) THEN
        PRINT *, 'Pb with x-dimension in file: ', cn_fharm
        STOP
     ENDIF
  ENDIF
  istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=npih)

  istatus = NF90_INQ_DIMID (ncid,'y',id_var)
  IF (istatus /= NF90_NOERR) THEN
     istatus=NF90_INQ_DIMID (ncid,'latitude',id_var)
     IF (istatus /= NF90_NOERR) THEN
        PRINT *, 'Pb with y-dimension in file: ', cn_fharm
        STOP
     ENDIF
  ENDIF
  istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=npjh)

  ALLOCATE(dtemp2d(npih,npjh))
  IF (l_moor) THEN
     ALLOCATE(dtempx(1,1,nconst), dtempy(1,1,nconst))
  ELSE
     npi=npih
     npj=npjh
     ALLOCATE(dtempx(npi,npj,nconst), dtempy(npi,npj,nconst))
  ENDIF

  DO jt=1, nconst
     istatus = NF90_INQ_VARID( ncid,TRIM(cname(jt))//'_x_elev',id_var)
     IF (istatus /= NF90_NOERR) THEN
        PRINT *, 'Can not find harmonic data ', TRIM(cname(jt))//'_x_elev in file: ', cn_fharm
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

     istatus = NF90_INQ_VARID( ncid,TRIM(cname(jt))//'_y_elev',id_var)
     IF (istatus /= NF90_NOERR) THEN
        PRINT *, 'Can not find harmonic data ', TRIM(cname(jt))//'_y_elev in file: ', cn_fharm
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

  ! Build residual
  !---------------

  it0 = FLOOR(dg_time(1)/86400.)*86400.
  istatus = NF90_CREATE(cn_fout, NF90_CLOBBER, ncid)
  IF (l_moor) THEN
     npiout=1
     npjout=1
     istatus = nf90_DEF_DIM(ncid,'time', NF90_UNLIMITED, id_dim(3))
     istatus = NF90_DEF_VAR(ncid,'time', NF90_FLOAT, id_dim(3), id_vart)
     istatus = NF90_PUT_ATT(ncid, id_vart, 'long_name','Time')
     istatus = NF90_PUT_ATT(ncid, id_vart, 'units', cl_units)
     istatus = nf90_DEF_DIM(ncid,'longitude',npiout, id_dim(1))
     istatus = NF90_DEF_VAR(ncid,'longitude', NF90_FLOAT, id_dim(1), id_varx)
     istatus = NF90_PUT_ATT(ncid, id_varx, 'long_name','Longitude')
     istatus = NF90_PUT_ATT(ncid, id_varx, 'units', 'degrees')
     istatus = nf90_DEF_DIM(ncid,'latitude', npjout, id_dim(2))
     istatus = NF90_DEF_VAR(ncid,'latitude', NF90_FLOAT, id_dim(2), id_vary)
     istatus = NF90_PUT_ATT(ncid, id_vary, 'long_name','latitude')
     istatus = NF90_PUT_ATT(ncid, id_vary, 'units', 'degrees')
  ELSE
     npiout=npi
     npjout=npj
     istatus = nf90_DEF_DIM(ncid,'time_counter', NF90_UNLIMITED, id_dim(3))
     istatus = NF90_DEF_VAR(ncid,'time_counter', NF90_FLOAT, id_dim(3), id_vart)
     istatus = NF90_PUT_ATT(ncid, id_vart, 'long_name','Time')
     istatus = NF90_PUT_ATT(ncid, id_vart, 'units', cl_units)
     istatus = nf90_DEF_DIM(ncid,'x',npiout, id_dim(1))
     istatus = nf90_DEF_DIM(ncid,'y',npjout, id_dim(2))
     istatus = NF90_DEF_VAR(ncid,'nav_lon', NF90_FLOAT, id_dim(1:2), id_varx)
     istatus = NF90_PUT_ATT(ncid, id_varx, 'long_name','Longitude')
     istatus = NF90_PUT_ATT(ncid, id_varx, 'units', 'degrees')
     istatus = NF90_DEF_VAR(ncid,'nav_lat', NF90_FLOAT, id_dim(1:2), id_vary)
     istatus = NF90_PUT_ATT(ncid, id_vary, 'long_name','latitude')
     istatus = NF90_PUT_ATT(ncid, id_vary, 'units', 'degrees')

  ENDIF
  istatus = nf90_DEF_VAR(ncid,TRIM(cn_var_nm), NF90_FLOAT, id_dim(1:3), id_var) 
  istatus = NF90_PUT_ATT(ncid, id_var, 'long_name', TRIM(cn_var_ln))
  istatus = NF90_PUT_ATT(ncid, id_var, 'units', TRIM(cn_var_unit))
  istatus = NF90_PUT_ATT(ncid, id_var, 'missing_value', 0.)
  istatus = NF90_ENDDEF(ncid)

  ALLOCATE(dres(npiout,npjout,1))

  istatus = NF90_PUT_VAR(ncid, id_vart,dg_time(1:npt))
  istatus = NF90_PUT_VAR(ncid, id_varx,dlon)
  istatus = NF90_PUT_VAR(ncid, id_vary,dlat)

  DO jt=1,npt 

     ijul0 = ijuli + FLOOR(dg_time(jt)/86400.) ! Julian day of first time dump in file
     CALL caldat(ijul0,imm,idd,iyy)
     dl_hfrac = ijuli+dg_time(jt)/86400.-ijul0
     !
     !          write(6,*)
     !         write(6,*)' time ',dg_time(jt)
     !         write(6,*)' ijuli ijul0 dl_hfrac ',ijuli,ijul0,dl_hfrac
     !         write(6,*) imm,idd,iyy
     !         write(6,*)
     !
     CALL tide_vuf( dvt(1:nconst), dut(1:nconst) , dft(1:nconst), cname(1:nconst) ,nconst, iyy, imm, idd, dl_hfrac)
     dres(:,:,:) = 0.d0
     DO jn=1, nconst
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

  FUNCTION julday(kmm,kid,kiyyy)
    ! ------------------------------------------------------------------
    !                 FUNCTION JULDAY
    !                 ***************
    !   PURPOSE:
    !   --------
    !    This routine returns the julian day number which begins at noon
    !  of the calendar date specified by month kmm, day kid, and year kiyyy.
    !  positive year signifies a.d.; negative, b.c.  (remember that the
    !  year after 1 b.c. was 1 a.d.)
    !  routine handles changeover to gregorian calendar on oct. 15, 1582.
    !
    !    METHOD:
    !    -------
    !     This routine comes directly from the Numerical Recipe Book,
    !   press et al., numerical recipes, cambridge univ. press, 1986.
    !
    !    ARGUMENTS:
    !    ----------
    !     kmm     : input, corresponding month
    !     kid     : input, corresponding day
    !     kiyyy   : input, corresponding year, positive if a.d, negative b.c.
    !      
    !     
    !   AUTHOR:
    !   ------
    !     1998: J.M. Molines for the Doctor form.
    !-----------------------------------------------------------------
    !
    ! ... Declarations
    !
    !     INTEGER jpgreg
    INTEGER, INTENT(IN) :: kmm, kid, kiyyy
    INTEGER, PARAMETER :: jpgreg = 15+31*(10+12*1582)
    INTEGER :: il_kiyyy, julday
    INTEGER iy, im, ia

    il_kiyyy = kiyyy

    ! ... Year 0 never existed ...
    IF (il_kiyyy.EQ.0) STOP 101
    !
    IF (il_kiyyy.LT.0) il_kiyyy = il_kiyyy + 1
    IF (kmm.GT.2) THEN
       iy = il_kiyyy
       im = kmm + 1
    ELSE
       iy = il_kiyyy - 1
       im = kmm + 13
    END IF
    !
    julday = INT(365.25*iy) + INT(30.6001*im) + kid + 1720995 
    IF (kid+31*(kmm+12*il_kiyyy).GE.jpgreg) THEN
       ia = INT(0.01*iy)
       julday = julday + 2 - ia + INT(0.25*ia) 
    END IF
    RETURN
  END FUNCTION julday

  SUBROUTINE caldat(kjulian,kmm,kid,kiyyy)
    ! -------------------------------------------------------------------
    !                   SUBROUTINE CALDAT
    !                   *****************
    !   PURPOSE:
    !   --------
    !    This routine convert a julian day in calendar date.
    !    It is the inverse of the function julday.  
    !
    !    METHOD:
    !    -------
    !     This routine comes directly from the Numerical Recipe Book,
    !   press et al., numerical recipes, cambridge univ. press, 1986.
    !
    !    ARGUMENTS:
    !    ----------
    !     kjulian : input julian day number
    !     kmm     : output, corresponding month
    !     kid     : output, corresponding day
    !     kiyyy   : output, corresponding year, positive if a.d, negative b.c.
    !      
    !    
    !   AUTHOR:
    !   ------
    !     1998: J.M. Molines for the Doctor form.
    !------------------------------------------------------------------------
    ! ... Declarations:
    !
    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: kjulian
    INTEGER, INTENT(OUT) :: kmm, kid, kiyyy

    INTEGER, PARAMETER   :: jpgreg = 2299161
    INTEGER :: ia, ialpha, ib, ic, id, ie
    !
    ! ... Cross over to Gregorian Calendar produces this correction:
    !
    IF ( kjulian .GE. jpgreg) THEN
       ialpha = INT ((( kjulian - 1867216) - 0.25)/36524.25 )
       ia     = kjulian +1 + ialpha -INT (0.25*ialpha)
    ELSE
       ia = kjulian
    END IF
    !
    ib = ia + 1524
    ic = INT (6680. + (( ib -2439870) - 122.1)/365.25 )
    id = 365* ic + INT (0.25*ic)
    ie = INT (( ib - id )/30.6001)
    !
    kid = ib - id - INT (30.6001*ie)
    kmm = ie -1
    IF ( kmm .GT. 12 ) kmm = kmm - 12
    kiyyy = ic - 4715
    IF ( kmm   .GT. 2 ) kiyyy = kiyyy - 1
    IF ( kiyyy .LE. 0 ) kiyyy = kiyyy - 1
    RETURN
  END SUBROUTINE caldat
END PROGRAM tid_predict
