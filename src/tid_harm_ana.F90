PROGRAM tid_harm_ana
  !!======================================================================
  !!                     ***  PROGRAM tid_harm_ana  ***
  !!=====================================================================
  !!  ** Purpose : Perform harmonic analysis of 2D fields, in particular
  !!         NEMO model output (default)
  !!
  !!  ** Method  : 
  !!
  !! History :  1.0  : 12/2018  : J.M. Molines : From Mercator from ???
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE tide
  USE surdetermine
  USE netcdf

  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER :: jpconst=20
  INTEGER(KIND=4)            :: jt, jn, ji, jj, jc, jfil
  INTEGER(KIND=4), DIMENSION(jpconst) :: id_varhx, id_varhy
  INTEGER(KIND=4), DIMENSION(3) :: id_dim
  INTEGER(KIND=4) :: id_var, id_vart, id_varx, id_vary, id_var_lon, id_var_lat
  INTEGER(KIND=4) :: ijul0, ijuli, ijule
  INTEGER(KIND=4) :: narg, iargc, init, ixtype
  INTEGER(KIND=4) :: npi, npj, npt, npiout, npjout
  INTEGER(KIND=4) :: ncid, istatus
  INTEGER(KIND=4) :: iyy, imm, idd, ihh, imn, isec
  INTEGER(KIND=4) :: iyye, imme, idde
  INTEGER(KIND=4) :: nhc, nhan, nsp, nun

  INTEGER(KIND=2), DIMENSION(:,:,:), ALLOCATABLE :: int2tabglo

  REAL(8) :: dn_tsamp ! Sampling period in seconds
  REAL(8) :: dscale_factor =1.d0, dadd_offset =0.d0, dmissval=1.d20
  REAL(8) :: dl_temp, dl_time, dl_amp, dl_ph, dl_hfrac, dn_ave
  REAL(8), DIMENSION(jpconst)              :: domega, dft, dut, dvt, dvt0
  REAL(8), DIMENSION(:),       ALLOCATABLE :: dg_time
  REAL(8), DIMENSION(:,:),     ALLOCATABLE :: dlon_r4, dlat_r4
  REAL(8), DIMENSION(:,:,:),   ALLOCATABLE :: dtmp_r4
  REAL(8), DIMENSION(:,:,:),   ALLOCATABLE :: dana_temp
  REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: dana_amp

  CHARACTER(len=4  ), DIMENSION(jpmax_harmo) :: cname
  CHARACTER(len=120)                         :: cn_fharm, cf_in
  CHARACTER(len=33 )                         :: ca_units, cl_dum
  CHARACTER(len=20 )                         :: cn_v_in, cn_v_out_x, cn_v_out_y

  LOGICAL :: ln_moor, l_exist, ln_short

  NAMELIST /analysis_param/ ln_moor, ln_short, &
       cn_v_in, cn_v_out_x, cn_v_out_y, cn_fharm ,&
       dn_tsamp, dn_ave 

  !!----------------------------------------------------------------------
  !! TIDAL_TOOLS , MEOM 2018
  !! $Id$
  !! Copyright (c) 2018, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/TIDAL_TOOLS_CeCILL.txt)
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  

  ! 0. Read the namelist
  ! -------------------------------
  INQUIRE( FILE = 'namelist', EXIST = l_exist )

  IF ( .NOT. l_exist ) STOP 'Namelist not found'

  OPEN( 20, file = 'namelist', &
       status = 'old'       , &
       form   = 'formatted' )
  READ( 20, NML = analysis_param )
  CLOSE( 20 )

  PRINT *, ' '
  PRINT *, 'Variable to be analysed = ', TRIM(cn_v_in)
  PRINT *, 'Sampling period = ', dn_tsamp
  PRINT *, 'Output file name = ', TRIM(cn_fharm)
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
  cname(12)='M6'
  cname(13)='S1'
  cname(14)='Msqm'
  cname(15)='Mtm'
  cname(16)='_2N2'
  cname(17)='MU2'
  cname(18)='NU2'
  cname(19)='L2'
  cname(20)='T2'

  narg = iargc()
  IF (narg<1) THEN
     PRINT *, 'Need at least one input file'
     PRINT *, 'tide_ana usage:' 
     PRINT *, './tide_ana file_list'
     PRINT *, 'EXEMPLE: ./tide_ana BISCAY-T05-y2008*_HF_gridT.nc'
     PRINT *, 'Output harmonic file in res_harm.nc'
     STOP
  ENDIF

  CALL getarg (1, cf_in)

  CALL read_file

  ! Set tidal waves pulsations and nodal corrections (at initial time)
  !-------------------------------------------------------------------
  CALL tide_pulse( domega(1:jpconst), cname(1:jpconst) ,jpconst)
  CALL tide_vuf( dvt0(1:jpconst), dut(1:jpconst) , dft(1:jpconst), cname(1:jpconst) ,jpconst, iyy, imm, idd, dl_hfrac)

  ! Perform harmonic analysis
  !--------------------------

  IF (ln_moor) THEN
     ALLOCATE(dtmp_r4(1,1,1))
     ALLOCATE(dana_temp(1,1,2*jpconst))
     ALLOCATE(dana_amp(1,1,jpconst,2))
     npiout=1
     npjout=1
  ELSE
     ALLOCATE(dtmp_r4(npi,npj,1))
     ALLOCATE(dana_temp(npi,npj,2*jpconst))
     ALLOCATE(dana_amp(npi,npj,jpconst,2))
     npiout=npi
     npjout=npj
  ENDIF

  dana_temp(:,:,:) = 0.d0
  dl_time = dn_ave * dn_tsamp
  nhan  = 0  ! Total number of time dumps
  nsp   = 0  ! Initialize length of sparse matrix

  !MBK
  WRITE(6,*)
  WRITE(6,*)' ztime : ',dl_time
  WRITE(6,*)' r_ave : ',dn_ave
  WRITE(6,*)' tsamp : ',dn_tsamp
  WRITE(6,*)
  CALL FLUSH(6)
  !MBK
  DO jfil=1,narg
     CALL getarg (jfil, cf_in)
     ! START LOOP ON INPUT FILES HERE
     ! open file again
     istatus = NF90_OPEN (cf_in, NF90_NOWRITE, ncid)
     IF (istatus /= NF90_NOERR) THEN
        PRINT *, 'Can not open file: ', cf_in
        STOP
     ELSE
        WRITE(6,*)''
        WRITE(6,'(a,2x,a)') 'OPEN FILE: ',TRIM(cf_in)
        CALL FLUSH(6)
     ENDIF

     ! get number of time dumps
     istatus=NF90_INQ_DIMID (ncid,'time',id_var)
     IF (istatus /= NF90_NOERR) THEN
        istatus=NF90_INQ_DIMID (ncid,'time_counter',id_var)
        IF (istatus /= NF90_NOERR) THEN
           PRINT *, 'Pb with time-dimension in file: ', cf_in
           STOP
        ENDIF
     ENDIF
     istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=npt)
     PRINT *, 'NUMBER OF TIME DUMPS TO READ IN FILE: ', npt

     istatus = NF90_INQ_VARID(ncid,TRIM(cn_v_in),id_var)
     IF (istatus /= NF90_NOERR) THEN
        PRINT *, 'CAN NOT FIND VARIABLE: ', TRIM(cn_v_in), ' IN FILE ', cf_in
        STOP
     ELSE
        istatus = NF90_GET_ATT(ncid, id_var, "missing_value", dmissval)
        istatus = NF90_INQUIRE_VARIABLE(ncid, id_var, cl_dum, ixtype)
        IF (ixtype==NF90_SHORT) THEN 
           ln_short=.TRUE.
           istatus = NF90_GET_ATT(ncid, id_var, "add_offset", dadd_offset)
           istatus = NF90_GET_ATT(ncid, id_var, "scale_factor", dscale_factor)
           IF (jfil==1) THEN
              IF (ln_moor)      ALLOCATE(int2tabglo(1,1,1))
              IF (.NOT.ln_moor) ALLOCATE(int2tabglo(npi,npj,1))
           ENDIF
        ENDIF
     ENDIF

     DO jt=1, npt
        !  Get variable
        IF (ln_short) THEN
           istatus = NF90_GET_VAR(ncid,id_var, int2tabglo, &
                start = (/ 1, 1, jt /),count = (/ npiout, npjout, 1 /))
           WHERE(int2tabglo .NE. dmissval) dtmp_r4 = int2tabglo*dscale_factor+dadd_offset
           WHERE(int2tabglo .EQ. dmissval) dtmp_r4 = 0.
        ELSE
           istatus = NF90_GET_VAR(ncid,id_var, dtmp_r4, &
                start = (/ 1, 1, jt /),count = (/ npiout, npjout, 1 /))
           WHERE(dtmp_r4 .EQ. dmissval) dtmp_r4 = 0.
        ENDIF

        dl_time = dl_time + dn_tsamp

        ! Update nodal corrections:
        ijul0 = ijuli + FLOOR(dl_time/86400.) ! Julian day of first time dump in file
        CALL caldat(ijul0,imm,idd,iyy)
        dl_hfrac = ijuli+dl_time/86400.-ijul0

        WRITE(6,*)' mm dd yy :',imm,idd,iyy
        WRITE(6,'(a,1x,f12.2,1x,a,1x,f12.2,1x,a,1x,f18.8)')' ztime    : ',dl_time,'  ijul0 : ',ijul0,' HFRAC : ', dl_hfrac
        CALL FLUSH(6)

        CALL tide_vuf( dvt(1:jpconst), dut(1:jpconst) , dft(1:jpconst), cname(1:jpconst) ,jpconst, iyy, imm, idd, dl_hfrac)

        nhan = nhan + 1
        nhc = 0
        DO jn = 1, jpconst
           DO jc = 1,2
              nhc = nhc + 1
              nsp = nsp + 1
              dl_temp =(     MOD(jc,2) * dft(jn) *COS(domega(jn)*dl_time + dvt0(jn) + dut(jn))  &
                   +(1.-MOD(jc,2))* dft(jn) *SIN(domega(jn)*dl_time + dvt0(jn) + dut(jn)))
              dana_temp(:,:,nhc) = dana_temp(:,:,nhc) + dl_temp*dtmp_r4(:,:,1)                   
              ISPARSE(nsp) = nhan
              JSPARSE(nsp) = nhc
              SPARSEVALUE(nsp) = dl_temp
           END DO
        END DO
     END DO
     istatus = NF90_CLOSE(ncid)
     ! STOP LOOP ON INPUT FILES HERE
  END DO


  ijule = ijul0 + FLOOR(dl_time/86400.) - 1
  CALL caldat(ijule,imme,idde,iyye)  

  PRINT *,''
  PRINT *,'START DATE (ddmmyy)        : ', idd, imm, iyy
  PRINT *,'END DATE (ddmmyy)          : ', idde, imme, iyye
  PRINT *,'NUMBER OF DAYS FOR ANALYSIS: ', ijule-ijul0 + 1
  PRINT *,'SAMPLING PERIOD (seconds)  : ', FLOOR(dn_tsamp)
  PRINT *,''

  ! Matrix inversion
  NBINCO   = 2*jpconst
  NBSPARSE = nsp
  CALL SUR_DETERMINE(1)

  ! Find solution for each point
  DO jj = 1, npjout
     DO ji = 1, npiout
        nun=0
        DO jn = 1, jpconst
           DO jc = 1,2
              nun = nun + 1
              TAB4(nun)=dana_temp(ji,jj,nun)
           ENDDO
        ENDDO

        CALL SUR_DETERMINE(0)

        ! Fill output array
        DO jn = 1, jpconst
           dana_amp(ji,jj,jn,1)= TAB7((jn-1)*2+1)
           dana_amp(ji,jj,jn,2)=-TAB7((jn-1)*2+2)
        END DO
     END DO
  END DO

  DEALLOCATE(dana_temp)

  ! WRITE output
  !-------------

  IF (ln_moor) THEN
     DO jn = 1, jpconst
        dl_amp = SQRT(dana_amp(1,1,jn,1)**2+dana_amp(1,1,jn,2)**2) ! Amplitude
        dl_ph  = ATAN2(-dana_amp(1,1,jn,2)/dl_amp,dana_amp(1,1,jn,1)/dl_amp)/rad ! Phase
        PRINT *, 'constituent name / Period (hour) / amplitude / phase (deg)'
        PRINT *, TRIM(cname(jn)), 2.*rpi/domega(jn)/3600., dl_amp, dl_ph
     END DO
  ENDIF

  CALL write_file

CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE read_file

    ! Read input file
    !----------------- 
    istatus = NF90_OPEN (cf_in, NF90_NOWRITE, ncid)
    IF (istatus /= NF90_NOERR) THEN
       PRINT *, 'Can not open file: ', cf_in
       STOP
    ENDIF

    istatus=NF90_INQ_DIMID (ncid,'x',id_var)
    IF (istatus == NF90_NOERR) THEN 
       ln_moor=.FALSE.
    ELSE
       ln_moor=.TRUE.
    ENDIF

    ! Get mooring indice on global grid

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
    istatus = NF90_GET_ATT(ncid, id_var, "units", ca_units)
    READ(ca_units,7000) iyy, imm, idd, ihh, imn, isec
7000 FORMAT('seconds since ', I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)

    ijuli = julday(imm,idd,iyy)             ! Julian day of time origin in data file
    ijul0 = ijuli + FLOOR(dg_time(1)/86400.d0) ! Julian day of first time dump in file
    !rbb
    PRINT*,dg_time(1)/86400.d0

    CALL caldat(ijul0,imm,idd,iyy)          ! Set initial date month/day/year 
    dl_hfrac=ihh+imn/60.d0+isec/3600.d0
    ijuli=ijul0

    DEALLOCATE(dg_time)

    IF (ln_moor) THEN
       ALLOCATE(dlon_r4(1,1), dlat_r4(1,1))
       istatus = NF90_INQ_VARID(ncid,'longitude',id_var)
       istatus = NF90_GET_VAR(ncid,id_var, dlon_r4)
       istatus = NF90_INQ_VARID(ncid,'latitude',id_var)
       istatus = NF90_GET_VAR(ncid,id_var, dlat_r4)
    ELSE
       istatus=NF90_INQ_DIMID (ncid,'x',id_var)
       istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=npi)
       istatus=NF90_INQ_DIMID (ncid,'y',id_var)
       istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=npj)
       ALLOCATE(dlon_r4(npi,npj), dlat_r4(npi,npj))
       istatus = NF90_INQ_VARID(ncid,'nav_lon',id_var)
       istatus = NF90_GET_VAR(ncid,id_var, dlon_r4)
       istatus = NF90_INQ_VARID(ncid,'nav_lat',id_var)
       istatus = NF90_GET_VAR(ncid,id_var, dlat_r4)
    ENDIF

    istatus = NF90_CLOSE(ncid)

  END SUBROUTINE read_file

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE write_file

    !rbb istatus = NF90_CREATE(cn_fharm, NF90_CLOBBER, ncid)
    istatus = NF90_CREATE(cn_fharm, NF90_64BIT_OFFSET, ncid)

    IF (ln_moor) THEN
       npiout=1
       npjout=1
       istatus = NF90_DEF_DIM(ncid,'longitude',npiout, id_dim(1))
       istatus = NF90_DEF_VAR(ncid,'longitude', NF90_FLOAT, id_dim(1), id_var_lon)
       istatus = NF90_PUT_ATT(ncid, id_var_lon, 'long_name','Longitude')
       istatus = NF90_PUT_ATT(ncid, id_var_lon, 'units', 'degrees')
       istatus = NF90_DEF_DIM(ncid,'latitude', npjout, id_dim(2))
       istatus = NF90_DEF_VAR(ncid,'latitude', NF90_FLOAT, id_dim(2), id_var_lat)
       istatus = NF90_PUT_ATT(ncid, id_var_lat, 'long_name','latitude')
       istatus = NF90_PUT_ATT(ncid, id_var_lat, 'units', 'degrees')
    ELSE
       npiout=npi
       npjout=npj
       istatus = NF90_DEF_DIM(ncid,'x',npiout, id_dim(1))
       istatus = NF90_DEF_DIM(ncid,'y',npjout, id_dim(2))
       istatus = NF90_DEF_VAR(ncid,'nav_lon', NF90_FLOAT, id_dim(1:2), id_var_lon)
       istatus = NF90_PUT_ATT(ncid, id_var_lon, 'long_name','Longitude')
       istatus = NF90_PUT_ATT(ncid, id_var_lon, 'units', 'degrees')
       istatus = NF90_DEF_VAR(ncid,'nav_lat', NF90_FLOAT, id_dim(1:2), id_var_lat)
       istatus = NF90_PUT_ATT(ncid, id_var_lat, 'long_name','latitude')
       istatus = NF90_PUT_ATT(ncid, id_var_lat, 'units', 'degrees')
    ENDIF

    DO jn = 1, jpconst
       istatus = NF90_DEF_VAR(ncid,TRIM(cname(jn))//TRIM(cn_v_out_x), NF90_FLOAT, id_dim(1:2), id_varhx(jn)) 
       istatus = NF90_PUT_ATT(ncid, id_varhx(jn), 'long_name',TRIM(cname(jn))//TRIM(cn_v_out_x) )
       istatus = NF90_PUT_ATT(ncid, id_varhx(jn), 'units', 'meters')
       istatus = NF90_PUT_ATT(ncid, id_varhx(jn), 'missing_value', 0.)

       istatus = NF90_DEF_VAR(ncid,TRIM(cname(jn))//TRIM(cn_v_out_y), NF90_FLOAT, id_dim(1:2), id_varhy(jn)) 
       istatus = NF90_PUT_ATT(ncid, id_varhy(jn), 'long_name',TRIM(cname(jn))//TRIM(cn_v_out_y) )
       istatus = NF90_PUT_ATT(ncid, id_varhy(jn), 'units', 'meters')
       istatus = NF90_PUT_ATT(ncid, id_varhy(jn), 'missing_value', 0.)
    END DO

    istatus = NF90_ENDDEF(ncid)

    istatus = NF90_PUT_VAR(ncid, id_var_lon,dlon_r4)
    istatus = NF90_PUT_VAR(ncid, id_var_lat,dlat_r4)
    DO jn = 1, jpconst
       istatus = NF90_PUT_VAR(ncid, id_varhx(jn),dana_amp(:,:,jn,1))
       istatus = NF90_PUT_VAR(ncid, id_varhy(jn),dana_amp(:,:,jn,2))
    END DO

    istatus = NF90_CLOSE(ncid)

  END SUBROUTINE write_file


  FUNCTION julday(kmm,kid,kiyyy)
    !CC ------------------------------------------------------------------
    !CC                 FUNCTION JULDAY
    !CC                 ***************
    !CC   PURPOSE:
    !CC   --------
    !CC    This routine returns the julian day number which begins at noon
    !CC  of the calendar date specified by month kmm, day kid, and year kiyyy.
    !CC  positive year signifies a.d.; negative, b.c.  (remember that the
    !CC  year after 1 b.c. was 1 a.d.)
    !CC  routine handles changeover to gregorian calendar on oct. 15, 1582.
    !CC
    !C    METHOD:
    !C    -------
    !C     This routine comes directly from the Numerical Recipe Book,
    !C   press et al., numerical recipes, cambridge univ. press, 1986.
    !C
    !C    ARGUMENTS:
    !C    ----------
    !C     kmm     : input, corresponding month
    !C     kid     : input, corresponding day
    !C     kiyyy   : input, corresponding year, positive if a.d, negative b.c.
    !C      
    !C     
    !C   AUTHOR:
    !C   ------
    !C     1998: J.M. Molines for the Doctor form.
    !CC -----------------------------------------------------------------
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
  !###
  SUBROUTINE caldat(kjulian,kmm,kid,kiyyy)
    !CC -------------------------------------------------------------------
    !CC                   SUBROUTINE CALDAT
    !CC                   *****************
    !CC   PURPOSE:
    !CC   --------
    !CC    This routine convert a julian day in calendar date.
    !CC    It is the inverse of the function julday.  
    !CC
    !C    METHOD:
    !C    -------
    !C     This routine comes directly from the Numerical Recipe Book,
    !C   press et al., numerical recipes, cambridge univ. press, 1986.
    !C
    !C    ARGUMENTS:
    !C    ----------
    !C     kjulian : input julian day number
    !C     kmm     : output, corresponding month
    !C     kid     : output, corresponding day
    !C     kiyyy   : output, corresponding year, positive if a.d, negative b.c.
    !C      
    !C    
    !C   AUTHOR:
    !C   ------
    !C     1998: J.M. Molines for the Doctor form.
    !CC------------------------------------------------------------------------
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

END PROGRAM tid_harm_ana
