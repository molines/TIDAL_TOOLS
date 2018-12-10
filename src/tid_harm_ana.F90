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

  INTEGER(KIND=4)  :: npconst=20
  INTEGER(KIND=4)            :: jt, jn, ji, jj, jc, jfil
  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: id_varhx, id_varhy, id_varA, id_varG
  INTEGER(KIND=4), DIMENSION(3) :: id_dim
  INTEGER(KIND=4) :: id_var, id_vart, id_varx, id_vary, id_var_lon, id_var_lat
  INTEGER(KIND=4) :: ijul0, ijuli, ijule,ijulb
  INTEGER(KIND=4) :: narg, iargc, ijarg, init, ixtype, inum=20
  INTEGER(KIND=4) :: npi, npj, npt, npiout, npjout, nfiles
  INTEGER(KIND=4) :: ncid, istatus
  INTEGER(KIND=4) :: iyy, imm, idd, ihh, imn, isec
  INTEGER(KIND=4) :: iyye, imme, idde
  INTEGER(KIND=4) :: iyyb, immb, iddb
  INTEGER(KIND=4) :: nhc, nhan, nsp, nun

  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: int2tabglo
  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: mask

  REAL(8) :: dn_tsamp ! Sampling period in seconds
  REAL(8) :: dscale_factor =1.d0, dadd_offset =0.d0, dmissval=1.d20
  REAL(8) :: dl_temp, dl_time, dl_amp, dl_ph, dl_hfrac, dn_ave
  REAL(8) :: dareatot, dcoef, drad

  REAL(8), DIMENSION(:),       ALLOCATABLE :: domega, dft, dut, dvt, dvt0, domega_deg_p_sec
  REAL(8), DIMENSION(:,:),     ALLOCATABLE :: dbeat_period
  REAL(8), DIMENSION(:),       ALLOCATABLE :: dg_time
  REAL(8), DIMENSION(:,:),     ALLOCATABLE :: darea2d, de1t, de2t
  REAL(8), DIMENSION(:,:),     ALLOCATABLE :: dlon_r4, dlat_r4
  REAL(8), DIMENSION(:,:),     ALLOCATABLE :: dtmp_r4
  REAL(8), DIMENSION(:,:,:),   ALLOCATABLE :: dana_temp
  REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: dana_amp

  CHARACTER(len=4  ), DIMENSION(jpmax_harmo) :: cname
  CHARACTER(len=120 ), DIMENSION(:), ALLOCATABLE :: cf_lst
  CHARACTER(len=120)                         :: cn_fharm, cf_in
  CHARACTER(len=120)                         :: cf_namli
  CHARACTER(len=120)                         :: cldum
  CHARACTER(len=120)                         :: cn_v_in, cn_v_out_x, cn_v_out_y
  CHARACTER(len=120)                         :: cdim_x,     cdim_y,     cdim_t
  CHARACTER(len=120)                         :: cdim_x_out, cdim_y_out, cdim_t_out
  CHARACTER(len=120)                         :: cv_lon,     cv_lat,    cv_time
  CHARACTER(len=120)                         :: cv_lon_out, cv_lat_out,cv_time_out
  CHARACTER(len=120)                         :: cf_hgr='mesh_hgr.nc'
  CHARACTER(len=120)                         :: cf_msk='mask.nc'
  CHARACTER(len=120)                         :: ca_units
  CHARACTER(len=120)                         :: ca_miss_in, ca_miss_out

  CHARACTER(LEN=1)                           :: cn_cgrid='T'  ! either one of T U V or F

  LOGICAL :: ln_moor
  LOGICAL :: l_exist
  LOGICAL :: ln_short
  LOGICAL :: lnc4=.FALSE.
  LOGICAL :: lzmean=.FALSE.
  LOGICAL :: lchk=.FALSE.
  LOGICAL :: ln_lonlat_2d_in=.FALSE. , ln_lonlat_2d_out=.TRUE.

  NAMELIST /input_file_format/ cdim_x, cdim_y, cdim_t, &
     &                         cv_lon, cv_lat, cv_time, &
     &                         ln_lonlat_2d_in
  NAMELIST /output_file_format/ cdim_x_out, cdim_y_out, cdim_t_out, &
     &                          cv_lon_out, cv_lat_out, cv_time_out, &
     &                         ln_lonlat_2d_out

  NAMELIST /analysis_param/ ln_moor, ln_short, &
       cn_v_in, cn_v_out_x, cn_v_out_y, cn_fharm ,&
       dn_tsamp, dn_ave 

  NAMELIST /constituents/ cname

  !!----------------------------------------------------------------------
  !! TIDAL_TOOLS , MEOM 2018
  !! $Id$
  !! Copyright (c) 2018, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/TIDAL_TOOLS_CeCILL.txt)
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  ! set default values ( to be included in namelist later )
  ! input files :
  cf_namli='namelist'

  cdim_x='lon'  ; cv_lon='lon'
  cdim_y='lat'  ; cv_lat='lat'
  cdim_t='time' ; cv_time='time'
  ca_miss_in='missing_value'
  ln_lonlat_2d_in=.FALSE.

  ! output file
  cn_fharm='res_harm.nc'
  cdim_x_out='x'  ; cv_lon_out='nav_lon'
  cdim_y_out='y'  ; cv_lat_out='nav_lat'
  cdim_t_out='time' ; cv_time_out='time_counter'
  ca_miss_out='_FillValue'
  ln_lonlat_2d_out=.TRUE.

  ! constituent must be choosen in the namelist
  cname(:)='none'

  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage :  tid_harm_ana -l LST-files [-o HARM-file] [-n NAMLIST-file] [-nc4]'
     PRINT *,'        ... [-zeromean]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Perform the harmonic analysis on the corresponding time series '
     PRINT *,'       represented by the list of files given as arguments.' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -l LST-files : a blank separated list of input files' 
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -n NAMLIST-file : Input namelist file. Default is ''namelist'' ' 
     PRINT *,'       -o HARM-file : Name of the output file with analysed harmonic '
     PRINT *,'           constituents. Default is  res_harm.nc, or set in the namelist.'
     PRINT *,'       -nc4 : output file is in Netcdf4/Hdf5 with chunking and deflation'
     PRINT *,'       -zeromean : subtract spatial mean of the field before analysis.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'        If -zeromean option, mesh_hgr and mask files are required.'
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : ', TRIM(cn_fharm) 
     PRINT *,'       '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      ' 
     PRINT *,'      '
     STOP
  ENDIF

  drad=ACOS(-1.d0)/180.d0

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-l'   ) ; CALL GetFileList 

        ! option
     CASE ( '-n'        ) ; CALL getarg(ijarg, cf_namli   ) ; ijarg=ijarg+1
     CASE ( '-o'        ) ; CALL getarg(ijarg, cn_fharm   ) ; ijarg=ijarg+1
     CASE ( '-nc4'      ) ; lnc4 = .TRUE.
     CASE ( '-zeromean' ) ; lzmean = .TRUE.
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  ! Sanity check
  DO jfil=1,nfiles
     lchk=lchk .OR. chkfile(cf_lst(jfil))
  ENDDO
  IF (lzmean) THEN
    lchk = lchk .OR. chkfile(cf_hgr)
    lchk = lchk .OR. chkfile(cf_msk)
    IF ( lchk ) STOP
  ENDIF

  ! 0. Read the namelist
  ! -------------------------------
  INQUIRE( FILE = cf_namli, EXIST = l_exist )

  IF ( .NOT. l_exist ) STOP 'Namelist not found'

  OPEN( inum, file = cf_namli, status = 'old' , form   = 'formatted' )
  READ( inum, NML = analysis_param )
  REWIND(inum)

  READ( inum, NML = constituents )
  CLOSE( inum )

  ! look for the number of required constituents for analysis
   npconst=COUNT ( (cname(:) /= 'none') )
   ALLOCATE ( id_varhx(npconst), id_varhy(npconst), id_varA(npconst), id_varG(npconst) )
   ALLOCATE ( domega(npconst), dft(npconst), dut(npconst), dvt(npconst), dvt0(npconst), domega_deg_p_sec(npconst) )
   ALLOCATE ( dbeat_period(npconst,npconst) )

  ! Some Control Prints
  PRINT *, ' '
  PRINT *, 'Variable to be analysed = ', TRIM(cn_v_in)
  PRINT *, 'Sampling period = ', dn_tsamp
  PRINT *, 'Output file name = ', TRIM(cn_fharm)
  PRINT *, ' '

  cf_in=cf_lst(1)

  CALL read_file

  ! Set tidal waves pulsations and nodal corrections (at initial time)
  !-------------------------------------------------------------------
  CALL tide_pulse( domega(1:npconst), cname(1:npconst) ,npconst)
  DO jn=1,npconst
       domega_deg_p_sec(jn)= domega(jn)*180.d0*3600.d0/acos(-1.d0)
       dbeat_period(jn,jn)=-1.d10   ! not used it is infinity
  ENDDO

  DO jn=1,npconst
   DO jc=jn+1,npconst
       dbeat_period(jn,jc)=360.d0/ABS(( domega_deg_p_sec(jn) - domega_deg_p_sec(jc) ))/24.d0
       dbeat_period(jc,jn)=dbeat_period(jn,jc)
!    PRINT *,TRIM(cname(jn))//' -- ',TRIM(cname(jc))//' : ',dbeat_period(jn,jc),' Days'
   ENDDO
  ENDDO
  PRINT *,' Maximum Beat period for the selected constituents :',MAXVAL(dbeat_period(:,:) )

   
  CALL tide_vuf( dvt0(1:npconst), dut(1:npconst) , dft(1:npconst), cname(1:npconst) ,npconst, iyy, imm, idd, dl_hfrac)

  ! Perform harmonic analysis
  !--------------------------

  IF (ln_moor) THEN
     ALLOCATE(dtmp_r4(1,1))
     ALLOCATE(dana_temp(1,1,2*npconst))
     ALLOCATE(dana_amp(1,1,npconst,4))
     npiout=1
     npjout=1
  ELSE
     ALLOCATE(dtmp_r4(npi,npj))
     ALLOCATE(dana_temp(npi,npj,2*npconst))
     ALLOCATE(dana_amp(npi,npj,npconst,4))
     npiout=npi
     npjout=npj
  ENDIF

! call to ZeroMeanInit must be done after the size of the domain is known
  IF (lzmean) THEN
    CALL ZeroMeanInit(cn_cgrid)   ! compute area2d and areatot
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

  DO jfil=1,nfiles
     cf_in=cf_lst(jfil)
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
     istatus=NF90_INQ_DIMID (ncid,cdim_t,id_var)
     IF (istatus /= NF90_NOERR) THEN
        PRINT *, 'Pb with time-dimension in file: ', cf_in
        STOP
     ENDIF
     istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=npt)
     PRINT *, 'NUMBER OF TIME DUMPS TO READ IN FILE: ', npt

     istatus = NF90_INQ_VARID(ncid,TRIM(cn_v_in),id_var)
     IF (istatus /= NF90_NOERR) THEN
        PRINT *, 'CAN NOT FIND VARIABLE: ', TRIM(cn_v_in), ' IN FILE ', cf_in
        STOP
     ELSE
        istatus = NF90_GET_ATT(ncid, id_var, ca_miss_in, dmissval)
        istatus = NF90_INQUIRE_VARIABLE(ncid, id_var, cldum, ixtype)
        IF (ixtype==NF90_SHORT) THEN 
           ln_short=.TRUE.
           istatus = NF90_GET_ATT(ncid, id_var, "add_offset", dadd_offset)
           istatus = NF90_GET_ATT(ncid, id_var, "scale_factor", dscale_factor)
           IF (jfil == 1) THEN
              IF (ln_moor) THEN 
                 ALLOCATE(int2tabglo(1,1))
              ELSE 
                 ALLOCATE(int2tabglo(npi,npj))
              ENDIF
           ENDIF
        ENDIF
     ENDIF

     DO jt=1, npt
        !  Get variable
        IF (ln_short) THEN
           istatus = NF90_GET_VAR(ncid,id_var, int2tabglo, &
                start = (/ 1, 1, jt /),count = (/ npiout, npjout, 1 /))
           WHERE(int2tabglo .NE. dmissval) dtmp_r4 = int2tabglo*dscale_factor+dadd_offset
           WHERE(int2tabglo .EQ. dmissval) dtmp_r4 = 0.d0
        ELSE
           istatus = NF90_GET_VAR(ncid,id_var, dtmp_r4, &
                start = (/ 1, 1, jt /),count = (/ npiout, npjout, 1 /))
           WHERE(dtmp_r4 .EQ. dmissval) dtmp_r4 = 0.d0
        ENDIF

        IF (lzmean ) CALL ZeroMean( dtmp_r4 )

        dl_time = dl_time + dn_tsamp  ! seconds

        ! Update nodal corrections:
        ijul0 = ijuli + FLOOR(dl_time/86400.) ! Julian day of first time dump in file
        CALL caldat(ijul0,imm,idd,iyy)
        dl_hfrac = ijuli+dl_time/86400.-ijul0

 !      IF ( MOD ((jt-1), 40 ) == 0 ) THEN
 !      WRITE(6,*)' mm dd yy :',imm,idd,iyy
 !      WRITE(6,'(a,1x,f12.2,1x,a,1x,i12,1x,a,1x,f18.8)')' ztime    : ',dl_time,'  ijul0 : ',ijul0,' HFRAC : ', dl_hfrac
 !      CALL FLUSH(6)
 !      ENDIF

        CALL tide_vuf( dvt(1:npconst), dut(1:npconst) , dft(1:npconst), cname(1:npconst) ,npconst, iyy, imm, idd, dl_hfrac)

        nhan = nhan + 1
        nhc = 0
        DO jn = 1, npconst
           DO jc = 1,2
              nhc = nhc + 1
              nsp = nsp + 1
              dl_temp =(     MOD(jc,2) * dft(jn) *COS(domega(jn)*dl_time + dvt0(jn) + dut(jn))  &
                   +(1.-MOD(jc,2))* dft(jn) *SIN(domega(jn)*dl_time + dvt0(jn) + dut(jn)))
              dana_temp(:,:,nhc) = dana_temp(:,:,nhc) + dl_temp*dtmp_r4(:,:)                   
              ISPARSE(nsp) = nhan
              JSPARSE(nsp) = nhc
              dsparsevalue(nsp) = dl_temp
           END DO
        END DO
     END DO
     istatus = NF90_CLOSE(ncid)
     ! STOP LOOP ON INPUT FILES HERE
  END DO

  ijule = ijulb + FLOOR(dl_time/86400.)
 print *,dl_time/86400.d0, ijul0, ijule
  CALL caldat(ijulb,immb,iddb,iyyb)  
  CALL caldat(ijule,imme,idde,iyye)  

  PRINT *,''
  PRINT *,'START DATE (ddmmyy)        : ', iddb, immb, iyyb
  PRINT *,'END DATE (ddmmyy)          : ', idde, imme, iyye
  PRINT *,'NUMBER OF DAYS FOR ANALYSIS: ', ijule-ijulb + 1
  PRINT *,'SAMPLING PERIOD (seconds)  : ', FLOOR(dn_tsamp)
  PRINT *,''

  ! Matrix inversion
  nbinco   = 2*npconst
  nbsparse = nsp
  CALL sur_determine(kinit=1)

  ! Find solution for each point
  DO jj = 1, npjout
     DO ji = 1, npiout
        nun=0
        DO jn = 1, npconst
           DO jc = 1,2
              nun = nun + 1
              dtab4(nun)=dana_temp(ji,jj,nun)
           ENDDO
        ENDDO

        CALL sur_determine(kinit=0)

        ! Fill output array
        DO jn = 1, npconst
           dana_amp(ji,jj,jn,1)= dtab7((jn-1)*2+1)
           dana_amp(ji,jj,jn,2)=-dtab7((jn-1)*2+2)
           dana_amp(ji,jj,jn,3)= SQRT(dana_amp(ji,jj,jn,1)**2+dana_amp(ji,jj,jn,2)**2)
           dana_amp(ji,jj,jn,4)= ATAN2(-dana_amp(ji,jj,jn,2),dana_amp(ji,jj,jn,1))/drad 
        END DO
     END DO
  END DO

  DEALLOCATE(dana_temp)

  ! WRITE output
  !-------------

  IF (ln_moor) THEN
     DO jn = 1, npconst
        dl_amp = SQRT(dana_amp(1,1,jn,1)**2+dana_amp(1,1,jn,2)**2) ! Amplitude
        dl_ph  = ATAN2(-dana_amp(1,1,jn,2)/dl_amp,dana_amp(1,1,jn,1)/dl_amp)/drad ! Phase
        PRINT *, 'constituent name / Period (hour) / amplitude / phase (deg)'
        PRINT *, TRIM(cname(jn)), 2.*rpi/domega(jn)/3600., dl_amp, dl_ph
     END DO
  ENDIF

  CALL write_file

CONTAINS

  SUBROUTINE read_file
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE read_file  ***
    !!
    !! ** Purpose :  Read input file 
    !!
    !! ** Method  :  Netcdf get var etc..
    !!
    !!----------------------------------------------------------------------
    ! Local variables :
    INTEGER(KIND=4) :: istatus
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dl_lon, dl_lat

    ! Read input file ( sanity check done in main)
    !----------------- 
    istatus = NF90_OPEN (cf_in, NF90_NOWRITE, ncid)
    istatus = NF90_INQ_DIMID (ncid,cdim_x,id_var)
    IF (istatus == NF90_NOERR) THEN 
       ln_moor=.FALSE.
    ELSE
       ln_moor=.TRUE.
    ENDIF

    istatus=NF90_INQ_DIMID (ncid,cdim_t,id_var)
    IF (istatus /= NF90_NOERR) THEN
       PRINT *, 'Pb with time-dimension in file: ', cf_in
       STOP
    ENDIF
    istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=npt)

    istatus = NF90_INQ_VARID(ncid,cv_time,id_var)
    IF (istatus /= NF90_NOERR) THEN
       PRINT *, 'Pb with time-variable in file: ', cf_in
       STOP
    ENDIF

    ALLOCATE(dg_time(npt))
    istatus = NF90_GET_VAR(ncid,id_var, dg_time)
    istatus = NF90_GET_ATT(ncid, id_var, "units", ca_units)
print *, ca_units
    IF ( index(ca_units,'hours') /= 0 ) THEN
       dcoef=24.d0
    ELSEIF ( index(ca_units,'seconds') /= 0 ) THEN
       dcoef=86400.d0
    ELSE 
       PRINT *, 'time unit not understood'
       STOP
    ENDIF
    READ(ca_units,7000)  iyy, imm, idd, ihh, imn, isec
7000 FORMAT('hours since ', I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)
!7000 FORMAT(a,a, I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)
    print *, iyy,imm,idd

    ijuli = julday(imm,idd,iyy)             ! Julian day of time origin in data file
    ijul0 = ijuli + FLOOR(dg_time(1)/dcoef) ! Julian day of first time dump in file
    !rbb
    PRINT*,dg_time(1)/dcoef, ijuli, ijul0

    CALL caldat(ijul0,imm,idd,iyy)          ! Set initial date month/day/year 
    immb=imm ; iddb=idd ; iyyb=iyy          ! save initial time (begin) for information
    ijulb=ijul0

    dl_hfrac=ihh+imn/60.d0+isec/3600.d0
    ijuli=ijul0

    DEALLOCATE(dg_time)

    IF (ln_moor) THEN
       ALLOCATE(dlon_r4(1,1), dlat_r4(1,1))
       istatus = NF90_INQ_VARID(ncid, cv_lon, id_var) ; istatus = NF90_GET_VAR(ncid, id_var, dlon_r4)
       istatus = NF90_INQ_VARID(ncid, cv_lat, id_var) ; istatus = NF90_GET_VAR(ncid, id_var, dlat_r4)
    ELSE
       istatus=NF90_INQ_DIMID (ncid, cdim_x, id_var)  ; istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=npi)
       istatus=NF90_INQ_DIMID (ncid, cdim_y, id_var)  ; istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=npj)
       IF ( ln_lonlat_2d_in ) THEN
          ALLOCATE( dlon_r4(npi,npj), dlat_r4(npi,npj) )
          istatus = NF90_INQ_VARID(ncid, cv_lon, id_var) ; istatus = NF90_GET_VAR(ncid, id_var, dlon_r4)
          istatus = NF90_INQ_VARID(ncid, cv_lat, id_var) ; istatus = NF90_GET_VAR(ncid, id_var, dlat_r4)
       ELSE
          ALLOCATE( dlon_r4(npi,1), dlat_r4(1,npj) ,dl_lon(npi), dl_lat(npj))
          istatus = NF90_INQ_VARID(ncid, cv_lon, id_var) ; istatus = NF90_GET_VAR(ncid, id_var, dl_lon)
          istatus = NF90_INQ_VARID(ncid, cv_lat, id_var) ; istatus = NF90_GET_VAR(ncid, id_var, dl_lat)
          dlon_r4(:,1)=dl_lon
          dlat_r4(1,:)=dl_lat
          DEALLOCATE (dl_lon, dl_lat )
       ENDIF
    ENDIF

    istatus = NF90_CLOSE(ncid)

  END SUBROUTINE read_file

  SUBROUTINE write_file
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE write_file  ***
    !!
    !! ** Purpose :  Create netcdf output file and fill it with tidal harmonics
    !!
    !! ** Method  :  NetCdf 
    !!
    !!----------------------------------------------------------------------
    ! local variables
    INTEGER(KIND=4) :: ji,jj
    INTEGER(KIND=4) :: istatus
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_tmp
    !!----------------------------------------------------------------------
    IF ( lnc4 ) THEN
       istatus = NF90_CREATE(cn_fharm, NF90_NETCDF4, ncid)
    ELSE
       istatus = NF90_CREATE(cn_fharm, NF90_64BIT_OFFSET, ncid)
    ENDIF

    IF (ln_moor) THEN
       npiout=1
       npjout=1
       istatus = NF90_DEF_DIM(ncid, cdim_x_out,      npiout, id_dim(1))
       istatus = NF90_DEF_VAR(ncid, cv_lon_out,      NF90_FLOAT, id_dim(1), id_var_lon)
       istatus = NF90_PUT_ATT(ncid, id_var_lon, 'long_name','Longitude')
       istatus = NF90_PUT_ATT(ncid, id_var_lon, 'units', 'degrees')
       istatus = NF90_DEF_DIM(ncid, cdim_y_out,      npjout, id_dim(2))
       istatus = NF90_DEF_VAR(ncid, cv_lat_out,      NF90_FLOAT, id_dim(2), id_var_lat)
       istatus = NF90_PUT_ATT(ncid, id_var_lat, 'long_name','Latitude')
       istatus = NF90_PUT_ATT(ncid, id_var_lat, 'units', 'degrees')
    ELSE
       npiout=npi
       npjout=npj
       istatus = NF90_DEF_DIM(ncid, cdim_x_out,     npiout, id_dim(1))
       istatus = NF90_DEF_DIM(ncid, cdim_y_out,     npjout, id_dim(2))
       IF ( ln_lonlat_2d_out ) THEN
         istatus = NF90_DEF_VAR(ncid, cv_lon_out,     NF90_FLOAT, id_dim(1:2), id_var_lon)
       ELSE
         istatus = NF90_DEF_VAR(ncid, cv_lon_out,     NF90_FLOAT, id_dim(1), id_var_lon)
       ENDIF
       istatus = NF90_PUT_ATT(ncid, id_var_lon, 'long_name','Longitude')
       istatus = NF90_PUT_ATT(ncid, id_var_lon, 'units', 'degrees')

       IF ( ln_lonlat_2d_out ) THEN
         istatus = NF90_DEF_VAR(ncid, cv_lat_out,     NF90_FLOAT, id_dim(1:2), id_var_lat)
       ELSE
         istatus = NF90_DEF_VAR(ncid, cv_lat_out,     NF90_FLOAT, id_dim(2), id_var_lat)
       ENDIF
       istatus = NF90_PUT_ATT(ncid, id_var_lat, 'long_name','Latitude')
       istatus = NF90_PUT_ATT(ncid, id_var_lat, 'units', 'degrees')
    ENDIF

    DO jn = 1, npconst
       istatus = NF90_DEF_VAR(ncid,TRIM(cname(jn))//TRIM(cn_v_out_x), NF90_FLOAT, id_dim(1:2), id_varhx(jn)) 
       istatus = NF90_PUT_ATT(ncid, id_varhx(jn), 'long_name',TRIM(cname(jn))//TRIM(cn_v_out_x) )
       istatus = NF90_PUT_ATT(ncid, id_varhx(jn), 'units', 'meters')
       istatus = NF90_PUT_ATT(ncid, id_varhx(jn), ca_miss_out, 0.)

       istatus = NF90_DEF_VAR(ncid,TRIM(cname(jn))//TRIM(cn_v_out_y), NF90_FLOAT, id_dim(1:2), id_varhy(jn)) 
       istatus = NF90_PUT_ATT(ncid, id_varhy(jn), 'long_name',TRIM(cname(jn))//TRIM(cn_v_out_y) )
       istatus = NF90_PUT_ATT(ncid, id_varhy(jn), 'units', 'meters')
       istatus = NF90_PUT_ATT(ncid, id_varhy(jn), ca_miss_out, 0.)

       istatus = NF90_DEF_VAR(ncid,TRIM(cname(jn))//'_A', NF90_FLOAT, id_dim(1:2), id_varA(jn))
       istatus = NF90_PUT_ATT(ncid, id_varA(jn), 'long_name',TRIM(cname(jn))//'_Amplitude')
       istatus = NF90_PUT_ATT(ncid, id_varA(jn), 'units', 'meters')
       istatus = NF90_PUT_ATT(ncid, id_varA(jn), ca_miss_out, 0.)

       istatus = NF90_DEF_VAR(ncid,TRIM(cname(jn))//'_G', NF90_FLOAT, id_dim(1:2), id_varG(jn))
       istatus = NF90_PUT_ATT(ncid, id_varG(jn), 'long_name',TRIM(cname(jn))//'_Phase')
       istatus = NF90_PUT_ATT(ncid, id_varG(jn), 'units', 'degrees')
       istatus = NF90_PUT_ATT(ncid, id_varG(jn), ca_miss_out, 0.)

    END DO

    istatus = NF90_ENDDEF(ncid)

    ! need to test if change of dimemsion between in and out
    IF ( ln_lonlat_2d_in  == ln_lonlat_2d_out ) THEN
      istatus = NF90_PUT_VAR(ncid, id_var_lon,dlon_r4)
      istatus = NF90_PUT_VAR(ncid, id_var_lat,dlat_r4)
    ELSE  ! 1d in and 2d out 
      ALLOCATE ( dl_tmp(npiout,npjout) )
      DO jj=1,npjout
        dl_tmp(:,jj) = dlon_r4(:,1)
      ENDDO
      istatus = NF90_PUT_VAR(ncid, id_var_lon,dl_tmp)
      DO ji=1,npiout
        dl_tmp(ji,:) = dlat_r4(1,:)
      ENDDO
      istatus = NF90_PUT_VAR(ncid, id_var_lat,dl_tmp)
      DEALLOCATE( dl_tmp )
    ENDIF
    DO jn = 1, npconst
       istatus = NF90_PUT_VAR(ncid, id_varhx(jn), dana_amp(:,:,jn,1))
       istatus = NF90_PUT_VAR(ncid, id_varhy(jn), dana_amp(:,:,jn,2))
       istatus = NF90_PUT_VAR(ncid, id_varA(jn),  dana_amp(:,:,jn,3))
       istatus = NF90_PUT_VAR(ncid, id_varG(jn),  dana_amp(:,:,jn,4))
    END DO

    istatus = NF90_CLOSE(ncid)

  END SUBROUTINE write_file


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


  SUBROUTINE GetFileList
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE GetFileList  ***
    !!
    !! ** Purpose :  Set up a file list given on the command line as 
    !!               blank separated list
    !!
    !! ** Method  :  Scan the command line until a '-' is found
    !!----------------------------------------------------------------------
    INTEGER (KIND=4)  :: ji
    INTEGER (KIND=4)  :: icur
    !!----------------------------------------------------------------------
    !!
    nfiles=0
    ! need to read a list of file ( number unknow ) 
    ! loop on argument till a '-' is found as first char
    icur=ijarg                          ! save current position of argument number
    DO ji = icur, narg                  ! scan arguments till - found
       CALL getarg ( ji, cldum )
       IF ( cldum(1:1) /= '-' ) THEN ; nfiles = nfiles+1
       ELSE                          ; EXIT
       ENDIF
    ENDDO
    ALLOCATE (cf_lst(nfiles) )
    DO ji = icur, icur + nfiles -1
       CALL getarg(ji, cf_lst( ji -icur +1 ) ) ; ijarg=ijarg+1
    END DO
  END SUBROUTINE GetFileList

  SUBROUTINE ZeroMeanInit (cd_cgrid)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ZeroMeanInit  ***
    !!
    !! ** Purpose :  Initialize the process of computing the mean value of
    !!               a 2D field 
    !!
    !! ** Method  :  Read mesh_hgr, mask file and computes the weights, fixed
    !!               in time :  aread2d(:,:) and areatot 
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cd_cgrid
    ! local variable
    INTEGER(KIND=4) :: icid, id, istatus
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dl_e1, dl_e2
    CHARACTER(LEN=20) :: cl_ve1, cl_ve2, cl_vmsk
    ! 
    !!----------------------------------------------------------------------
    ! choose variable name according to cd_cgrid  position
    SELECT CASE ( cd_cgrid )
    CASE ( 'T' )
       cl_ve1 = ' e1t'
       cl_ve2 = ' e2t'
       cl_vmsk= 'tmask'
    CASE ( 'U' )
       cl_ve1 = ' e1u'
       cl_ve2 = ' e2u'
       cl_vmsk= 'umask'
    CASE ( 'V' )
       cl_ve1 = ' e1v'
       cl_ve2 = ' e2v'
       cl_vmsk= 'vmask'
    CASE ( 'F' )
       cl_ve1 = ' e1f'
       cl_ve2 = ' e2f'
       cl_vmsk= 'fmask'
    END SELECT

    ! Allocate arrays
    ALLOCATE ( darea2d(npi,npj), mask(npi,npj) )
    ALLOCATE ( dl_e1(npi,npj), dl_e2(npi,npj) )
    ! read horizontal metrics
    istatus = NF90_OPEN(cf_hgr, NF90_NOWRITE, icid)

    istatus = NF90_INQ_VARID(icid, cl_ve1, id ) 
    istatus = NF90_GET_VAR(icid, id, dl_e1, start=(/1,1,1/), count=(/npi,npj,1/) )

    istatus = NF90_INQ_VARID(icid, cl_ve2, id ) 
    istatus = NF90_GET_VAR(icid, id, dl_e2, start=(/1,1,1/), count=(/npi,npj,1/) )
    istatus = NF90_CLOSE(icid)

    ! read mask
    istatus = NF90_OPEN(cf_msk, NF90_NOWRITE, icid)

    istatus = NF90_INQ_VARID(icid, cl_vmsk, id ) 
    istatus = NF90_GET_VAR(icid, id, mask, start=(/1,1,1/), count=(/npi,npj,1/) )
    istatus = NF90_CLOSE(icid)

    ! compute weights
    darea2d(:,:) = dl_e1(:,:) * dl_e2(:,:) * mask(:,:)
    dareatot     = SUM( darea2d )

   DEALLOCATE (dl_e1, dl_e2 )
  END SUBROUTINE ZeroMeanInit

  SUBROUTINE ZeroMean(dd_tab)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ZeroMean  ***
    !!
    !! ** Purpose :  Set the mean value of the input/output array to zero
    !!
    !! ** Method  :  Ccompute mean value (weighted) and subtract it from the field 
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=8), DIMENSION(:,:), INTENT(inout) :: dd_tab
    ! local variables
    REAL(KIND=8)                                :: dl_mean
    !!----------------------------------------------------------------------
    dl_mean =  SUM( dd_tab(:,:)*darea2d(:,:) ) / dareatot
    dd_tab(:,:) = (dd_tab(:,:) - dl_mean ) * mask(:,:) 
  END SUBROUTINE ZeroMean

  LOGICAL FUNCTION chkfile (cd_file, ld_verbose )
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION chkfile  ***
    !!
    !! ** Purpose :  Check if cd_file exists.
    !!               Return false if it exists, true if it does not
    !!               Do nothing is filename is 'none'
    !!
    !! ** Method  : Doing it this way allow statements such as
    !!              IF ( chkfile( cf_toto) ) STOP 99  ! missing file
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*),  INTENT(in) :: cd_file
    LOGICAL, OPTIONAL, INTENT(in) :: ld_verbose

    INTEGER(KIND=4)               :: ierr
    LOGICAL                       :: ll_exist, ll_verbose
    !!----------------------------------------------------------------------
    IF ( PRESENT(ld_verbose) ) THEN
       ll_verbose = ld_verbose
    ELSE
       ll_verbose = .TRUE.
    ENDIF

    IF ( TRIM(cd_file) /= 'none')  THEN
       INQUIRE (file = TRIM(cd_file), EXIST=ll_exist)

       IF (ll_exist) THEN
          chkfile = .false.
       ELSE
          IF ( ll_verbose ) PRINT *, ' File ',TRIM(cd_file),' is missing '
          chkfile = .true.
       ENDIF
    ELSE
       chkfile = .false.  ! 'none' file is not checked
    ENDIF

  END FUNCTION chkfile

END PROGRAM tid_harm_ana
