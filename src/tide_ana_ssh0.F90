   PROGRAM tide_ana
 
   USE tide
   USE surdetermine
   USE netcdf

   IMPLICIT NONE

   LOGICAL :: moor, l_exist, ln_short
   INTEGER, PARAMETER :: nconst=20
   REAL(8) :: tsamp ! Sampling period in seconds
   CHARACTER(len=20) :: name_var, name_var_out_x, name_var_out_y
   INTEGER :: jul0, juli, jule
   INTEGER :: narg, iargc, init, xtype
   INTEGER, DIMENSION(nconst) :: id_varhx, id_varhy

   REAL(8) :: scale_factor =1., add_offset =0.,missval=1.e20
   REAL(8) :: ztemp, ztime, zamp, zph, HFRAC, r_ave
   REAL(8) :: zmean, areatot
   REAL(8), DIMENSION(nconst) :: omega, ft, ut, vt, vt0
   REAL(8), DIMENSION(:), ALLOCATABLE :: ig_time
   REAL(8), DIMENSION(:,:), ALLOCATABLE :: temp2d_r4, lon_r4, lat_r4
   REAL(8), DIMENSION(:,:), ALLOCATABLE :: area2d, e1t, e2t, mask 
   REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: tmp_r4
   REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: ana_temp
   REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: ana_amp
   INTEGER(KIND=2), DIMENSION(:,:,:), ALLOCATABLE :: int2tabglo

   CHARACTER(len=4), DIMENSION(jpmax_harmo) :: tname
   CHARACTER(len=120) :: ncfile_harm, ncfile_in
   CHARACTER(len=33) :: units, tmpch
   INTEGER :: t, jpi, jpj, nt, ncid, istatus, i0, j0, n, id_var, id_vart, id_varx, id_vary, jpiout, jpjout, jpih, jpjh
   INTEGER, DIMENSION(3) :: id_dim
   INTEGER :: yy, mm, dd, hh, ig_min, sec, ji, jj, nhc, jc, nhan, ksp, kun, yye, mme, dde, kfile, id_var_lon, id_var_lat

   NAMELIST / analysis_param / moor, &
                               ln_short, &
                               name_var, &
                               name_var_out_x, &
                               name_var_out_y, &
                               tsamp, &
                               r_ave, &
                               ncfile_harm

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
   PRINT *, 'Variable to be analysed = ', TRIM(name_var)
   PRINT *, 'Sampling period = ', tsamp
   PRINT *, 'Output file name = ', TRIM(ncfile_harm)
   PRINT *, ' '


         tname(1 )='M2'
         tname(2 )='N2'
         tname(3 )='S2'
         tname(4 )='K2'
         tname(5 )='Q1'
         tname(6 )='K1'
         tname(7 )='O1'
         tname(8 )='P1'
         tname(9 )='Mf'
         tname(10)='Mm'
         tname(11)='M4'
         tname(12)='M6'
         tname(13)='S1'
         tname(14)='Msqm'
         tname(15)='Mtm'
         tname(16)='_2N2'
         tname(17)='MU2'
         tname(18)='NU2'
         tname(19)='L2'
         tname(20)='T2'

         narg = iargc()
         IF (narg<1) THEN
            PRINT *, 'Need at least one input file'
            PRINT *, 'tide_ana usage:' 
            PRINT *, './tide_ana file_list'
            PRINT *, 'EXEMPLE: ./tide_ana BISCAY-T05-y2008*_HF_gridT.nc'
            PRINT *, 'Output harmonic file in res_harm.nc'
            STOP
         ENDIF 

         CALL getarg (1, ncfile_in)

         CALL read_file

         ! Set tidal waves pulsations and nodal corrections (at initial time)
         !-------------------------------------------------------------------
         CALL tide_pulse( omega(1:nconst), tname(1:nconst) ,nconst)
         CALL tide_vuf( vt0(1:nconst), ut(1:nconst) , ft(1:nconst), tname(1:nconst) ,nconst, yy, mm, dd, HFRAC)
 
         ! Perform harmonic analysis
         !--------------------------

         IF (moor) THEN
           ALLOCATE(tmp_r4(1,1,1))
           ALLOCATE(ana_temp(1,1,2*nconst))
           ALLOCATE(ana_amp(1,1,nconst,2))
           jpiout=1
           jpjout=1
         ELSE
           ALLOCATE(tmp_r4(jpi,jpj,1))
           ALLOCATE(ana_temp(jpi,jpj,2*nconst))
           ALLOCATE(ana_amp(jpi,jpj,nconst,2))
           jpiout=jpi
           jpjout=jpj
         ENDIF
         ALLOCATE( area2d(jpi,jpj), e1t(jpi,jpj), e2t(jpi,jpj), mask(jpi,jpj) ) 

         istatus = NF90_OPEN ('mesh_hgr.nc', NF90_NOWRITE, ncid)
         IF (istatus /= NF90_NOERR) THEN
           PRINT *, 'Can not open file: mesh_hgr.nc'
           STOP
         ELSE
           WRITE(6,*)''
           WRITE(6,'(a,2x,a)') 'OPEN FILE: mesh_hgr.nc'
           CALL FLUSH(6)
         ENDIF
         istatus = NF90_INQ_VARID(ncid,'e1t',id_var)
         istatus = NF90_GET_VAR(ncid,id_var, e1t, &
              start = (/ 1, 1, 1 /),count = (/ jpiout, jpjout, 1 /))
         istatus = NF90_INQ_VARID(ncid,'e2t',id_var)
         istatus = NF90_GET_VAR(ncid,id_var, e2t, &
              start = (/ 1, 1, 1 /),count = (/ jpiout, jpjout, 1 /))
         istatus = NF90_CLOSE(ncid)

         istatus = NF90_OPEN ('mask.nc', NF90_NOWRITE, ncid)
         IF (istatus /= NF90_NOERR) THEN
           PRINT *, 'Can not open file: mask.nc'
           STOP
         ELSE
           WRITE(6,*)''
           WRITE(6,'(a,2x,a)') 'OPEN FILE: mask.nc'
           CALL FLUSH(6)
         ENDIF
         istatus = NF90_INQ_VARID(ncid,'tmask',id_var)
         istatus = NF90_GET_VAR(ncid,id_var, mask, &
              start = (/ 1, 1, 1 /),count = (/ jpiout, jpjout, 1 /))
         istatus = NF90_CLOSE(ncid)

         area2d(:,:) = e1t(:,:) * e2t(:,:) * mask(:,:)
         areatot     = SUM( area2d(:,:) )

         ana_temp(:,:,:) = 0.e0
         ztime = r_ave * tsamp
         nhan  = 0  ! Total number of time dumps
         ksp   = 0  ! Initialize length of sparse matrix

!MBK
         write(6,*)
         write(6,*)' ztime : ',ztime
         write(6,*)' r_ave : ',r_ave
         write(6,*)' tsamp : ',tsamp
         write(6,*)
         CALL FLUSH(6)
!MBK
         DO kfile=1,narg
         CALL getarg (kfile, ncfile_in)
         ! START LOOP ON INPUT FILES HERE
         ! open file again
         istatus = NF90_OPEN (ncfile_in, NF90_NOWRITE, ncid)
         IF (istatus /= NF90_NOERR) THEN
           PRINT *, 'Can not open file: ', ncfile_in
           STOP
         ELSE
           WRITE(6,*)''
           WRITE(6,'(a,2x,a)') 'OPEN FILE: ',TRIM(ncfile_in)
           CALL FLUSH(6)
         ENDIF

         ! get number of time dumps
         istatus=NF90_INQ_DIMID (ncid,'time',id_var)
         IF (istatus /= NF90_NOERR) THEN
            istatus=NF90_INQ_DIMID (ncid,'time_counter',id_var)
            IF (istatus /= NF90_NOERR) THEN
              PRINT *, 'Pb with time-dimension in file: ', ncfile_in
              STOP
            ENDIF
         ENDIF
         istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=nt)
         PRINT *, 'NUMBER OF TIME DUMPS TO READ IN FILE: ', nt

         istatus = NF90_INQ_VARID(ncid,TRIM(name_var),id_var)
         IF (istatus /= NF90_NOERR) THEN
            PRINT *, 'CAN NOT FIND VARIABLE: ', TRIM(name_var), ' IN FILE ', ncfile_in
            STOP
         ELSE
            istatus = NF90_GET_ATT(ncid, id_var, "missing_value", missval)
            istatus = NF90_INQUIRE_VARIABLE(ncid, id_var, tmpch, xtype)
            IF (xtype==NF90_SHORT) THEN 
               ln_short=.TRUE.
               istatus = NF90_GET_ATT(ncid, id_var, "add_offset", add_offset)
               istatus = NF90_GET_ATT(ncid, id_var, "scale_factor", scale_factor)
               IF (kfile==1) THEN
                  IF (moor)      ALLOCATE(int2tabglo(1,1,1))
                  IF (.NOT.moor) ALLOCATE(int2tabglo(jpi,jpj,1))
               ENDIF
            ENDIF    
         ENDIF

         DO t=1, nt
         !  Get variable
            IF (ln_short) THEN
              istatus = NF90_GET_VAR(ncid,id_var, int2tabglo, &
              start = (/ 1, 1, t /),count = (/ jpiout, jpjout, 1 /))
              WHERE(int2tabglo .NE. missval) tmp_r4 = int2tabglo*scale_factor+add_offset
              WHERE(int2tabglo .EQ. missval) tmp_r4 = 0.
            ELSE
              istatus = NF90_GET_VAR(ncid,id_var, tmp_r4, &
              start = (/ 1, 1, t /),count = (/ jpiout, jpjout, 1 /))
              WHERE(tmp_r4 .EQ. missval) tmp_r4 = 0.
!rbb
              zmean = SUM( tmp_r4(:,:,1)*area2d(:,:) ) / areatot
              write(6,*)' mean_2D :', zmean
              tmp_r4(:,:,1) = ( tmp_r4(:,:,1) - zmean ) * mask(:,:)
!rbb end
            ENDIF

            ztime = ztime + tsamp

            ! Update nodal corrections:
            jul0 = juli + FLOOR(ztime/86400.) ! Julian day of first time dump in file
            CALL caldat(jul0,mm,dd,yy)
            HFRAC = juli+ztime/86400.-jul0

         write(6,*)' mm dd yy :',mm,dd,yy
         write(6,'(a,1x,f12.2,1x,a,1x,f12.2,1x,a,1x,f18.8)')' ztime    : ',ztime,'  jul0 : ',jul0,' HFRAC : ', HFRAC
         CALL FLUSH(6)

            CALL tide_vuf( vt(1:nconst), ut(1:nconst) , ft(1:nconst), tname(1:nconst) ,nconst, yy, mm, dd, HFRAC)

            nhan = nhan + 1
            nhc = 0
            DO n = 1, nconst
               DO jc = 1,2
                  nhc = nhc + 1
                  ksp = ksp + 1
                  ztemp =(     MOD(jc,2) * ft(n) *COS(omega(n)*ztime + vt0(n) + ut(n))  &
                          +(1.-MOD(jc,2))* ft(n) *SIN(omega(n)*ztime + vt0(n) + ut(n)))
                  ana_temp(:,:,nhc) = ana_temp(:,:,nhc) + ztemp*tmp_r4(:,:,1)                   
                  ISPARSE(ksp) = nhan
                  JSPARSE(ksp) = nhc
                  SPARSEVALUE(ksp) = ztemp
               END DO
            END DO
         END DO
         istatus = NF90_CLOSE(ncid)
         ! STOP LOOP ON INPUT FILES HERE
         END DO


         jule = jul0 + FLOOR(ztime/86400.) - 1
         CALL caldat(jule,mme,dde,yye)  

         print *,''
         print *,'START DATE (ddmmyy)        : ', dd, mm, yy
         print *,'END DATE (ddmmyy)          : ', dde, mme, yye
         print *,'NUMBER OF DAYS FOR ANALYSIS: ', jule-jul0 + 1
         print *,'SAMPLING PERIOD (seconds)  : ', FLOOR(tsamp)
         print *,''

         ! Matrix inversion
         NBINCO   = 2*nconst
         NBSPARSE = ksp
         CALL SUR_DETERMINE(1)

         ! Find solution for each point
         DO jj = 1, jpjout
            DO ji = 1, jpiout
               kun=0
               DO n = 1, nconst
                 DO jc = 1,2
                    kun = kun + 1
                    TAB4(kun)=ana_temp(ji,jj,kun)
                 ENDDO
               ENDDO

               CALL SUR_DETERMINE(0)

               ! Fill output array
               DO n = 1, nconst
                 ana_amp(ji,jj,n,1)= TAB7((n-1)*2+1)
                 ana_amp(ji,jj,n,2)=-TAB7((n-1)*2+2)
               END DO
            END DO
         END DO

         DEALLOCATE(ana_temp)

         ! WRITE output
         !-------------
        
         IF (moor) THEN
            DO n = 1, nconst
              zamp = SQRT(ana_amp(1,1,n,1)**2+ana_amp(1,1,n,2)**2) ! Amplitude
              zph  = ATAN2(-ana_amp(1,1,n,2)/zamp,ana_amp(1,1,n,1)/zamp)/rad ! Phase
              PRINT *, 'constituent name / Period (hour) / amplitude / phase (deg)'
              PRINT *, TRIM(tname(n)), 2.*rpi/omega(n)/3600., zamp, zph
            END DO
         ENDIF

         CALL write_file

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE read_file

         ! Read input file
         !----------------- 
         istatus = NF90_OPEN (ncfile_in, NF90_NOWRITE, ncid)
         IF (istatus /= NF90_NOERR) THEN
           PRINT *, 'Can not open file: ', ncfile_in
           STOP
         ENDIF

         istatus=NF90_INQ_DIMID (ncid,'x',id_var)
         IF (istatus == NF90_NOERR) THEN 
            moor=.FALSE.
         ELSE
            moor=.TRUE.
         ENDIF

         ! Get mooring indice on global grid

         istatus=NF90_INQ_DIMID (ncid,'time',id_var)
         IF (istatus /= NF90_NOERR) THEN
            istatus=NF90_INQ_DIMID (ncid,'time_counter',id_var)
            IF (istatus /= NF90_NOERR) THEN
              PRINT *, 'Pb with time-dimension in file: ', ncfile_in
              STOP
            ENDIF
         ENDIF
         istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=nt)

         istatus = NF90_INQ_VARID(ncid,'time',id_var)
         IF (istatus /= NF90_NOERR) THEN
           istatus = NF90_INQ_VARID(ncid,'time_counter',id_var)
            IF (istatus /= NF90_NOERR) THEN
              PRINT *, 'Pb with time-dimension in file: ', ncfile_in
              STOP
            ENDIF
         ENDIF
         ALLOCATE(ig_time(nt))
         istatus = NF90_GET_VAR(ncid,id_var, ig_time)
         istatus = NF90_GET_ATT(ncid, id_var, "units", units)
         READ(units,7000) yy, mm, dd, hh, ig_min, sec
7000 FORMAT('seconds since ', I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)

         juli = julday(mm,dd,yy)             ! Julian day of time origin in data file
         jul0 = juli + FLOOR(ig_time(1)/86400.) ! Julian day of first time dump in file
!rbb
print*,ig_time(1)/86400.

         CALL caldat(jul0,mm,dd,yy)          ! Set initial date month/day/year 
         HFRAC=hh+ig_min/60.+sec/3600.
         juli=jul0

         DEALLOCATE(ig_time)

         IF (moor) THEN
           ALLOCATE(lon_r4(1,1), lat_r4(1,1))
           istatus = NF90_INQ_VARID(ncid,'longitude',id_var)
           istatus = NF90_GET_VAR(ncid,id_var, lon_r4)
           istatus = NF90_INQ_VARID(ncid,'latitude',id_var)
           istatus = NF90_GET_VAR(ncid,id_var, lat_r4)
         ELSE
           istatus=NF90_INQ_DIMID (ncid,'x',id_var)
           istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=jpi)
           istatus=NF90_INQ_DIMID (ncid,'y',id_var)
           istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=jpj)
           ALLOCATE(lon_r4(jpi,jpj), lat_r4(jpi,jpj))
           istatus = NF90_INQ_VARID(ncid,'nav_lon',id_var)
           istatus = NF90_GET_VAR(ncid,id_var, lon_r4)
           istatus = NF90_INQ_VARID(ncid,'nav_lat',id_var)
           istatus = NF90_GET_VAR(ncid,id_var, lat_r4)
         ENDIF
        
         istatus = NF90_CLOSE(ncid)

END SUBROUTINE read_file

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE write_file

         istatus = NF90_CREATE(ncfile_harm, NF90_CLOBBER, ncid)

         IF (moor) THEN
            jpiout=1
            jpjout=1
            istatus = NF90_DEF_DIM(ncid,'longitude',jpiout, id_dim(1))
            istatus = NF90_DEF_VAR(ncid,'longitude', NF90_FLOAT, id_dim(1), id_var_lon)
            istatus = NF90_PUT_ATT(ncid, id_var_lon, 'long_name','Longitude')
            istatus = NF90_PUT_ATT(ncid, id_var_lon, 'units', 'degrees')
            istatus = NF90_DEF_DIM(ncid,'latitude', jpjout, id_dim(2))
            istatus = NF90_DEF_VAR(ncid,'latitude', NF90_FLOAT, id_dim(2), id_var_lat)
            istatus = NF90_PUT_ATT(ncid, id_var_lat, 'long_name','latitude')
            istatus = NF90_PUT_ATT(ncid, id_var_lat, 'units', 'degrees')
         ELSE
            jpiout=jpi
            jpjout=jpj
            istatus = NF90_DEF_DIM(ncid,'x',jpiout, id_dim(1))
            istatus = NF90_DEF_DIM(ncid,'y',jpjout, id_dim(2))
            istatus = NF90_DEF_VAR(ncid,'nav_lon', NF90_FLOAT, id_dim(1:2), id_var_lon)
            istatus = NF90_PUT_ATT(ncid, id_var_lon, 'long_name','Longitude')
            istatus = NF90_PUT_ATT(ncid, id_var_lon, 'units', 'degrees')
            istatus = NF90_DEF_VAR(ncid,'nav_lat', NF90_FLOAT, id_dim(1:2), id_var_lat)
            istatus = NF90_PUT_ATT(ncid, id_var_lat, 'long_name','latitude')
            istatus = NF90_PUT_ATT(ncid, id_var_lat, 'units', 'degrees')
         ENDIF

         DO n = 1, nconst
            istatus = NF90_DEF_VAR(ncid,TRIM(tname(n))//TRIM(name_var_out_x), NF90_FLOAT, id_dim(1:2), id_varhx(n)) 
            istatus = NF90_PUT_ATT(ncid, id_varhx(n), 'long_name',TRIM(tname(n))//TRIM(name_var_out_x) )
            istatus = NF90_PUT_ATT(ncid, id_varhx(n), 'units', 'meters')
            istatus = NF90_PUT_ATT(ncid, id_varhx(n), 'missing_value', 0.)

            istatus = NF90_DEF_VAR(ncid,TRIM(tname(n))//TRIM(name_var_out_y), NF90_FLOAT, id_dim(1:2), id_varhy(n)) 
            istatus = NF90_PUT_ATT(ncid, id_varhy(n), 'long_name',TRIM(tname(n))//TRIM(name_var_out_y) )
            istatus = NF90_PUT_ATT(ncid, id_varhy(n), 'units', 'meters')
            istatus = NF90_PUT_ATT(ncid, id_varhy(n), 'missing_value', 0.)
         END DO

         istatus = NF90_ENDDEF(ncid)

         istatus = NF90_PUT_VAR(ncid, id_var_lon,lon_r4)
         istatus = NF90_PUT_VAR(ncid, id_var_lat,lat_r4)
         DO n = 1, nconst
            istatus = NF90_PUT_VAR(ncid, id_varhx(n),ana_amp(:,:,n,1))
            istatus = NF90_PUT_VAR(ncid, id_varhy(n),ana_amp(:,:,n,2))
         END DO

         istatus = NF90_CLOSE(ncid)

END SUBROUTINE write_file


      function julday(kmm,kid,kiyyy)
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
      IF (il_kiyyy.eq.0) stop 101
!
      IF (il_kiyyy.lt.0) il_kiyyy = il_kiyyy + 1
      IF (kmm.gt.2) THEN
         iy = il_kiyyy
         im = kmm + 1
      ELSE
         iy = il_kiyyy - 1
         im = kmm + 13
      END IF
!
      julday = int(365.25*iy) + int(30.6001*im) + kid + 1720995 
      IF (kid+31*(kmm+12*il_kiyyy).ge.jpgreg) THEN
         ia = int(0.01*iy)
         julday = julday + 2 - ia + int(0.25*ia) 
      END IF
      RETURN
      END FUNCTION julday
!###
      subroutine caldat(kjulian,kmm,kid,kiyyy)
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

      END PROGRAM
