         PROGRAM tide_predict
 
         USE tide
         USE netcdf

         IMPLICIT NONE

         LOGICAL :: moor=.TRUE., l_exist
         INTEGER, PARAMETER :: nconst=11
         CHARACTER(len=20)  :: name_var
         CHARACTER(len=255) :: long_name_var, var_units
         INTEGER :: julday, jul0, juli, t0
         REAL(8), DIMENSION(nconst) :: omega, ft, ut, vt, vt0
         REAL(8), DIMENSION(:), ALLOCATABLE :: rg_time
         REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: res
         REAL(8), DIMENSION(:,:), ALLOCATABLE :: temp2d_r4, lon_r4, lat_r4
         REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: tempx, tempy

         REAL(8) :: zamp, zph, HFRAC
         CHARACTER(len=4), DIMENSION(jpmax_harmo) :: tname
         CHARACTER(len=120) :: ncfile_harm, ncfile_in, ncfile_out  
         CHARACTER(len=33) :: ttout, units     
         INTEGER :: narg, iargc
         INTEGER :: t, jpi, jpj, nt, ncid, istatus, i0, j0, n, id_var, id_vart, id_varx, id_vary, jpiout, jpjout, jpih, jpjh
         INTEGER, DIMENSION(3) :: id_dim
         INTEGER :: yy, mm, dd, hh, ig_min, sec, ji, jj

   NAMELIST / prediction_param / name_var, long_name_var, var_units,  ncfile_harm, ncfile_out

   ! 0. Read the namelist
   ! -------------------------------
   INQUIRE( FILE = 'namelist', EXIST = l_exist )

   IF ( .NOT. l_exist ) STOP 'Namelist not found'

   OPEN( 20, file = 'namelist',status = 'old' ,  form   = 'formatted' )
   READ( 20, NML = prediction_param )
   CLOSE( 20 )

   PRINT *, ' '
   PRINT *, 'Variable to be predicted = ', TRIM(name_var)
   PRINT *, 'Variable long name = ', TRIM(long_name_var)
   PRINT *, 'Variable units = ', TRIM(var_units)
   PRINT *, 'Harmonic data file = ', TRIM(ncfile_harm)
   PRINT *, 'Prediction output file = ', TRIM(ncfile_out)
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
!         tname(12)='M6'

         narg = iargc()

         IF (narg<1) THEN
            PRINT *, 'WRONG NUMBER OF ARGUMENTS: need at least 1'
            PRINT *, 'tide_predict usage:'
            PRINT *, 'Compute Sea level tidal prediction corresponding to a given model file'
            PRINT *, './tide_predict file_in'
            PRINT *, 'file_in:   input netcdf file (mooring or 2d)'
            STOP
         ENDIF

         CALL getarg (1, ncfile_in)


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
         IF (moor) THEN
            istatus = NF90_GET_ATT(ncid, NF90_GLOBAL, "i-indice", i0)
            IF (istatus /= NF90_NOERR) THEN
              PRINT *, 'Can not find mooring i-indice global attribute in: ', ncfile_in
!              STOP
            ENDIF
            istatus = NF90_GET_ATT(ncid, NF90_GLOBAL, "j-indice", j0)
            IF (istatus /= NF90_NOERR) THEN
              PRINT *, 'Can not find mooring j-indice global attribute in: ', ncfile_in
!              STOP
            ENDIF
         ENDIF

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
         ALLOCATE(rg_time(nt))
         istatus = NF90_GET_VAR(ncid,id_var, rg_time)
         istatus = NF90_GET_ATT(ncid, id_var, "units", units)
         READ(units,7000) yy, mm, dd, hh, ig_min, sec
7000 FORMAT('seconds since ', I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)

         juli = julday(mm,dd,yy)             ! Julian day of time origin in data file
         jul0 = juli + FLOOR(rg_time(1)/86400.) ! Julian day of first time dump in file
         HFRAC=hh+ig_min/60.+sec/3600.
         CALL caldat(jul0,mm,dd,yy)
!
          write(6,*)
          write(6,*)' juli jul0 HFRAC ',juli,jul0,HFRAC
          write(6,*) mm,dd,yy
          write(6,*)
!
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

         WRITE(6,*)'INPUT FILE START DATE (ddmmyy) ', dd, mm, yy
         WRITE(6,*) 'Number OF time dumps ', nt
         WRITE(6,*) 'Number OF harmonics  ', nconst
         IF (moor) THEN
            WRITE(6,*)'Mooring file indices ', i0, j0
         ENDIF
         CALL flush(6)

         istatus = NF90_CLOSE(ncid)

         ! Set tidal waves pulsations and nodal corrections (at initial time)
         !-------------------------------------------------------------------
         CALL tide_pulse( omega(1:nconst), tname(1:nconst) ,nconst)
         CALL tide_vuf( vt(1:nconst), ut(1:nconst) , ft(1:nconst), tname(1:nconst) ,nconst, yy, mm, dd, HFRAC)
         vt0 =vt
 
         ! Read harmonic data file
         !-------------------------      
         istatus = NF90_OPEN (ncfile_harm, NF90_NOWRITE, ncid)
         IF (istatus /= NF90_NOERR) THEN
           PRINT *, 'Can not open file: ', ncfile_harm
           STOP
         ENDIF

         istatus=NF90_INQ_DIMID (ncid,'x',id_var)
         IF (istatus /= NF90_NOERR) THEN
           istatus=NF90_INQ_DIMID (ncid,'longitude',id_var)
           IF (istatus /= NF90_NOERR) THEN
              PRINT *, 'Pb with x-dimension in file: ', ncfile_harm
              STOP
           ENDIF
         ENDIF
         istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=jpih)

         istatus = NF90_INQ_DIMID (ncid,'y',id_var)
         IF (istatus /= NF90_NOERR) THEN
           istatus=NF90_INQ_DIMID (ncid,'latitude',id_var)
           IF (istatus /= NF90_NOERR) THEN
              PRINT *, 'Pb with y-dimension in file: ', ncfile_harm
              STOP
           ENDIF
         ENDIF
         istatus = NF90_INQUIRE_DIMENSION(ncid, id_var, len=jpjh)

         ALLOCATE(temp2d_r4(jpih,jpjh))
         IF (moor) THEN
            ALLOCATE(tempx(1,1,nconst), tempy(1,1,nconst))
         ELSE
            jpi=jpih
            jpj=jpjh
            ALLOCATE(tempx(jpi,jpj,nconst), tempy(jpi,jpj,nconst))
         ENDIF

         DO t=1, nconst
           istatus = NF90_INQ_VARID( ncid,TRIM(tname(t))//'_x_elev',id_var)
           IF (istatus /= NF90_NOERR) THEN
             PRINT *, 'Can not find harmonic data ', TRIM(tname(t))//'_x_elev in file: ', ncfile_harm
             STOP
           ENDIF
           istatus = NF90_GET_VAR(ncid,id_var, temp2d_r4)
           
           IF (moor) THEN 
              IF ((jpih==1).AND.(jpjh==1)) THEN   
                 tempx(1,1,t)=temp2d_r4(1,1)
              ELSE       
                 tempx(1,1,t)=temp2d_r4(i0,j0)
              ENDIF
           ELSE
              tempx(:,:,t)=temp2d_r4(:,:)
           ENDIF

           istatus = NF90_INQ_VARID( ncid,TRIM(tname(t))//'_y_elev',id_var)
           IF (istatus /= NF90_NOERR) THEN
             PRINT *, 'Can not find harmonic data ', TRIM(tname(t))//'_y_elev in file: ', ncfile_harm
             STOP
           ENDIF
           istatus = NF90_GET_VAR(ncid,id_var, temp2d_r4)
           IF (moor) THEN     
              IF ((jpih==1).AND.(jpjh==1)) THEN  
                 tempy(1,1,t)=temp2d_r4(1,1)
              ELSE 
                 tempy(1,1,t)=temp2d_r4(i0,j0)
              ENDIF
           ELSE
              tempy(:,:,t)=temp2d_r4(:,:)
           ENDIF
         END DO
         istatus = NF90_CLOSE(ncid)
         DEALLOCATE(temp2d_r4)

         ! Build residual
         !---------------

         t0 = FLOOR(rg_time(1)/86400.)*86400.
         istatus = NF90_CREATE(ncfile_out, NF90_CLOBBER, ncid)
         IF (moor) THEN
            jpiout=1
            jpjout=1
            istatus = nf90_DEF_DIM(ncid,'time', NF90_UNLIMITED, id_dim(3))
            istatus = NF90_DEF_VAR(ncid,'time', NF90_FLOAT, id_dim(3), id_vart)
            istatus = NF90_PUT_ATT(ncid, id_vart, 'long_name','Time')
            istatus = NF90_PUT_ATT(ncid, id_vart, 'units', units)
            istatus = nf90_DEF_DIM(ncid,'longitude',jpiout, id_dim(1))
            istatus = NF90_DEF_VAR(ncid,'longitude', NF90_FLOAT, id_dim(1), id_varx)
            istatus = NF90_PUT_ATT(ncid, id_varx, 'long_name','Longitude')
            istatus = NF90_PUT_ATT(ncid, id_varx, 'units', 'degrees')
            istatus = nf90_DEF_DIM(ncid,'latitude', jpjout, id_dim(2))
            istatus = NF90_DEF_VAR(ncid,'latitude', NF90_FLOAT, id_dim(2), id_vary)
            istatus = NF90_PUT_ATT(ncid, id_vary, 'long_name','latitude')
            istatus = NF90_PUT_ATT(ncid, id_vary, 'units', 'degrees')
         ELSE
            jpiout=jpi
            jpjout=jpj
            istatus = nf90_DEF_DIM(ncid,'time_counter', NF90_UNLIMITED, id_dim(3))
            istatus = NF90_DEF_VAR(ncid,'time_counter', NF90_FLOAT, id_dim(3), id_vart)
            istatus = NF90_PUT_ATT(ncid, id_vart, 'long_name','Time')
            istatus = NF90_PUT_ATT(ncid, id_vart, 'units', units)
            istatus = nf90_DEF_DIM(ncid,'x',jpiout, id_dim(1))
            istatus = nf90_DEF_DIM(ncid,'y',jpjout, id_dim(2))
            istatus = NF90_DEF_VAR(ncid,'nav_lon', NF90_FLOAT, id_dim(1:2), id_varx)
            istatus = NF90_PUT_ATT(ncid, id_varx, 'long_name','Longitude')
            istatus = NF90_PUT_ATT(ncid, id_varx, 'units', 'degrees')
            istatus = NF90_DEF_VAR(ncid,'nav_lat', NF90_FLOAT, id_dim(1:2), id_vary)
            istatus = NF90_PUT_ATT(ncid, id_vary, 'long_name','latitude')
            istatus = NF90_PUT_ATT(ncid, id_vary, 'units', 'degrees')

         ENDIF
         istatus = nf90_DEF_VAR(ncid,TRIM(name_var), NF90_FLOAT, id_dim(1:3), id_var) 
         istatus = NF90_PUT_ATT(ncid, id_var, 'long_name', TRIM(long_name_var))
         istatus = NF90_PUT_ATT(ncid, id_var, 'units', TRIM(var_units))
         istatus = NF90_PUT_ATT(ncid, id_var, 'missing_value', 0.)
         istatus = NF90_ENDDEF(ncid)

         ALLOCATE(res(jpiout,jpjout,1))

         istatus = NF90_PUT_VAR(ncid, id_vart,rg_time(1:nt))
         istatus = NF90_PUT_VAR(ncid, id_varx,lon_r4)
         istatus = NF90_PUT_VAR(ncid, id_vary,lat_r4)

         DO t=1,nt 

            jul0 = juli + FLOOR(rg_time(t)/86400.) ! Julian day of first time dump in file
            CALL caldat(jul0,mm,dd,yy)
            HFRAC = juli+rg_time(t)/86400.-jul0
!
!          write(6,*)
!         write(6,*)' time ',rg_time(t)
!         write(6,*)' juli jul0 HFRAC ',juli,jul0,HFRAC
!         write(6,*) mm,dd,yy
!         write(6,*)
!
            CALL tide_vuf( vt(1:nconst), ut(1:nconst) , ft(1:nconst), tname(1:nconst) ,nconst, yy, mm, dd, HFRAC)
            res(:,:,:) = 0.    
            DO n=1, nconst
               DO ji=1,jpiout
                  DO jj=1,jpjout
                     zamp = SQRT(tempx(ji,jj,n)**2 + tempy(ji,jj,n)**2)
                     zph = ATAN2(tempy(ji,jj,n),tempx(ji,jj,n)) + vt0(n) + ut(n)
                     res(ji,jj,1) = res(ji,jj,1) + ft(n) * zamp * COS(zph) * COS(omega(n)*(rg_time(t)-t0)) &
                          &                      - ft(n) * zamp * SIN(zph) * SIN(omega(n)*(rg_time(t)-t0))
                  END DO
               END DO
               write(6,'(a,1x,4(f18.8,2x))') TRIM(tname(n)),ft(n),zamp,SIN(zph),rg_time(t)-t0
            END DO
            istatus = NF90_PUT_VAR(ncid, id_var,res,start = (/ 1, 1, t /),count = (/ jpiout, jpjout, 1 /))
         END DO

         istatus = NF90_CLOSE(ncid)

         END 

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
      parameter (jpgreg=15+31*(10+12*1582))
      INTEGER kmm, kid, kiyyy, julday
      INTEGER iy, im, ia
! ... Year 0 never existed ...
      IF (kiyyy.eq.0) stop 101
!
      IF (kiyyy.lt.0) kiyyy = kiyyy + 1
      IF (kmm.gt.2) THEN
         iy = kiyyy
         im = kmm + 1
      ELSE
         iy = kiyyy - 1
         im = kmm + 13
      END IF
!
      julday = int(365.25*iy) + int(30.6001*im) + kid + 1720995 
      IF (kid+31*(kmm+12*kiyyy).ge.jpgreg) THEN
         ia = int(0.01*iy)
         julday = julday + 2 - ia + int(0.25*ia) 
      END IF
      RETURN
      END 
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
        INTEGER jpgreg
        PARAMETER (jpgreg = 2299161)
        INTEGER kjulian, kmm, kid, kiyyy
        INTEGER ia,ialpha,ib,ic,id,ie
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
        END
