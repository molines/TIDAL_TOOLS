MODULE tide
   !!=================================================================================
   !!                       ***  MODULE  tide  ***
   !! Compute nodal modulations corrections and pulsations
   !!=================================================================================
   !!---------------------------------------------------------------------------------
   !!   OPA 9.0 , LODYC-IPSL  (2003)
   !!---------------------------------------------------------------------------------

   IMPLICIT NONE
   PUBLIC

   INTEGER, PARAMETER ::   &
       jpmax_harmo = 20,   &      ! maximum number of harmonic
        wp=8

   REAL(wp) ::             &  !:
      rpi = 3.141592653589793_wp           ,  &  !: pi
      rad = 3.141592653589793_wp / 180._wp       !: conversion from degre into radian

   REAL(wp) :: sh_T, sh_s,sh_h, sh_p, sh_p1, &
               sh_xi, sh_nu, sh_nuprim, sh_nusec, sh_R, &
               sh_I, sh_x1ra

   TYPE tide0
     CHARACTER(LEN=4)  :: name_tide
     REAL(wp) :: equitide
     INTEGER  :: nutide
     INTEGER  :: nt, ns, nh, np, np1, shift
     INTEGER  :: nksi, nnu0, nnu1, nnu2, R 
     INTEGER  :: formula 
   END TYPE tide0

   TYPE(tide0), DIMENSION(jpmax_harmo) :: Wave


CONTAINS

   SUBROUTINE tide_pulse( omega, tname ,nc)

      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE tide_pulse  ***
      !!                      
      !! ** Purpose : Compute tidal frequencies
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      CHARACTER (LEN=4), DIMENSION(nc), INTENT( in ) ::   &
         tname      ! Names of tidal constituents

      INTEGER, INTENT( in ) :: &
         nc         ! Total number of tidal constituents

      REAL (wp), DIMENSION(nc), INTENT( out ) ::   &
         omega      ! pulsation in radians/s

      !! * Local declarations
      INTEGER :: jh, jd
      INTEGER, DIMENSION(nc) :: indice
      REAL(wp) :: rl_scale  =  36525*24.0
      REAL(wp) :: omega_T=  13149000.0
      REAL(wp) :: omega_s=    481267.892
      REAL(wp) :: omega_h=     36000.76892
      REAL(wp) :: omega_p=      4069.0322056
      REAL(wp) :: omega_n=      1934.1423972
      REAL(wp) :: omega_p1=        1.719175
      !!----------------------------------------------------------------------
#  include "tide.h90"

      ! Find constituent parameters:
      indice(:) = 0
      DO jh=1,nc
        DO jd=1,jpmax_harmo
          IF (TRIM(tname(jh)) .eq. Wave(jd)%name_tide) THEN
            indice(jh) = jd
            EXIT
          END IF
        END DO
        IF (indice(jh).EQ.0) THEN
          PRINT *, 'ERROR in tide_pulse routine'
          PRINT *, 'Can not find tidal constituent:', tname(jh)
          PRINT *, 'tidal cpts list: M2 N2 2N2 S2 K2 K1 O1 Q1 P1 M4 Mf Mm Msqm Mtm'
          STOP
        END IF
      END DO

      DO jh=1,nc
        omega(jh) = omega_T * Wave(indice(jh))%nT &
                  + omega_s * Wave(indice(jh))%ns &
                  + omega_h * Wave(indice(jh))%nh &
                  + omega_p * Wave(indice(jh))%np &
                  + omega_p1* Wave(indice(jh))%np1
        omega(jh) = (omega(jh)/rl_scale)*rad/3600.
      END DO

   END SUBROUTINE tide_pulse

   SUBROUTINE tide_vuf( vt, ut , ft, tname ,nc , yy, mm, dd, HFRAC)
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE tide_vuf  ***
      !!                      
      !! ** Purpose : Compute nodal modulation corrections
      !!
      !! ** Outputs :
      !!          vt: Phase of tidal potential relative to Greenwich (radians)
      !!          ut: Phase correction u due to nodal motion (radians)
      !!          ft: Nodal correction factor
      !!
      !! ** Inputs :
      !!          tname: array of constituents names (dimension<=nc) 
      !!             nc: number of constituents
      !!             yy: current year
      !!             mm:    "    month
      !!             dd:    "    day
      !!             hh:    "    hour
      !!             min:   "    minutes
      !!             sec:   "    seconds
      !!   
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT(in) :: &
         yy, mm, dd ! Current date
      REAL(wp), INTENT(in) :: HFRAC ! Day fraction
      CHARACTER (LEN=4), DIMENSION(nc), INTENT( in ) ::   &
         tname      ! Names of tidal constituents
      INTEGER, INTENT( in ) :: &
         nc         ! Total number of tidal constituents
      REAL (wp), DIMENSION(nc), INTENT( out ) ::   &
         vt, &      !
         ut, &      !
         ft         !
      !! * Local declarations
      INTEGER :: jh, jd
      INTEGER :: X
      INTEGER, DIMENSION(nc) :: indice
      REAL(wp) :: DYR, DAYJ, DDAY
      REAL(wp) :: cosI,p,q,t2,t4,sin2I,s2,tgI2,P1,sh_tgn2,at1,at2
      REAL(wp) :: sh_N

      !!----------------------------------------------------------------------
#  include "tide.h90"
! astronomic angle

      X=AINT((yy-1901.)/4.)
      DYR=yy-1900.

      DAYJ=DAYJUL(yy,mm,dd)
      DDAY=DAYJ+X-1.

!      HFRAC=hh+min/60.+sec/3600.

!  Sh_n Longitude of ascending lunar node
!----------------------------------------------------------------------
      sh_N=(259.1560564-19.328185764*DYR-.0529539336*DDAY-.0022064139*HFRAC)*rad
! T mean solar angle (Greenwhich time)
!----------------------------------------------------------------------
      sh_T=(180.+HFRAC*(360./24.))*rad
! h mean solar Longitude
!----------------------------------------------------------------------
      sh_h=(280.1895014-.238724988*DYR+.9856473288*DDAY+.0410686387*HFRAC)*rad
! s mean lunar Longitude
!----------------------------------------------------------------------
      sh_s=(277.0256206+129.38482032*DYR+13.176396768*DDAY+.549016532*HFRAC)*rad
! p1 Longitude of solar perigee
!----------------------------------------------------------------------
      sh_p1=(281.2208569+.01717836*DYR+.000047064*DDAY+.000001961*HFRAC)*rad
! p Longitude of lunar perigee
!----------------------------------------------------------------------
      sh_p=(334.3837214+40.66246584*DYR+.111404016*DDAY+.004641834*HFRAC)*rad

      sh_N =mod(sh_N ,2*rpi)
      sh_s =mod(sh_s ,2*rpi)
      sh_h =mod(sh_h, 2*rpi)
      sh_p =mod(sh_p, 2*rpi)
      sh_p1=mod(sh_p1,2*rpi)

      cosI=0.913694997 -0.035692561 *cos(sh_N)

      sh_I=acos(cosI)

      sin2I=sin(sh_I)
      sh_tgn2=tan(sh_N/2.0)
  
      at1=atan(1.01883*sh_tgn2)
      at2=atan(0.64412*sh_tgn2)
  
      sh_xi=-at1-at2+sh_N

      if (sh_N > rpi) sh_xi=sh_xi-2.0*rpi

      sh_nu=at1-at2

! For constituents l2 k1 k2
!----------------------------------------------------------------------
      tgI2=tan(sh_I/2.0)
      P1=sh_p-sh_xi
  
      t2=tgI2*tgI2
      t4=t2*t2
      sh_x1ra=sqrt(1.0-12.0*t2*cos(2.0*P1)+36.0*t4)
  
      p=sin(2.0*P1)
      q=1.0/(6.0*t2)-cos(2.0*P1)
      sh_R=atan(p/q)
  
      p=sin(2.0*sh_I)*sin(sh_nu)
      q=sin(2.0*sh_I)*cos(sh_nu)+0.3347
      sh_nuprim=atan(p/q)
  
      s2=sin(sh_I)*sin(sh_I)
      p=s2*sin(2.0*sh_nu)
      q=s2*cos(2.0*sh_nu)+0.0727
      sh_nusec=0.5*atan(p/q)

      ! Find constituent parameters:
      indice(:) = 0
      DO jh=1,nc
        DO jd=1,jpmax_harmo
          IF (TRIM(tname(jh)) .eq. Wave(jd)%name_tide) THEN
            indice(jh) = jd
            EXIT
          END IF
        END DO
        IF (indice(jh).EQ.0) THEN
          PRINT *, 'ERROR in tide_pulse routine'
          PRINT *, 'Can not find tidal constituent:', tname(jh)
          PRINT *, 'tidal cpts list: M2 N2 2N2 S2 K2 K1 O1 Q1 P1 M4 Mf Mm Msqm Mtm'
          STOP
        END IF
      END DO

       DO jh =1,nc
!  Phase of the tidal potential relative to the Greenwhich 
!  meridian (e.g. the position of the fictuous celestial body). Units are
!  radian:
         vt(jh) = sh_T *Wave(indice(jh))%nT        &
                 +sh_s *Wave(indice(jh))%ns        &
                 +sh_h *Wave(indice(jh))%nh        &
                 +sh_p *Wave(indice(jh))%np        &
                 +sh_p1*Wave(indice(jh))%np1       &
                 +Wave(indice(jh))%shift*rad
!
!  Phase correction u due to nodal motion. Units are radian:
         ut(jh) = sh_xi    *Wave(indice(jh))%nksi  &
                 +sh_nu    *Wave(indice(jh))%nnu0  &
                 +sh_nuprim*Wave(indice(jh))%nnu1  &
                 +sh_nusec *Wave(indice(jh))%nnu2  &
                 +sh_R     *Wave(indice(jh))%R

!  Nodal correction factor:
         ft(jh) = nodal_factort(Wave(indice(jh))%formula)
       END DO

   END SUBROUTINE tide_vuf
     
   recursive function nodal_factort(formula) result (f)
   !!----------------------------------------------------------------------

   INTEGER, INTENT(IN) :: formula
   REAL(wp) :: f
   REAL(wp) :: s,f1,f2

   SELECT CASE (formula)

!!  formule 0, solar waves
  
    case ( 0 )
      f=1.0
   
!! formule 1, compound waves (78 x 78)
  
    case ( 1 )
      f=nodal_factort(78)
      f=f*f

!! formule 2, compound waves (78 x 0)  ===  (78) 
  
    case ( 2 )
      f1=nodal_factort(78)
      f=nodal_factort(0)
      f=f1*f
  
!! formule 4,  compound waves (78 x 235) 
  
    case ( 4 )
      f1=nodal_factort(78)
      f=nodal_factort(235)
      f=f1*f
  
!! formule 5,  compound waves (78 *78 x 235)
  
    case ( 5 )
      f1=nodal_factort(78)
      f=nodal_factort(235)
      f=f*f1*f1
  
!! formule 6,  compound waves (78 *78 x 0)
  
    case ( 6 )
      f1=nodal_factort(78)
      f=nodal_factort(0)
      f=f*f1*f1 
  
!! formule 7, compound waves (75 x 75)
  
    case ( 7 )
      f=nodal_factort(75)
      f=f*f
      
!! formule 8,  compound waves (78 x 0 x 235)
  
    case ( 8 )
      f=nodal_factort(78)
      f1=nodal_factort(0)
      f2=nodal_factort(235)
      f=f*f1*f2
  
!! formule 9,  compound waves (78 x 0 x 227)
  
    case ( 9 )
      f=nodal_factort(78)
      f1=nodal_factort(0)
      f2=nodal_factort(227)
      f=f*f1*f2
  
!! formule 10,  compound waves (78 x 227)
  
    case ( 10 )
      f=nodal_factort(78)
      f1=nodal_factort(227)
      f=f*f1
  
!! formule 11,  compound waves (75 x 0)
  
    case ( 11 )
      f=nodal_factort(75)
      f1=nodal_factort(0)
      f=f*f1
      
!! formule 12,  compound waves (78 x 78 x 78 x 0) 
  
    case ( 12 )
      f1=nodal_factort(78)
      f=nodal_factort(0)
      f=f*f1*f1*f1
      
!! formule 13, compound waves (78 x 75)
  
    case ( 13 )
      f1=nodal_factort(78)
      f=nodal_factort(75)
      f=f*f1
  
!! formule 14, compound waves (235 x 0)  ===  (235)
  
    case ( 14 )
      f=nodal_factort(235)
      f1=nodal_factort(0)
      f=f*f1
  
!! formule 15, compound waves (235 x 75) 
  
    case ( 15 )
      f=nodal_factort(235)
      f1=nodal_factort(75)
      f=f*f1
  
!! formule 16, compound waves (78 x 0 x 0)  ===  (78)
  
    case ( 16 )
      f=nodal_factort(78)
      f1=nodal_factort(0)
      f=f*f1*f1
      
!! formule 17,  compound waves (227 x 0) 
  
    case ( 17 )
      f1=nodal_factort(227)
      f=nodal_factort(0)
      f=f*f1
      
!! formule 18,  compound waves (78 x 78 x 78 )
  
    case ( 18 ) 
      f1=nodal_factort(78)
      f=f1*f1*f1
  
!! formule 19, compound waves (78 x 0 x 0 x 0)  ===  (78)
  
    case ( 19 )
      f=nodal_factort(78)
      f1=nodal_factort(0)
      f=f*f1*f1*f1

!! formule 73
  
    case ( 73 )
      s=sin(sh_I)
      f=(2./3.-s*s)/0.5021

!! formule 74
  
    case ( 74 )
      s=sin(sh_I)
      f=s*s/0.1578
  
!! formule 75
  
    case ( 75 )
      s=cos (sh_I/2)
      f=sin (sh_I)*s*s/0.3800

!! formule 76
  
    case ( 76 )
      f=sin (2*sh_I)/0.7214
  
!! formule 77
  
    case ( 77 )
      s=sin (sh_I/2)
      f=sin (sh_I)*s*s/0.0164
  
!! formule 78
  
    case ( 78 )
      s=cos (sh_I/2)
      f=s*s*s*s/0.9154

!! formule 79
    
    case ( 79 )
      s=sin(sh_I)
      f=s*s/0.1565
  
!! formule 144
  
    case ( 144 )
      s=sin (sh_I/2)
      f=(1-10*s*s+15*s*s*s*s)*cos(sh_I/2)/0.5873;

!! formule 149
  
    case ( 149 )
      s=cos (sh_I/2);
      f=s*s*s*s*s*s/0.8758

!! formule 215
  
    case ( 215 )
      s=cos (sh_I/2)
      f=s*s*s*s/0.9154*sh_x1ra
  
!! formule 227 
  
    case ( 227 )
      s=sin (2*sh_I)
      f=sqrt (0.8965*s*s+0.6001*s*cos (sh_nu)+0.1006)

!! formule 235 
   
    case ( 235 )
      s=sin (sh_I)
      f=sqrt (19.0444*s*s*s*s+2.7702*s*s*cos (2*sh_nu)+.0981)

    END SELECT
end function nodal_factort

      FUNCTION DAYJUL(YR,MONTH,DAY)
!
!*** THIS ROUTINE COMPUTES THE JULIAN DAY (AS A REAL VARIABLE)
!
      INTEGER :: YR,MONTH,DAY
      INTEGER ::  DAYT(12),DAYS(12)
      INTEGER :: DINC,I
      REAL(wp) :: DAYJUL,YRLP     

      DATA DAYT/0.,31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334./
      DAYS(1)=0.
      DAYS(2)=31.
      DINC=0.
      YRLP=MOD((YR-1900.),4.)
      IF(YRLP .EQ. 0.) DINC=1.
      DO I=3,12
       DAYS(I)=DAYT(I)+DINC
      END DO
      DAYJUL=DAYS(MONTH)+DAY

      END FUNCTION DAYJUL

END MODULE tide
