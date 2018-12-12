MODULE utils
  IMPLICIT NONE
  !!======================================================================
  !!                     ***  MODULE  utils  ***
  !! Provides common routines/function for the tidal tools
  !!=====================================================================
  !! History : 1.0  !  12/2018    J.M.Molines  set up the module
  !!----------------------------------------------------------------------

  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  !! julday          : convert calendar date to julian day
  !! caldat          : convert julian day to calendar date 
  !! chkfile         : check existance of a file
  PUBLIC julday
  PUBLIC caldat
  PUBLIC chkfile

  !!----------------------------------------------------------------------
  !! TIDAL_TOOLS , MEOM 2018
  !! $Id$
  !! Copyright (c) 2018, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/TIDAL_TOOLSCeCILL.txt)
  !!----------------------------------------------------------------------

CONTAINS
  FUNCTION julday(kmm,kid,kiyyy)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION julday  ***
    !!
    !! ** Purpose :  This function returns the julian day number which 
    !!               begins at noon  of the calendar date specified by 
    !!               month kmm, day kid, and year kiyyy.  
    !!               Positive year signifies a.d.; negative, b.c.  
    !!               (remember that the year after 1 b.c. was 1 a.d.)
    !!               This function  handles changeover to gregorian 
    !!               calendar on oct. 15, 1582.
    !!
    !! ** Method  :  This routine comes directly from the Numerical Recipe Book. 
    !!
    !! ** References : numerical recipes, cambridge univ. press, 1986.
    !!----------------------------------------------------------------------
    INTEGER, INTENT(IN) :: kmm, kid, kiyyy  ! month, day, year for the input date
    ! local varaiables
    INTEGER, PARAMETER :: jpgreg = 15+31*(10+12*1582)
    INTEGER :: il_kiyyy, julday
    INTEGER iy, im, ia
    !!----------------------------------------------------------------------
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
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE caldat  ***
    !!
    !! ** Purpose :  This routine convert a julian day in calendar date
    !!               It is the inverse of the function julday 
    !!
    !! ** Method  :  This routine comes directly from the Numerical Recipes Book
    !!
    !! References : Numerical recipes, Cambridge Univ. Press, 1986
    !!----------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: kjulian             ! input julian day number
    INTEGER, INTENT(OUT) :: kmm, kid, kiyyy     ! output for month, day, year
                                                !  year is > 0 a.d or < 0 if b.c.
    ! local variables
    INTEGER, PARAMETER   :: jpgreg = 2299161
    INTEGER :: ia, ialpha, ib, ic, id, ie
    !------------------------------------------------------------------------
    ! ... Cross over to Gregorian Calendar produces this correction:
    !
    IF ( kjulian >=  jpgreg) THEN
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
    IF ( kmm >  12 ) kmm = kmm - 12
    kiyyy = ic - 4715
    IF ( kmm   >  2 ) kiyyy = kiyyy - 1
    IF ( kiyyy <= 0 ) kiyyy = kiyyy - 1
    RETURN
  END SUBROUTINE caldat

  LOGICAL FUNCTION chkfile (cd_file, ld_verbose )
    !!----------------------------------------------------------------------
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
    LOGICAL, OPTIONAL, INTENT(in) :: ld_verbose   ! if false, do not emit warning

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
END MODULE utils
