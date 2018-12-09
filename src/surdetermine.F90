MODULE surdetermine 
  !!======================================================================
  !!                     ***  MODULE  surdetermin  ***
  !! Provide  routines for inversion of sparse matrix
  !!=====================================================================
  !! History : 1.  12/2018 : J.M. Moline from TIDE Mercator adn ???
  !!----------------------------------------------------------------------

  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  !!  SUR_DETERMINE : 
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4), PARAMETER :: jp_bincomax = 50 
  INTEGER(KIND=4), PARAMETER :: jp_dimsparse = jp_bincomax*4000*24

  PUBLIC  SUR_DETERMINE
  INTEGER(KIND=4), PUBLIC, DIMENSION(jp_dimsparse) :: JSPARSE , ISPARSE
  INTEGER(KIND=4), PUBLIC                          :: nbsparse, nbinco

  REAL(KIND=8),    PUBLIC, DIMENSION(jp_dimsparse) :: dsparsevalue
  REAL(KIND=8),    PUBLIC, DIMENSION(jp_bincomax ) :: dtab4, dtab7

  PRIVATE
  INTEGER(KIND=4), SAVE, DIMENSION(jp_bincomax)             :: ijpos1
  REAL(KIND=8),    SAVE, DIMENSION(jp_bincomax,jp_bincomax) :: ds_tab3, ds_pilier
  REAL(KIND=8),    SAVE, DIMENSION(jp_bincomax)             :: ds_pivot

  !!---------------------------------------------------------------------------------
  !! TIDAL_TOOLS , MEOM 2018
  !! $Id$
  !! Copyright (c) 2018, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/TIDAL_TOOLSCeCILL.txt)
  !!----------------------------------------------------------------------

CONTAINS

  SUBROUTINE SUR_DETERMINE(kinit)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE SUR_DETERMINE  ***
    !!
    !! ** Purpose :  Sparse Matrix inverter 
    !!
    !! ** Method  :    ???
    !!
    !! References :   ???
    !!----------------------------------------------------------------------
    INTEGER(KIND=4), INTENT(in) :: kinit
    ! Local variables
    INTEGER(KIND=4)                         :: ji_sd, ji1_sd, jii_sd, jj_sd, j1_sd, j2_sd   ! Loop index
    INTEGER(KIND=4), DIMENSION(jp_bincomax) :: ijpos2, ipivot      

    REAL(KIND=8)                            :: dl_val1, dl_val2, dl_x1
    REAL(KIND=8),    DIMENSION(jp_bincomax) :: dl_tabx, dl_col1, dl_col2
    !---------------------------------------------------------------------------------
    IF ( kinit == 1 ) THEN
       IF( nbsparse  > jp_dimsparse) STOP 'surdetermine init error 1'
       IF(nbinco     > jp_bincomax ) THEN
          PRINT*,'nbinco      =',nbinco
          PRINT*,'jp_bincomax =',jp_bincomax
          STOP '==> surdetermine init erreur2'
       END IF

       ds_tab3(:,:)=0.e0

       DO j1_sd=1,nbsparse
          DO j2_sd=1,nbsparse 
!! JMM 2 next line stupid ??
             ISPARSE(j2_sd)=ISPARSE(j2_sd)
             JSPARSE(j2_sd)=JSPARSE(j2_sd)
!!

             IF(ISPARSE(j2_sd)  == ISPARSE(j1_sd)) THEN
                  ds_tab3(JSPARSE(j1_sd),JSPARSE(j2_sd))= &
                  &     ds_tab3(JSPARSE(j1_sd),JSPARSE(j2_sd)) + dsparsevalue(j1_sd)*dsparsevalue(j2_sd)
             ENDIF
          END DO
       END DO

       DO jj_sd=1,nbinco
          ijpos1(jj_sd)=jj_sd
          ijpos2(jj_sd)=jj_sd
       ENDDO


       DO ji_sd=1,nbinco
          ! Look for the larger non null pivot
          dl_val1=ABS(ds_tab3(ji_sd,ji_sd))

          ipivot(ji_sd)=ji_sd
          DO jj_sd=ji_sd,nbinco
             dl_val2=ABS(ds_tab3(ji_sd,jj_sd))
             IF(dl_val2  >= dl_val1) THEN
                ipivot(ji_sd)=jj_sd
                dl_val1=dl_val2
             ENDIF
          END DO

          DO ji1_sd=1,nbinco
             dl_col1(ji1_sd) = ds_tab3(ji1_sd,ji_sd)
             dl_col2(ji1_sd) = ds_tab3(ji1_sd,ipivot(ji_sd))

             ds_tab3(ji1_sd,ji_sd)         = dl_col2(ji1_sd)
             ds_tab3(ji1_sd,ipivot(ji_sd)) = dl_col1(ji1_sd)
          END DO

          ijpos2(ji_sd)         = ijpos1(ipivot(ji_sd))
          ijpos2(ipivot(ji_sd)) = ijpos1(ji_sd)

          ijpos1(ji_sd)         = ijpos2(ji_sd)
          ijpos1(ipivot(ji_sd)) = ijpos2(ipivot(ji_sd))

          ds_pivot(ji_sd)       = ds_tab3(ji_sd,ji_sd)

          DO jj_sd=1,nbinco
             ds_tab3(ji_sd,jj_sd) = ds_tab3(ji_sd,jj_sd)/ds_pivot(ji_sd)
          END DO

          DO jii_sd=ji_sd+1,nbinco
             ds_pilier(jii_sd,ji_sd) = ds_tab3(jii_sd,ji_sd)
             DO jj_sd=1,nbinco
                ds_tab3(jii_sd,jj_sd)= &
                  &  ds_tab3(jii_sd,jj_sd)-ds_tab3(ji_sd,jj_sd)*ds_pilier(jii_sd,ji_sd)
             END DO
          END DO

       END DO
    ELSE ! End kinit==1

       !
       DO ji_sd=1,nbinco
          dtab4(ji_sd)  =  dtab4(ji_sd)/ds_pivot(ji_sd)
          DO jii_sd=ji_sd+1,nbinco
             dtab4(jii_sd)  = dtab4(jii_sd)-dtab4(ji_sd)*ds_pilier(jii_sd,ji_sd)
          END DO
       END DO

       !  Solve the system
       dl_tabx(nbinco) = dtab4(nbinco)/ds_tab3(nbinco,nbinco)
       ji_sd  = nbinco
       DO ji_sd = nbinco-1,1,-1
          dl_x1=0.
          DO jj_sd=ji_sd+1,nbinco
             dl_x1 = dl_x1 + dl_tabx(jj_sd)*ds_tab3(ji_sd,jj_sd)
          END DO
          dl_tabx(ji_sd) = dtab4(ji_sd) - dl_x1
       END DO

       DO jj_sd=1,nbinco
          dtab7(ijpos1(jj_sd)) = dl_tabx(jj_sd)
       END DO
    ENDIF

  END SUBROUTINE sur_determine

END MODULE surdetermine
