MODULE surdetermine 

  IMPLICIT NONE
  PUBLIC

  INTEGER, PARAMETER :: NBINCOMAX = 50 
  INTEGER, PARAMETER :: DIMSPARSE = NBINCOMAX*4000*24
  INTEGER, PARAMETER :: wp=8

  INTEGER :: NBSPARSE, NBINCO
  REAL(wp), DIMENSION(DIMSPARSE) :: SPARSEVALUE
  INTEGER, DIMENSION(DIMSPARSE) :: JSPARSE , ISPARSE

  INTEGER, SAVE, DIMENSION(NBINCOMAX) :: JPOS1
  REAL(wp), DIMENSION(NBINCOMAX) :: TAB4, TAB7
  REAL(wp), SAVE, DIMENSION(NBINCOMAX,NBINCOMAX) :: TAB3, PILIER
  REAL(wp), SAVE, DIMENSION(NBINCOMAX) :: PIVOT

  !!---------------------------------------------------------------------------------
  !!
  !!---------------------------------------------------------------------------------

CONTAINS

  SUBROUTINE SUR_DETERMINE(init)

    INTEGER, INTENT(in) :: init

    INTEGER  :: &
         I_SD, I1_SD, II_SD, J_SD, K1_SD, K2_SD
    REAL(wp) :: VALEUR1, VALEUR2, X1
    REAL(wp), DIMENSION(NBINCOMAX) :: TABX, COL1, COL2
    INTEGER, DIMENSION(NBINCOMAX) :: JPOS2, JPIVOT      
    !---------------------------------------------------------------------------------

    IF (init==1) THEN
       IF(NBSPARSE.GT.DIMSPARSE) STOP 'surdetermine erreur1'
       IF(NBINCO.GT.NBINCOMAX)THEN
          PRINT*,'NBINCO   =',NBINCO
          PRINT*,'NBINCOMAX=',NBINCOMAX
          STOP 'DONC dans surdetermine erreur2'
       END IF

       TAB3(:,:)=0.e0

       DO K1_SD=1,NBSPARSE
          DO K2_SD=1,NBSPARSE 
             ISPARSE(K2_SD)=ISPARSE(K2_SD)
             JSPARSE(K2_SD)=JSPARSE(K2_SD)
             IF(ISPARSE(K2_SD).EQ.ISPARSE(K1_SD)) &
                  TAB3(JSPARSE(K1_SD),JSPARSE(K2_SD))= &
                  TAB3(JSPARSE(K1_SD),JSPARSE(K2_SD))  &
                  +SPARSEVALUE(K1_SD)*SPARSEVALUE(K2_SD)
          END DO
       END DO

       DO J_SD=1,NBINCO
          JPOS1(J_SD)=J_SD
          JPOS2(J_SD)=J_SD
       ENDDO


       DO I_SD=1,NBINCO
          ! recherche du plus grand pivot non nul
          VALEUR1=ABS(TAB3(I_SD,I_SD))

          JPIVOT(I_SD)=I_SD
          DO J_SD=I_SD,NBINCO
             VALEUR2=ABS(TAB3(I_SD,J_SD))
             IF(VALEUR2.GE.VALEUR1) THEN
                JPIVOT(I_SD)=J_SD
                VALEUR1=VALEUR2
             ENDIF
          END DO

          DO I1_SD=1,NBINCO
             COL1(I1_SD)=TAB3(I1_SD,I_SD)
             COL2(I1_SD)=TAB3(I1_SD,JPIVOT(I_SD))
             TAB3(I1_SD,I_SD)=COL2(I1_SD)
             TAB3(I1_SD,JPIVOT(I_SD))=COL1(I1_SD)
          END DO

          JPOS2(I_SD)=JPOS1(JPIVOT(I_SD))
          JPOS2(JPIVOT(I_SD))=JPOS1(I_SD)
          JPOS1(I_SD)=JPOS2(I_SD)
          JPOS1(JPIVOT(I_SD))=JPOS2(JPIVOT(I_SD))


          !-------------------------------
          PIVOT(I_SD)=TAB3(I_SD,I_SD)
          DO J_SD=1,NBINCO
             TAB3(I_SD,J_SD)=TAB3(I_SD,J_SD)/PIVOT(I_SD)
          END DO
          !-------------------------------

          !-------------------------------
          DO II_SD=I_SD+1,NBINCO
             PILIER(II_SD,I_SD)=TAB3(II_SD,I_SD)
             DO J_SD=1,NBINCO
                TAB3(II_SD,J_SD)= &
                     TAB3(II_SD,J_SD)-TAB3(I_SD,J_SD)*PILIER(II_SD,I_SD)
             END DO
          END DO

       END DO
    ELSE ! End init==1

       !
       DO I_SD=1,NBINCO
          TAB4(I_SD)=TAB4(I_SD)/PIVOT(I_SD)
          DO II_SD=I_SD+1,NBINCO
             TAB4(II_SD)=TAB4(II_SD)-TAB4(I_SD)*PILIER(II_SD,I_SD)
          END DO
       END DO

       !  resolution du systeme:
       TABX(NBINCO)=TAB4(NBINCO)/TAB3(NBINCO,NBINCO)
       I_SD=NBINCO
       DO I_SD=NBINCO-1,1,-1
          X1=0.
          DO J_SD=I_SD+1,NBINCO
             X1=X1+TABX(J_SD)*TAB3(I_SD,J_SD)
          END DO
          TABX(I_SD)=TAB4(I_SD)-X1
       END DO

       DO J_SD=1,NBINCO
          TAB7(JPOS1(J_SD))=TABX(J_SD)
       END DO
    ENDIF

  END SUBROUTINE SUR_DETERMINE

END MODULE surdetermine
