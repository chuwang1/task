! trloop.f90

MODULE trloop

  PRIVATE
  PUBLIC tr_loop
  PUBLIC tr_scaling_adjust

CONTAINS

!     ***********************************************************

!           MAIN ROUTINE FOR TRANSPORT CALCULATION

!     ***********************************************************

  SUBROUTINE tr_loop(ierr)

      USE TRCOMM
      USE trbpsd, ONLY: tr_bpsd_put, tr_bpsd_get,plasmaf
      USE trexec
      USE trprof, ONLY: tr_reset_density
      USE libitp
      USE equnit
      IMPLICIT NONE
      INTEGER,INTENT(OUT):: IERR
      INTEGER:: nr

      ierr=0
      IF(NT.GE.NTMAX) GOTO 9000

      CALL tr_eval(NT,IERR)
      IF(IERR.NE.0) GOTO 9000

      RIP=RIPS
      IF(NTMAX.NE.0) DIPDT=(RIPE-RIPS)/(DBLE(NTMAX)*DT)
      write(6,'(A,1P4E12.4)') "**RIP,RIPS,RIPE,DIP=",RIP,RIPS,RIPE,DIPDT

      call tr_bpsd_get(ierr)
      if(ierr.ne.0) GOTO 9000

 1000 CONTINUE

      CALL tr_exec(IERR)
      IF(IERR.NE.0) GOTO 9000

      ! Reset density from profile if model_nevolve=1
      IF(model_nevolve.EQ.1) CALL tr_reset_density

      DO nr=1,nrmax
         QPINV(nr)=(4.D0*PI**2*RDPVRHOG(nr))/(TTRHOG(nr)*ARRHOG(nr))
      END DO

!     /* Sawtooth Oscillation */
      Q0=FCTR(RG(1),RG(2),QP(1),QP(2))
      IF(Q0.LT.1.D0) TST=TST+DT

      IF(TST+0.5D0*DT.GT.TPRST) THEN
         CALL TRSAWT
         TST=0.D0
      ENDIF

      call tr_bpsd_put(IERR)
      
      if(ierr.ne.0) GOTO 9000
      CALL tr_scaling_adjust(NT)
      NT=NT+1

!     *** SET GEOMETRY VIA TASK/EQ ***

      IF(MODELG.EQ.5 .AND. EQRELOAD.EQ.1) THEN
         CALL eq_load(MODELG,KNAMEQ,IERR)
         IF(IERR.NE.0) THEN
            WRITE(6,*) 'XX eq_load(reload): ierr=',ierr,' knameq=',TRIM(KNAMEQ)
            RETURN
         ENDIF
         EQRELOAD=0
      ENDIF

      IF(NTEQIT.NE.0) THEN
         IF(MOD(NT,NTEQIT).EQ.0) THEN
            if(modelg.eq.8) THEN
!               call equ_calc
            endif
            IF(modelg.eq.9) THEN
               call eq_calc
            endif
         ENDIF
      ENDIF
      call tr_bpsd_get(IERR)
      if(ierr.ne.0) return

      CALL tr_eval(NT,IERR)
      IF(IERR.NE.0) GOTO 9000

      IF(NT.LT.NTMAX) GOTO 1000

9000  CONTINUE
      RIPS=RIPE
      RETURN
    END SUBROUTINE tr_loop

  SUBROUTINE tr_scaling_adjust(NTSTEP_IN)

      USE TRCOMM
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NTSTEP_IN
      REAL(rkind) :: ERR_TAUE, ERR_TAUE_RAW, C_SCALING_NEW

      IF(MDLKAI < 180 .OR. MDLKAI > 189) RETURN

      SELECT CASE(MDLKAI)
      CASE(180)
         TAUE_TARGET = TAUE89
      CASE(181)
         TAUE_TARGET = TAUE98
      CASE DEFAULT
         TAUE_TARGET = TAUE98
      END SELECT

      IF(TAUE_TARGET < 1.D-10) TAUE_TARGET = 1.D-2
      ERR_TAUE_RAW = (TAUE2 - TAUE_TARGET) / TAUE_TARGET
      ERR_TAUE = MIN(MAX(ERR_TAUE_RAW, -0.2D0), 0.2D0)

      C_SCALING_NEW = C_SCALING * (1.D0 + ALPHA_RELAX * ERR_TAUE)
      IF(C_SCALING_NEW > C_SCALING * 1.02D0) C_SCALING_NEW = C_SCALING * 1.02D0
      IF(C_SCALING_NEW < C_SCALING * 0.98D0) C_SCALING_NEW = C_SCALING * 0.98D0
      C_SCALING_NEW = MIN(MAX(C_SCALING_NEW, C_SCALING_MIN), C_SCALING_MAX)

      L_SCALING_CONVERGED = ABS(ERR_TAUE) < 0.01D0
      C_SCALING = C_SCALING_NEW

      IF(MOD(NTSTEP_IN,100).EQ.0) THEN
         WRITE(6,'(A,I6,A,F8.4,A,F8.4,A,F8.4,A,F7.2,A)') &
            'SCALING: NT=', NTSTEP_IN, &
            ' tau_sim=', TAUE2, &
            ' tau_tgt=', TAUE_TARGET, &
            ' C=', C_SCALING, &
            ' err=', ERR_TAUE_RAW*100.D0, '%'
      ENDIF
  END SUBROUTINE tr_scaling_adjust
  END MODULE trloop
