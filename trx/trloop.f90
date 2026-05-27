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
      USE trfixed, ONLY: tr_prep_arfixed
      USE libitp
      USE equnit
      IMPLICIT NONE
      INTEGER,INTENT(OUT):: IERR
      INTEGER:: nr

      ierr=0
      IF(NT.GE.NTMAX) GOTO 9000

      ! Check if Ar density needs to be loaded (for 'c' continue command)
      IF(model_arfixed.EQ.1 .AND. MAXVAL(ANAR(1:NRMAX)).LT.1.D-20) THEN
         CALL tr_prep_arfixed
      END IF

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

      ! Auto-adjust C_SCALING for scaling-based transport model (MDLKAI=180-189)
      CALL tr_scaling_adjust(NT)

      NT=NT+1

!     *** SET GEOMETRY VIA TASK/EQ ***

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

!     ***********************************************************

!           AUTO-ADJUST C_SCALING FOR SCALING-BASED TRANSPORT

!     ***********************************************************

  SUBROUTINE tr_scaling_adjust(NTSTEP_IN)

      USE TRCOMM
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NTSTEP_IN
      REAL(rkind) :: ERR_TAUE, ERR_TAUE_RAW, C_SCALING_NEW

      ! Only apply for scaling-based transport models (MDLKAI=180-189)
      IF(MDLKAI < 180 .OR. MDLKAI > 189) RETURN

      ! Select target confinement time based on MDLKAI
      SELECT CASE(MDLKAI)
      CASE(180)
         TAUE_TARGET = TAUE89
      CASE(181)
         TAUE_TARGET = TAUE98
      CASE(182)
         TAUE_TARGET = H_FACTOR_USER * TAUE89
      CASE(183)
         TAUE_TARGET = H_FACTOR_USER * TAUE98
      CASE DEFAULT
         TAUE_TARGET = TAUE98
      END SELECT

      ! Prevent division by zero
      IF(TAUE_TARGET < 1.D-10) TAUE_TARGET = 1.D-2

      ! Calculate error: e = (tau_sim - tau_target) / tau_target
      ERR_TAUE_RAW = (TAUE2 - TAUE_TARGET) / TAUE_TARGET
      ERR_TAUE = ERR_TAUE_RAW

      ! Limit maximum adjustment per step to prevent numerical instability
      ! Max adjustment factor is 20% per step
      ERR_TAUE = MIN(MAX(ERR_TAUE, -0.2D0), 0.2D0)

      ! Adjust C using relaxation method
      ! If tau_sim > tau_target: need more transport (increase C)
      ! If tau_sim < tau_target: need less transport (decrease C)
      C_SCALING_NEW = C_SCALING * (1.D0 + ALPHA_RELAX * ERR_TAUE)

      ! Limit maximum change per step to 2% to prevent oscillation
      IF(C_SCALING_NEW > C_SCALING * 1.02D0) C_SCALING_NEW = C_SCALING * 1.02D0
      IF(C_SCALING_NEW < C_SCALING * 0.98D0) C_SCALING_NEW = C_SCALING * 0.98D0

      ! Apply limits
      C_SCALING_NEW = MIN(MAX(C_SCALING_NEW, C_SCALING_MIN), C_SCALING_MAX)

      ! Check convergence
      IF(ABS(ERR_TAUE) < 0.01D0) THEN  ! 1% tolerance
         L_SCALING_CONVERGED = .TRUE.
      ELSE
         L_SCALING_CONVERGED = .FALSE.
      ENDIF

      ! Update C_SCALING
      C_SCALING = C_SCALING_NEW

      ! Output diagnostic information every 100 steps
      IF(MOD(NTSTEP_IN, 100) == 0) THEN
         WRITE(6,'(A,I6,A,F8.4,A,F8.4,A,F8.4,A,F7.2,A)') &
            'SCALING: NT=', NTSTEP_IN, &
            ' tau_sim=', TAUE2, &
            ' tau_tgt=', TAUE_TARGET, &
            ' C=', C_SCALING, &
            ' err=', ERR_TAUE_RAW*100.D0, '%'
      ENDIF

      RETURN
  END SUBROUTINE tr_scaling_adjust

  END MODULE trloop
