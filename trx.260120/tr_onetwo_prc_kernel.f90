MODULE tr_onetwo_prc_kernel

!------------------------------------------------------------------------------
! Pure ONETWO-style cyclotron/synchrotron radiation kernel.
!
! This module deliberately has no dependency on trcomm.  It accepts the same
! physical inputs used by PyMak's calculate_cyclotron_radiation_onetwo routine
! and returns PRC in TR internal units, W/m**3.
!------------------------------------------------------------------------------

  USE bpsd_kinds, ONLY: rkind
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: onetwo_prc_point
  PUBLIC :: onetwo_prc_profile

  REAL(KIND=rkind), PARAMETER :: PI_CGS = &
       3.141592653589793238462643383279502884197_rkind
  REAL(KIND=rkind), PARAMETER :: E_CHARGE_CGS = &
       4.803204712570263E-10_rkind
  REAL(KIND=rkind), PARAMETER :: M_E_CGS = &
       9.1093837015E-28_rkind
  REAL(KIND=rkind), PARAMETER :: C_CGS = &
       2.99792458E10_rkind
  REAL(KIND=rkind), PARAMETER :: KEV_TO_ERG = &
       1.602176634E-9_rkind
  REAL(KIND=rkind), PARAMETER :: KEV_TO_J = &
       1.602176634E-16_rkind
  REAL(KIND=rkind), PARAMETER :: T_TO_GAUSS = 1.0E4_rkind
  REAL(KIND=rkind), PARAMETER :: M_TO_CM = 1.0E2_rkind
  REAL(KIND=rkind), PARAMETER :: NE_1E20_M3_TO_CM3 = 1.0E14_rkind
  REAL(KIND=rkind), PARAMETER :: CM3_TO_M3 = 1.0E6_rkind
  REAL(KIND=rkind), PARAMETER :: TINY_DENOM = 1.0E-300_rkind

CONTAINS

  PURE SUBROUTINE onetwo_prc_point(ne_1e20_m3, te_kev, bt_t, a_m, r0_m, &
       wall_reflection, prc_w_m3, phi_bar)

    REAL(KIND=rkind), INTENT(IN) :: ne_1e20_m3
    REAL(KIND=rkind), INTENT(IN) :: te_kev
    REAL(KIND=rkind), INTENT(IN) :: bt_t
    REAL(KIND=rkind), INTENT(IN) :: a_m
    REAL(KIND=rkind), INTENT(IN) :: r0_m
    REAL(KIND=rkind), INTENT(IN) :: wall_reflection
    REAL(KIND=rkind), INTENT(OUT) :: prc_w_m3
    REAL(KIND=rkind), INTENT(OUT), OPTIONAL :: phi_bar

    REAL(KIND=rkind) :: wall_reflection_clamped
    REAL(KIND=rkind) :: bt_gauss
    REAL(KIND=rkind) :: a_cm
    REAL(KIND=rkind) :: r0_cm
    REAL(KIND=rkind) :: ne_cm3
    REAL(KIND=rkind) :: te_erg
    REAL(KIND=rkind) :: mec2_erg
    REAL(KIND=rkind) :: omega_ce
    REAL(KIND=rkind) :: omega_pe_sq
    REAL(KIND=rkind) :: rm
    REAL(KIND=rkind) :: chi
    REAL(KIND=rkind) :: phi_local
    REAL(KIND=rkind) :: phi_sqrt_arg
    REAL(KIND=rkind) :: vdotsq
    REAL(KIND=rkind) :: qsync_kev_cm3_s
    REAL(KIND=rkind) :: qsync_kev_m3_s

    prc_w_m3 = 0.0_rkind
    IF (PRESENT(phi_bar)) phi_bar = 0.0_rkind

    IF (ne_1e20_m3 <= 0.0_rkind) RETURN
    IF (te_kev <= 0.0_rkind) RETURN
    IF (ABS(bt_t) <= 0.0_rkind) RETURN
    IF (a_m <= 0.0_rkind) RETURN
    IF (r0_m <= 0.0_rkind) RETURN

    wall_reflection_clamped = MIN(MAX(wall_reflection, 0.0_rkind), 1.0_rkind)

    bt_gauss = ABS(bt_t) * T_TO_GAUSS
    a_cm = a_m * M_TO_CM
    r0_cm = r0_m * M_TO_CM
    ne_cm3 = ne_1e20_m3 * NE_1E20_M3_TO_CM3
    te_erg = te_kev * KEV_TO_ERG
    mec2_erg = M_E_CGS * C_CGS * C_CGS

    omega_ce = E_CHARGE_CGS * bt_gauss / (M_E_CGS * C_CGS)
    omega_pe_sq = 4.0_rkind * PI_CGS * ne_cm3 * E_CHARGE_CGS * &
         E_CHARGE_CGS / M_E_CGS

    rm = a_cm / MAX(r0_cm, TINY_DENOM)
    chi = rm * SQRT(MAX(mec2_erg / MAX(te_erg, TINY_DENOM), 0.0_rkind))

    phi_sqrt_arg = C_CGS * omega_ce / MAX(a_cm * omega_pe_sq, TINY_DENOM)
    phi_sqrt_arg = phi_sqrt_arg * (1.0_rkind - wall_reflection_clamped) * &
         (1.0_rkind + chi)

    phi_local = 60.0_rkind * (MAX(te_erg, 0.0_rkind) / mec2_erg)**1.5_rkind
    phi_local = phi_local * SQRT(MAX(phi_sqrt_arg, 0.0_rkind))
    IF (PRESENT(phi_bar)) phi_bar = phi_local

    vdotsq = omega_ce * omega_ce * 2.0_rkind * te_erg / M_E_CGS

    qsync_kev_cm3_s = ne_cm3 / 1.5_rkind * E_CHARGE_CGS * E_CHARGE_CGS / &
         (C_CGS**3)
    qsync_kev_cm3_s = qsync_kev_cm3_s * vdotsq * phi_local / KEV_TO_ERG
    qsync_kev_m3_s = qsync_kev_cm3_s * CM3_TO_M3

    prc_w_m3 = qsync_kev_m3_s * KEV_TO_J

  END SUBROUTINE onetwo_prc_point


  PURE SUBROUTINE onetwo_prc_profile(nrmax, ne_1e20_m3, te_kev, bt_t, a_m, &
       r0_m, wall_reflection, prc_w_m3, phi_bar)

    INTEGER, INTENT(IN) :: nrmax
    REAL(KIND=rkind), INTENT(IN) :: ne_1e20_m3(nrmax)
    REAL(KIND=rkind), INTENT(IN) :: te_kev(nrmax)
    REAL(KIND=rkind), INTENT(IN) :: bt_t
    REAL(KIND=rkind), INTENT(IN) :: a_m
    REAL(KIND=rkind), INTENT(IN) :: r0_m
    REAL(KIND=rkind), INTENT(IN) :: wall_reflection
    REAL(KIND=rkind), INTENT(OUT) :: prc_w_m3(nrmax)
    REAL(KIND=rkind), INTENT(OUT), OPTIONAL :: phi_bar(nrmax)

    INTEGER :: nr
    REAL(KIND=rkind) :: phi_local

    DO nr = 1, nrmax
       CALL onetwo_prc_point(ne_1e20_m3(nr), te_kev(nr), bt_t, a_m, r0_m, &
            wall_reflection, prc_w_m3(nr), phi_local)
       IF (PRESENT(phi_bar)) phi_bar(nr) = phi_local
    END DO

  END SUBROUTINE onetwo_prc_profile

END MODULE tr_onetwo_prc_kernel
