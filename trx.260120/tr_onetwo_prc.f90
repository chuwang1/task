MODULE tr_onetwo_prc_mod

!-----------------------------------------
! Interface of TASK/TR with ONETWO PRC
!-----------------------------------------

CONTAINS

SUBROUTINE tr_onetwo_prc

USE bpsd_kinds
USE trcomm, ONLY: nrmax,rn,rt,prc,bb,ra,rr,refrad
USE tr_onetwo_prc_kernel, ONLY: onetwo_prc_profile

IMPLICIT NONE

REAL(KIND=rkind) :: phi_bar(nrmax)

CALL onetwo_prc_profile(nrmax, rn(1:nrmax,1), rt(1:nrmax,1), &
     bb, ra, rr, refrad, prc(1:nrmax), phi_bar)

RETURN

END SUBROUTINE tr_onetwo_prc

END MODULE tr_onetwo_prc_mod
