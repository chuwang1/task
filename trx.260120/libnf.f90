! libnf.f90

MODULE libnf_local
  USE bpsd_kinds
  INTEGER:: id_nf_local,mode_nf_local
  REAL(rkind):: pm_local,temperature_local
END MODULE libnf_local

MODULE libnf
  USE bpsd_kinds
  USE bpsd_constants
  USE trcomm
  IMPLICIT NONE

  ! Fusion model
  !   model_pnf=0 : no fusion reaction
  !   model_pnf=1 : D + T -> He4 + n              nnfmax=1  nsmax=4 DT
  !   model_pnf=2 : D + D -> T + p                nnfmax=4  nsmax=6 DD1 DD2
  !                 D + D -> He3 + n                                DD3
  !                 D + T -> He4 + n                                DT
  !   model_pnf=3 : D + D -> T + p                nnfmax=6  nsmax=6 DD1 DD2
  !                 D + D -> He3 + n                                DD3
  !                 D + T -> He4 + n                                DT
  !                 D + He3 -> He4 + p                              DHe31 DHe32
  !   model_pnf=4 : D + D -> T + p                nnfmax=13 nsmax=7 DD1 DD2
  !                 D + D -> He3 + n                                DD3
  !                 D + T -> He4 + n                                DT
  !                 D + He3 -> He4 + p                              DH31 DHe32
  !                 T + T -> He4 + 2n                               TT
  !                 T + He3 -> He4 + p + n                          THe31 THe32
  !                 T + He3 -> He4 + D                              THe33 THe34
  !                 T + He3 -> He5 + p -> He4 + n + p               THe35 THe36
  
  ! Fusion reaction id
  
  INTEGER,PARAMETER,PUBLIC:: id_nf_DT=1   ! mode_nf=1 D + T   -> <He4> +  n
  INTEGER,PARAMETER,PUBLIC:: id_nf_DD=2   ! mode_nf=1 D + D   -> <T>   +  p
                                          ! mode_nf=2 D + D   ->  T    + <p> 
                                          ! mode_nf=3 D + D   -> <He3> +  n
  INTEGER,PARAMETER,PUBLIC:: id_nf_DHe3=3 ! mode_nf=1 D + He3 -> <He4> +  p
                                          ! mode_nf=2 D + He3 ->  He4  + <p>
  INTEGER,PARAMETER,PUBLIC:: id_nf_TT=4   ! mode_nf=1 T + T   -> <He4> + 2n
  INTEGER,PARAMETER,PUBLIC:: id_nf_THe3=5 ! mode_nf=1 T + He3 -> <He4> +  p  +n
                                          ! mode_nf=2 T + He3 ->  He4  + <p> +n
                                          ! mode_nf=3 T + He3 -> <He4> +  D
                                          ! mode_nf=4 T + He3 ->  He4  + <D>
                                          ! mode_nf=5 T + He3 -> <He5> +  p
                                          ! mode_nf=6 T + He3 ->  He5  + <p>
  INTEGER,PARAMETER,DIMENSION(5):: &
       mode_nf_max=(/1,3,2,1,6/)
  
  REAL(rkind),DIMENSION(4,10):: &
       usvnf_dd,usvnf_dt,usvnf_dhe3,usvnf_tt,usvnf_the3
  REAL(rkind),DIMENSION(10):: &
       tempa_log
  REAL(rkind),DIMENSION(10):: &
       svnf_dt_log,svnf_dd_log,svnf_dhe3_log,svnf_tt_log,svnf_the3_log
  REAL(rkind),DIMENSION(10):: &
       dsvnf

    ! *** reaction variables ***

  INTEGER,DIMENSION(13):: &
       id_nf_nnf,mode_nf_nnf
  Integer,DIMENSION(5,6):: &
       ns1_idnf,ns2_idnf,nsp_idnf
  REAL(rkind),DIMENSION(5,6)::  &
       wgt_idnf,eng_idnf,enn_idnf

  ! *** library subroutines ***

  PUBLIC set_usigmav_nf  ! set_usigmav_nf
  PUBLIC sigma_nf        ! sigma_nf(id_nf,mode_nf,energy) Reaction rate fit
  PUBLIC sigmav_nf       ! sigmav_nf(id_nf,mode_nf,temperature) Maxwellian fit
  PUBLIC sigmav_nf_local ! sigmavf_nf(X) Maxwellian local
  PUBLIC sigmav_nf_int   ! sigmav_nf(id_nf,mode_nf,temperature) Maxwellian int

CONTAINS
  
  ! --- set spline coefficients for reaction rate sigmav

  SUBROUTINE set_usigmav_nf
    USE plcomm
    USE libspl1d
    IMPLICIT NONE
    ! temperature range for reaction rate sigmav
    
    REAL(rkind),PARAMETER,DIMENSION(10):: &
         tempa=(/ 1.D0, 2.D0, 5.D0, 10.D0, 20.D0, &
                 50.D0, 100.D0, 200.D0, 500.D0, 1000.D0 /)
  
    ! reaction rate sigmav for DD
    REAL(rkind),PARAMETER,DIMENSION(10):: &
         svnf_dd=(/ 1.5D-22, 5.4D-21, 1.8D-19, 1.2D-18, 5.2D-18, &
                    2.1D-17, 4.5D-17, 8.8D-17, 1.8D-16, 2.2D-16 /)
    ! reaction rate sigmav for DT
    REAL(rkind),PARAMETER,DIMENSION(10):: &
         svnf_dt=(/ 5.5D-21, 2.6D-19, 1.3D-17, 1.1D-16, 4.2D-16, &
                    8.7D-16, 8.5D-16, 6.3D-16, 3.7D-16, 2.7D-16 /)
    ! reaction rate sigmav for DHe3
    REAL(rkind),PARAMETER,DIMENSION(10):: &
         svnf_dhe3=(/ 1.0D-26, 1.4D-23, 6.7D-21, 2.3D-19, 3.8D-18, &
                      5.4D-17, 1.6D-16, 2.4D-16, 2.3D-16, 1.8D-16 /)
    ! reaction rate sigmav for TT
    REAL(rkind),PARAMETER,DIMENSION(10):: &
         svnf_tt=(/ 3.3D-22, 7.1D-21, 1.4D-19, 7.2D-19, 2.5D-18, &
                    8.7D-18, 1.9D-17, 4.2D-17, 8.4D-17, 8.0D-17 /)
    ! reaction rate sigmav for THe3
    REAL(rkind),PARAMETER,DIMENSION(10):: &
         svnf_the3=(/ 1.0D-28, 1.0D-25, 2.1D-22, 1.2D-20, 2.6D-19, &
                      5.3D-18, 2.7D-17, 9.2D-17, 2.9D-16, 5.2D-16 /)

    INTEGER:: ntemp,id,ierr

    SELECT CASE(model_pnf)
    CASE(0) ! no fusion reaction
       nnfmax=0
       RETURN
    CASE(1) ! DT
       nnfmax=1
    CASE(2) ! DT+DD 
       nnfmax=4
    CASE(3) ! DT+DD+DHe3
       nnfmax=6
    CASE(4) ! DT+DD+DHe3+TT+THe3
       nnfmax=13
    CASE DEFAULT
       WRITE(6,*) 'XX Error libnf: undefined model_pnf: model_pnf=',model_pnf
       STOP
    END SELECT

    id_nf_nnf(1)=  id_nf_dt
    mode_nf_nnf(1)=1
    id_nf_nnf(2)=  id_nf_dd
    mode_nf_nnf(2)=1
    id_nf_nnf(3)=  id_nf_dd
    mode_nf_nnf(3)=2
    id_nf_nnf(4)=  id_nf_dd
    mode_nf_nnf(4)=3
    id_nf_nnf(5)=  id_nf_dhe3
    mode_nf_nnf(5)=1
    id_nf_nnf(6)=  id_nf_dhe3
    mode_nf_nnf(6)=2
    id_nf_nnf(7)=   id_nf_tt
    mode_nf_nnf(7)= 1
    id_nf_nnf(8)=   id_nf_the3
    mode_nf_nnf(8)= 1
    id_nf_nnf(9)=   id_nf_the3
    mode_nf_nnf(9)= 2
    id_nf_nnf(10)=  id_nf_the3
    mode_nf_nnf(10)=3
    id_nf_nnf(11)=  id_nf_the3
    mode_nf_nnf(11)=4
    id_nf_nnf(12)=  id_nf_the3
    mode_nf_nnf(12)=5
    id_nf_nnf(13)=  id_nf_the3
    mode_nf_nnf(13)=6

    SELECT CASE(model_pnf)
    CASE(1)
       IF(NS_D*NS_T*NS_He4.EQ.0) THEN
          IF(NS_D.EQ.0) WRITE(6,*)   'XX Error: libnf: NS_D=0'
          IF(NS_T.EQ.0) WRITE(6,*)   'XX Error: libnf: NS_T=0'
          IF(NS_He4.EQ.0) WRITE(6,*) 'XX Error: libnf: NS_He4=0'
          STOP
       END IF
    CASE(2,3,4)
       IF(NS_D*NS_T*NS_He4*NS_H*NS_He3.EQ.0) THEN
          IF(NS_D.EQ.0)   WRITE(6,*) 'XX Error: libnf: NS_D=0'
          IF(NS_T.EQ.0)   WRITE(6,*) 'XX Error: libnf: NS_T=0'
          IF(NS_H.EQ.0)   WRITE(6,*) 'XX Error: libnf: NS_H=0'
          IF(NS_He4.EQ.0) WRITE(6,*) 'XX Error: libnf: NS_He4=0'
          IF(NS_He3.EQ.0) WRITE(6,*) 'XX Error: libnf: NS_He3=0'
          STOP
       END IF
    CASE DEFAULT
       WRITE(6,*) 'XX Error libnb: undefined model_pnf: model_pnf=',model_pnf
       STOP
    END SELECT

    ns1_idnf(id_nf_dt,1)=NS_D
    ns2_idnf(id_nf_dt,1)=NS_T
    wgt_idnf(id_nf_dt,1)=1.0D0
    nsp_idnf(id_nf_dt,1)=NS_He4
    eng_idnf(id_nf_dt,1)=3.5D3*RKEV
    enn_idnf(id_nf_dt,1)=14.1D3*RKEV

    ns1_idnf(id_nf_dd,1)=NS_D
    ns2_idnf(id_nf_dd,1)=NS_D
    wgt_idnf(id_nf_dd,1)=0.5D0
    nsp_idnf(id_nf_dd,1)=NS_T
    eng_idnf(id_nf_dd,1)=1.01D3*RKEV
    enn_idnf(id_nf_dd,1)=0.D0
    
    ns1_idnf(id_nf_dd,2)=NS_D
    ns2_idnf(id_nf_dd,2)=NS_D
    wgt_idnf(id_nf_dd,2)=0.5D0
    nsp_idnf(id_nf_dd,2)=NS_H
    eng_idnf(id_nf_dd,2)=3.02D3*RKEV
    enn_idnf(id_nf_dd,3)=0.D0

    ns1_idnf(id_nf_dd,3)=NS_D
    ns2_idnf(id_nf_dd,3)=NS_D
    wgt_idnf(id_nf_dd,3)=0.5D0
    nsp_idnf(id_nf_dd,3)=NS_He3
    eng_idnf(id_nf_dd,3)=0.82D3*RKEV
    enn_idnf(id_nf_dd,3)=2.45D3*RKEV

    IF(model_pnf.GE.3) THEN
       ns1_idnf(id_nf_dhe3,1)=NS_D
       ns2_idnf(id_nf_dhe3,1)=NS_He3
       wgt_idnf(id_nf_dhe3,1)=1.D0
       nsp_idnf(id_nf_dhe3,1)=NS_He4
       eng_idnf(id_nf_dhe3,1)=3.6D3*RKEV
       enn_idnf(id_nf_dhe3,1)=0.D0

       ns1_idnf(id_nf_dhe3,2)=NS_D
       ns2_idnf(id_nf_dhe3,2)=NS_He3
       wgt_idnf(id_nf_dhe3,2)=1.D0
       nsp_idnf(id_nf_dhe3,2)=NS_H
       eng_idnf(id_nf_dhe3,2)=14.7D3*RKEV
       enn_idnf(id_nf_dhe3,2)=0.D0

       ns1_idnf(id_nf_tt,1)=NS_T
       ns2_idnf(id_nf_tt,1)=NS_T
       wgt_idnf(id_nf_tt,1)=1.D0
       nsp_idnf(id_nf_tt,1)=NS_He4
       eng_idnf(id_nf_tt,1)=1.25D3*RKEV  ! 11.3MeV*0.25/2.25 (not uniq)
       enn_idnf(id_nf_tt,1)=10.05D3*RKEV ! 11.3Mev*2.00/2.25 (not uniq)

       ns1_idnf(id_nf_the3,1)=NS_T
       ns2_idnf(id_nf_the3,1)=NS_He3
       wgt_idnf(id_nf_the3,1)=0.51D0
       nsp_idnf(id_nf_the3,1)=NS_He4
       eng_idnf(id_nf_the3,1)=1.34D3*RKEV ! 12.1MeV*0.25/2.25 (not uniq)
       enn_idnf(id_nf_the3,1)=5.38D3*RKEV ! 12.1Mev*1.00/2.25 (not uniq)

       ns1_idnf(id_nf_the3,2)=NS_T
       ns2_idnf(id_nf_the3,2)=NS_He3
       wgt_idnf(id_nf_the3,2)=0.51D0
       nsp_idnf(id_nf_the3,2)=NS_H
       eng_idnf(id_nf_the3,2)=5.38D3*RKEV ! 12.1MeV*1.0/2.25
       enn_idnf(id_nf_the3,2)=0.D0

       ns1_idnf(id_nf_the3,3)=NS_T
       ns2_idnf(id_nf_the3,3)=NS_He3
       wgt_idnf(id_nf_the3,3)=0.43D0
       nsp_idnf(id_nf_the3,3)=NS_He4
       eng_idnf(id_nf_the3,3)=4.8D3*RKEV
       enn_idnf(id_nf_the3,3)=0.D0

       ns1_idnf(id_nf_the3,4)=NS_T
       ns2_idnf(id_nf_the3,4)=NS_He3
       wgt_idnf(id_nf_the3,4)=0.43D0
       nsp_idnf(id_nf_the3,4)=NS_D
       eng_idnf(id_nf_the3,4)=9.58D3*RKEV
       enn_idnf(id_nf_the3,4)=0.D0

       ns1_idnf(id_nf_the3,5)=NS_T
       ns2_idnf(id_nf_the3,5)=NS_He3
       wgt_idnf(id_nf_the3,5)=0.06D0
       nsp_idnf(id_nf_the3,5)=NS_He4
       eng_idnf(id_nf_the3,5)=0.38D3*RKEV ! 1.89 MeV * 0.25/1.25
       enn_idnf(id_nf_the3,5)=1.51D3*RKEV ! 1.89 MeV * 1.00/1.25

       ns1_idnf(id_nf_the3,6)=NS_T
       ns2_idnf(id_nf_the3,6)=NS_He3
       wgt_idnf(id_nf_the3,6)=0.06D0
       nsp_idnf(id_nf_the3,6)=NS_H
       eng_idnf(id_nf_the3,6)=9.46D3*RKEV
       enn_idnf(id_nf_the3,6)=0.D0
    END IF

    DO ntemp=1,10
       tempa_log(ntemp)=LOG10(tempa(ntemp))
       svnf_dt_log(ntemp)=LOG10(svnf_dt(ntemp))
       svnf_dd_log(ntemp)=LOG10(svnf_dd(ntemp))
       svnf_dhe3_log(ntemp)=LOG10(svnf_dhe3(ntemp))
       svnf_tt_log(ntemp)=LOG10(svnf_tt(ntemp))
       svnf_the3_log(ntemp)=LOG10(svnf_the3(ntemp))
    END DO
    
    CALL SPL1D(tempa_log,svnf_dt_log,  dsvnf,usvnf_dt,  10,0,ierr)
    id=1
    IF(ierr.EQ.0) THEN
       CALL SPL1D(tempa_log,svnf_dd_log,  dsvnf,usvnf_dd,  10,0,ierr)
       id=2
    ENDIF
    IF(ierr.EQ.0) THEN
       CALL SPL1D(tempa_log,svnf_dhe3_log,dsvnf,usvnf_dhe3,10,0,ierr)
       id=3
    END IF
    IF(ierr.EQ.0) THEN
       CALL SPL1D(tempa_log,svnf_tt_log,  dsvnf,usvnf_tt,  10,0,ierr)
       id=4
    END IF
    IF(ierr.EQ.0) THEN
       CALL SPL1D(tempa_log,svnf_the3_log,dsvnf,usvnf_the3,10,0,ierr)
       id=5
    END IF
    IF(ierr.NE.0) THEN
       WRITE(6,'(A,I4)') 'XX SPL1D error in set_usvnf: id=',id
       STOP
    END IF

    RETURN
  END SUBROUTINE set_usigmav_nf

  ! --- cross section of nuclear fusion reaction ---
  ! ---     in barn (10^{-28}m^{-2})
  ! ---     as a function of energy in keV
  ! ---     NRL Plasma Formulary 2019

  FUNCTION sigma_nf(id_nf,mode_nf,energy)

    IMPLICIT NONE
    
    ! Duane coef (NRL Formulary 2019)
    REAL(rkind),PARAMETER,DIMENSION(5,6):: &
         Duane=reshape((/46.097D0, 372.D0, 4.36D-4,  1.220D0, 0.D0, &
                         47.88D0,  482.D0, 3.08D-4,  1.177D0, 0.D0, &
                         45.95D0,  5.02D4, 1.368D-2, 1.076D0, 409.D0, &
                         89.27D0,  2.59D4, 3.98D-3,  1.297D0, 647.D0, &
                         38.39D0,  448.D0, 1.02D-3,  2.09D0,  0.D0, &
                        123.1D0, 1.125D4, 0.D0,     0.D0,    0.D0/), &
                       (/5,6/))

    INTEGER,INTENT(IN):: id_nf,mode_nf
    REAL(rkind),INTENT(IN):: energy
    INTEGER:: id
    REAL(rkind):: sigma_nf

!    WRITE(6,*) '@@@ point 71:',id_nf,mode_nf,energy
    IF(id_nf.LT.1.OR.id_nf.GT.6) THEN
       WRITE(6,'(A,I4)') &
            'XX sigma_nf: input error: undefined id_nf: ',id_nf
       STOP
    ENDIF
    IF(mode_nf.LT.1.OR.mode_nf.GT.mode_nf_max(id_nf)) THEN
       WRITE(6,'(A,I4)') &
            'XX sigma_nf: input error: undefined mode_nf: ',mode_nf
       STOP
    ENDIF
    IF(energy.LE.0.D0) THEN
       WRITE(6,*) 'XX sigma_duane: input error: non-positive energy: ',energy
       STOP
    ENDIF

    SELECT CASE(id_nf)
    CASE(id_nf_DT) ! DT
       id=3
    CASE(id_nf_DD) ! DD
       SELECT CASE(mode_nf)
       CASE(1)
          id=1
       CASE(2)
          id=2
       END SELECT
    CASE(id_nf_DHe3) ! DHe3
       id=4
    CASE(id_nf_TT)   ! TT
       id=5
    CASE(id_nf_THe3) ! THe3
       id=6
    END SELECT
    
!    WRITE(6,*) '@@@ point 72:',id
    sigma_nf=(Duane(5,id) &
         +Duane(2,id) &
         /((Duane(4,id)-Duane(3,id)*energy)**2+1.D0)) &
         /(energy*(EXP(Duane(1,id)/SQRT(energy))-1.D0))
!    WRITE(6,*) '@@@ point 73:',sigma_nf
    RETURN
  END FUNCTION sigma_nf
  
  ! --- reaction rate of nuclear fusion: sigmav  ---
  ! ---     as a function of temperature in keV

  FUNCTION sigmav_nf(id_nf,temperature)

    USE libspl1d
    IMPLICIT NONE
    
    INTEGER,INTENT(IN):: id_nf
    REAL(rkind),INTENT(IN):: temperature
    REAL(rkind):: sigmav_nf,temperature_log,temperature_local,sigmav_log
    INTEGER:: ierr
    
    IF(id_nf.LT.1.OR.id_nf.GT.6) THEN
       WRITE(6,'(A,I4)') 'XX sigmav_nf: input error: undefined id_nf: ',id_nf
       STOP
    END IF

    temperature_local=temperature
    IF(temperature.LT.1.D0) temperature_local=1.D0
    IF(temperature.GT.1000.D0) temperature_local=1000.D0

    temperature_log=LOG10(temperature_local)
    SELECT CASE(id_nf)
    CASE(id_nf_dt)
       CALL SPL1DF(temperature_log,sigmav_log,tempa_log,usvnf_dt,10,ierr)
    CASE(id_nf_dd)
       CALL SPL1DF(temperature_log,sigmav_log,tempa_log,usvnf_dd,10,ierr)
    CASE(id_nf_dhe3)
       CALL SPL1DF(temperature_log,sigmav_log,tempa_log,usvnf_dhe3,10,ierr)
    CASE(id_nf_tt)
       CALL SPL1DF(temperature_log,sigmav_log,tempa_log,usvnf_tt,10,ierr)
    CASE(id_nf_the3)
       CALL SPL1DF(temperature_log,sigmav_log,tempa_log,usvnf_the3,10,ierr)
    END SELECT
    IF(ierr.NE.0) THEN
       WRITE(6,'(A,2I4)')     'XX SPL1DF error in sigmav_nf: id_nf,ierr=', &
            id_nf,ierr
       WRITE(6,'(A,ES12.4)') '       temperature=',temperature
       STOP
    END IF
    sigmav_nf=10.D0**sigmav_log*1.D-6   ! cm^3/s -> m^3/s
  END FUNCTION sigmav_nf

  ! --- sigmav_nf by integeral over energy ---
  !         <sigmav>=\int_0^\infty 4\pi v^2 dv sigma v f(v)
  !                  f(v)=(m/2\pi T)^{3/2} exp(-mv^2/2T)   normalized to 1
  !                  E=mv^2/2
  !                  v=SQRT{2E/m}
  !                  dv=(1/2) SQRT{2/mE} dE
  !                  4\pi v^2 dv=4\pi (2E/m) (1/2) SQRT{2/mE} dE
  !                             =4\pi SQRT{E^2/m^2 2/mE} dE
  !                             =4\pi SQRT{2E/m^3} dE
  !                  f(E)=(m/2\pi T)^{3/2} exp(-E/T)
  !                <f(E)>=\int_0^\infty 4\pi v^2 dv f(v) 
  !                      =\int_0^\infty 4\pi SQRT{2E/m^3} dE
  !                               (m/2\pi T)^{3/2} exp(-E/T)
  !                      =\int_0^\infty SQRT{16 \pi^2 2E/m^3 m^3/8 \pi^3 T^3}
  !                               dE exp(-E/T)
  !                      =\int_0^\infty SQRT{4E/\pi T} exp(-E/T) dE/T
  !                      =\int_0^\infty SQRT{4X/\pi} exp(-X) dX
  !                      = (2/SQRT{\pi}) \int_0^\infty SQRT{X} exp(-X) dX
  !                      = 1
  !          <sigmav>=\int_0^\infty 4\pi SQRT{2E/m^3} dE sigma SQRT{2E/m} f(E)
  !                  =\int_0^\infty 4\pi 2E/m^2 dE sigma
  !                                   (m/2\pi T)^{3/2} exp(-E/T)
  !                  =\int_0^\infty SQRT{8 E^2/m\pi T} sigma exp(-E/T) dE/T
  !                  =\int_0^\infty SQRT{8T/m\pi} (E/T) sigma exp(-E/T) dE/T
  !                  =\int_0^\infty SQRT{8T/m\pi} X sigma exp(-X) dX

  ! --- sigmav for energy --- X=energy/temperature

  FUNCTION sigmav_nf_local(X)
    USE plcomm
    USE libnf_local
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: X
    REAL(rkind):: sigmav_nf_local
    REAL(rkind):: energy,velocity,sigma

    energy=temperature_local*X
    velocity=SQRT(2.D0*energy*RKEV/pm_local)
    sigma=sigma_nf(id_nf_local,mode_nf_local,energy)

!    sigmav_nf_local=SQRT(8.D0*temperature_local/(pi*pm_local)) &
!         *X*sigma*EXP(-X)*1.D-6*1.D-28      ! cm^3/s -> m^3/s, barn
    sigmav_nf_local=SQRT(2.D0/pi)*EXP(-X)
    WRITE(6,*) '@@@ point 1.9:',X,sigmav_nf_local

!    WRITE(6,'(A,5ES12.4)') '@@@ point 2:', &
!         X,energy,velocity,sigma,sigmav_nf_local
    RETURN
  END FUNCTION sigmav_nf_local


  FUNCTION sigmav_nf_int(id_nf,mode_nf,temperature)

    USE plcomm
    USE libnf_local
    USE libde
    IMPLICIT NONE
    INTEGER,INTENT(IN):: id_nf,mode_nf
    REAL(rkind),INTENT(IN):: temperature
    REAL(rkind):: sigmav_nf_int
    REAL(rkind):: error_int,H0,EPS
    INTEGER:: ILST

    IF(id_nf.LT.1.OR.id_nf.GT.6) THEN
       WRITE(6,'(A,I4)') 'XX sigmav_nf: input error: undefined id_nf: ',id_nf
       STOP
    END IF

    id_nf_local=id_nf
    mode_nf_local=mode_nf
    temperature_local=temperature
    IF(temperature_local.LT.1.D0) temperature_local=1.D0
    IF(temperature_local.GT.1000.D0) temperature_local=1000.D0
    pm_local=PA(nsp_idnf(id_nf_local,mode_nf_local))*AMP

    H0=1.D-2
    EPS=1.D-4
    ilst=1
    WRITE(6,'(A,2ES12.4)') '@@@ point 1:',temperature_local,pm_local
    
    CALL DEHIFE(sigmav_nf_int,error_int,H0,eps,ilst,sigmav_nf_local, &
         'sigmav_nf_int')

    RETURN
  END FUNCTION sigmav_nf_int

END MODULE libnf
      
