! libnfPB.f90

MODULE libnfPB
  
  USE bpsd_kinds
  USE bpsd_constants
  IMPLICIT NONE
  PRIVATE
  PUBLIC sigma_nf_PB_NS
  PUBLIC sigma_nf_PB_SW
  PUBLIC S_nf_PB_NS
  PUBLIC S_nf_PB_SW

CONTAINS
  
  ! --- p-B reaction ---
  !        Ref. A. Tantori and F Belloni, Nucl. Fusion 63 (2023) 086001 (9pp)
  !     cross section --- sigma
  !     astrophysics factor --- S

  FUNCTION sigma_nf_PB_NS(E_Mev)

    USE plcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: E_MeV ! Energy in MeV
    REAL(rkind):: sigma_nf_PB_NS

    sigma_nf_PB_NS=S_nf_PB_NS(E_MeV*1.D3)/E_MeV*EXP(-SQRT(22.589D0/E_MeV))
    RETURN
  END FUNCTION sigma_nf_PB_NS
    
  FUNCTION sigma_nf_PB_SW(E_Mev)

    USE plcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: E_MeV ! Energy in MeV
    REAL(rkind):: sigma_nf_PB_SW

    sigma_nf_PB_SW=S_nf_PB_SW(E_MeV*1.D3)/E_MeV*EXP(-SQRT(22.589D0/E_MeV))
    RETURN
  END FUNCTION sigma_nf_PB_SW
    
  FUNCTION S_nf_PB_NS(E_keV)

    USE plcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: E_keV ! Energy in keV
    REAL(rkind):: S_nf_PB_NS
    ! constants all in MeV
    REAL(rkind),PARAMETER:: C_0=197.D0
    REAL(rkind),PARAMETER:: C_1=0.240D0
    REAL(rkind),PARAMETER:: C_2=2.31D-4
    REAL(rkind),PARAMETER:: A_L=1.82D4
    REAL(rkind),PARAMETER:: E_L=148.0D-3
    REAL(rkind),PARAMETER:: delE_L=2.35D-3
    REAL(rkind),PARAMETER:: D_0=330.D0
    REAL(rkind),PARAMETER:: D_1=66.1D0
    REAL(rkind),PARAMETER:: D_2=-20.3D0
    REAL(rkind),PARAMETER:: D_5=-1.58D0
    REAL(rkind),PARAMETER:: A_0=2.57D6
    REAL(rkind),PARAMETER:: A_1=5.67D5
    REAL(rkind),PARAMETER:: A_2=1.34D5
    REAL(rkind),PARAMETER:: A_3=5.68D5
    REAL(rkind),PARAMETER:: E_0=581.3D-3
    REAL(rkind),PARAMETER:: E_1=1083.D-3
    REAL(rkind),PARAMETER:: E_2=2405.D-3
    REAL(rkind),PARAMETER:: E_3=3344.D-3
    REAL(rkind),PARAMETER:: delE_0=85.7D-3
    REAL(rkind),PARAMETER:: delE_1=234.D-3
    REAL(rkind),PARAMETER:: delE_2=138.D-3
    REAL(rkind),PARAMETER:: delE_3=309.D-3
    REAL(rkind),PARAMETER:: B=4.38D0
  
    REAL(rkind):: E_MeV, E_n, S

    E_MeV=E_kev*0.001D0 ! Energy in  MeV

    IF(E_MeV.LE.0.4D0) THEN ! S_1
       S=C_0+C_1*E_keV+C_2*E_keV**2+A_L*1.D-6/((E_MeV-E_L)**2+delE_L**2)
    ELSE IF(E_MeV.LE.0.642D0) THEN ! S_2
       E_n=1.D1*(E_MeV-0.400D0)
       S=D_0+D_1*E_n+D_2*E_n**2+D_5*E_n**5
    ELSE IF(E_MeV.LE.3.5D0) THEN ! S_3
       S=B+A_0*1.D-6/((E_MeV-E_0)**2+delE_0**2) &
          +A_1*1.D-6/((E_MeV-E_1)**2+delE_1**2) &
          +A_2*1.D-6/((E_MeV-E_2)**2+delE_2**2) &
          +A_3*1.D-6/((E_MeV-E_3)**2+delE_3**2)
    ELSE ! out of range
       S=B+A_0*1.D-6/((3.5D0-E_0)**2+delE_0**2) &
          +A_1*1.D-6/((3.5D0-E_1)**2+delE_1**2) &
          +A_2*1.D-6/((3.5D0-E_2)**2+delE_2**2) &
          +A_3*1.D-6/((3.5D0-E_3)**2+delE_3**2)
    END IF
    S_nf_PB_NS=S
    RETURN
  END FUNCTION S_nf_PB_NS

  FUNCTION S_nf_PB_SW(E_keV)

    USE plcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: E_keV ! Energy in keV
    REAL(rkind):: S_nf_PB_SW
    ! constants all in MeV
    REAL(rkind),PARAMETER:: C_0=197.D0
    REAL(rkind),PARAMETER:: C_1=0.269D0
    REAL(rkind),PARAMETER:: C_2=2.54D-4
    REAL(rkind),PARAMETER:: D_0=346.D0
    REAL(rkind),PARAMETER:: D_1=150.D0
    REAL(rkind),PARAMETER:: D_2=-59.9D0
    REAL(rkind),PARAMETER:: D_5=-0.460D0
    REAL(rkind),PARAMETER:: A_0=1.98D6
    REAL(rkind),PARAMETER:: A_1=3.89D6
    REAL(rkind),PARAMETER:: A_2=1.36D6
    REAL(rkind),PARAMETER:: A_3=3.71D6
    REAL(rkind),PARAMETER:: E_0=640.9D-3
    REAL(rkind),PARAMETER:: E_1=1211.D-3
    REAL(rkind),PARAMETER:: E_2=2340.D-3
    REAL(rkind),PARAMETER:: E_3=3294.D-3
    REAL(rkind),PARAMETER:: delE_0=85.5D-3
    REAL(rkind),PARAMETER:: delE_1=414.D-3
    REAL(rkind),PARAMETER:: delE_2=221.D-3
    REAL(rkind),PARAMETER:: delE_3=351.D-3
    REAL(rkind),PARAMETER:: B=0.381D0
  
    REAL(rkind):: E_MeV, E_n, S

    E_MeV=E_kev*0.001D0 ! Energy in  MeV

    IF(E_MeV.LE.0.4D0) THEN ! S_1
       S=C_0+C_1*E_keV+C_2*E_keV**2
    ELSE IF(E_MeV.LE.0.668D0) THEN ! S_2
       E_n=1.D1*(E_MeV-0.400D0)
       S=D_0+D_1*E_n+D_2*E_n**2+D_5*E_n**5
    ELSE IF(E_MeV.LE.9.76D0) THEN ! S_3
       S=B+A_0*1.D-6/((E_MeV-E_0)**2+delE_0**2) &
          +A_1*1.D-6/((E_MeV-E_1)**2+delE_1**2) &
          +A_2*1.D-6/((E_MeV-E_2)**2+delE_2**2) &
          +A_3*1.D-6/((E_MeV-E_3)**2+delE_3**2)
    ELSE ! out of range
       S=B+A_0*1.D-6/((9.76D0-E_0)**2+delE_0**2) &
          +A_1*1.D-6/((9.76D0-E_1)**2+delE_1**2) &
          +A_2*1.D-6/((9.76D0-E_2)**2+delE_2**2) &
          +A_3*1.D-6/((9.76D0-E_3)**2+delE_3**2)
    END IF
    S_nf_PB_SW=S
    RETURN
  END FUNCTION S_nf_PB_SW

END MODULE libnfPB
