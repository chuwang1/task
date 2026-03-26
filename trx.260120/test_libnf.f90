PROGRAM test_libnf
  USE plcomm
  USE plinit
  USE libnf
  USE libgrf
  IMPLICIT NONE
  REAL(rkind),ALLOCATABLE:: xg(:),yg(:,:)
  REAL(rkind):: Es,Ee,delE,En
  INTEGER:: nmax,n

  CALL GSOPEN
  CALL pl_init
  model_pnf=3
  NS_D=2
  NS_T=3
  NS_He4=4
  NS_He3=5
  NS_H=6
  CALL set_usigmav_nf

  CALL PAGES
!  CALL linbnf_gr1(1)
!  CALL linbnf_gr2(2)
  CALL linbnf_gr3(3)
  CALL linbnf_gr4(4)
  CALL PAGEE

  CALl GSCLOS
  STOP
  
CONTAINS

  ! --- sigma (E) ---
  
  SUBROUTINE linbnf_gr1(idg)

    USE plcomm
    USE libnf
    USE libgrf
    IMPLICIT NONE
    INTEGER,INTENT(IN):: idg
    REAL(rkind),ALLOCATABLE:: xg(:),yg(:,:)
    REAL(rkind):: Es,Ee,delE,En
    INTEGER:: nmax,n

    Es=LOG10(10.D0)  ! start energy in keV
    Ee=LOG10(10000.D0)  ! end energy in keV
    nmax=201
    delE=(Ee-Es)/(nmax-1)

    ALLOCATE(xg(nmax),yg(nmax,6))

    DO n=1,nmax
       En=Es+(n-1)*delE
       xg(n)=En
       yg(n,1)=LOG10(sigma_nf(id_nf_DT,1,  10.D0**En))
       yg(n,2)=LOG10(sigma_nf(id_nf_DD,1,  10.D0**En))
       yg(n,3)=LOG10(sigma_nf(id_nf_DD,2,  10.D0**En))
       yg(n,4)=LOG10(sigma_nf(id_nf_DHe3,1,10.D0**En))
       yg(n,5)=LOG10(sigma_nf(id_nf_TT,1,  10.D0**En))
       yg(n,6)=LOG10(sigma_nf(id_nf_THe3,1,10.D0**En))
    END DO

    CALL GRD1D(idg,xg,yg,nmax,nmax,6,&
         '@sigma [barn=10^{-28}m^2] vs E [keV]@',3,fmin=LOG10(1.D-3))

    RETURN
  END SUBROUTINE linbnf_gr1

  ! --- sigma v (E) ---
  
  SUBROUTINE linbnf_gr2(idg)

    USE plcomm
    USE libnf
    USE libgrf
    IMPLICIT NONE
    INTEGER,INTENT(IN):: idg
    REAL(rkind),ALLOCATABLE:: xg(:),yg(:,:)
    REAL(rkind):: Es,Ee,delE,En
    INTEGER:: nmax,n

    Es=LOG10(1.D0)  ! start energy in keV
    Ee=LOG10(1000.D0)  ! end energy in keV
    nmax=201
    ALLOCATE(xg(nmax),yg(nmax,5))
    delE=(Ee-Es)/(nmax-1)

    DO n=1,nmax
       En=Es+(n-1)*delE
       xg(n)=En
       yg(n,1)=LOG10(sigmav_nf(id_nf_DT,  10.D0**En))
       yg(n,2)=LOG10(sigmav_nf(id_nf_DD,  10.D0**En))
       yg(n,3)=LOG10(sigmav_nf(id_nf_DHe3,10.D0**En))
       yg(n,4)=LOG10(sigmav_nf(id_nf_TT,  10.D0**En))
       yg(n,5)=LOG10(sigmav_nf(id_nf_THe3,10.D0**En))
    END DO

    CALL GRD1D(idg,xg,yg,nmax,nmax,5, &
         '@sigmav [m^3sec^{-1}] vs T [keV]@',3,fmin=-27.D0)
    RETURN
  END SUBROUTINE linbnf_gr2

  ! --- sigma v f (E) ---
  
  SUBROUTINE linbnf_gr3(idg)

    USE plcomm
    USE libnf
    USE libgrf
    USE libnf_local
    IMPLICIT NONE
    INTEGER,INTENT(IN):: idg
    REAL(rkind),ALLOCATABLE:: xg(:),yg(:,:)
    REAL(rkind):: Es,Ee,delE,xl
    INTEGER:: nmax,n

    Es=LOG10(1.D0)  ! start energy in keV
    Ee=LOG10(1000.D0)  ! end energy in keV
    nmax=201
    ALLOCATE(xg(nmax),yg(nmax,5))
    delE=(Ee-Es)/(nmax-1)

    id_nf_local=id_nf_DD
    mode_nf_local=1
    temperature_local=1000.D0
    pm_local=PA(nsp_idnf(id_nf_local,mode_nf_local))*AMP
    
    DO n=1,nmax
       xl=DBLE(n)/DBLE(nmax)
       xg(n)=xl
       id_nf_local=id_nf_DT
       yg(n,1)=LOG10(sigmav_nf_local(xl))
       id_nf_local=id_nf_DD
       yg(n,2)=LOG10(sigmav_nf_local(xl))
       id_nf_local=id_nf_DHe3
       yg(n,3)=LOG10(sigmav_nf_local(xl))
       id_nf_local=id_nf_TT
       yg(n,4)=LOG10(sigmav_nf_local(xl))
       id_nf_local=id_nf_THe3
       yg(n,5)=LOG10(sigmav_nf_local(xl))
    END DO

    CALL GRD1D(idg,xg,yg,nmax,nmax,5, &
         '@sigmavf vs E [keV] for 10keV@',3,fmin=-27.D0)

    RETURN
  END SUBROUTINE linbnf_gr3
  
  ! --- \int sigma v f (E) dE---
  
  SUBROUTINE linbnf_gr4(idg)

    USE plcomm
    USE libnf
    USE libgrf
    IMPLICIT NONE
    INTEGER,INTENT(IN):: idg
    REAL(rkind),ALLOCATABLE:: xg(:),yg(:,:)
    REAL(rkind):: Es,Ee,delE,En
    INTEGER:: nmax,n

    Es=LOG10(1.D0)  ! start energy in keV
    Ee=LOG10(1000.D0)  ! end energy in keV
    nmax=201
    ALLOCATE(xg(nmax),yg(nmax,6))
    delE=(Ee-Es)/(nmax-1)

    DO n=1,nmax
       En=Es+(n-1)*delE
       xg(n)=En
       yg(n,1)=LOG10(sigmav_nf_int(id_nf_DT,  1,10.D0**En))
       yg(n,2)=LOG10(sigmav_nf_int(id_nf_DD,  1,10.D0**En))
       yg(n,3)=LOG10(sigmav_nf_int(id_nf_DHe3,1,10.D0**En))
       yg(n,4)=LOG10(sigmav_nf_int(id_nf_TT,  1,10.D0**En))
       yg(n,5)=LOG10(sigmav_nf_int(id_nf_THe3,1,10.D0**En))
       WRITE(6,'(I4,6ES12.4)') n,xg(n),yg(n,1),yg(n,2),yg(n,3),yg(n,4),yg(n,5)
    END DO

    CALL GRD1D(idg,xg,yg,nmax,nmax,5, &
         '@sigmav [m^3sec^{-1}] vs T [keV]@',3,fmin=-27.D0)

    RETURN
  END SUBROUTINE linbnf_gr4
END PROGRAM test_libnf

