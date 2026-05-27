! trfixed.f90

MODULE trfixed

  USE trcomm,ONLY: rkind


  !  *** Define fixed profile of density and temperature ***

  INTEGER:: ntime_nfixed_max                ! number of time points
  INTEGER:: ndata_nfixed_max                ! number of coef data
  REAL(rkind):: rho_min_nfixed,rho_max_nfixed  ! range of fixed profile
  REAL(rkind),ALLOCATABLE:: time_nfixed(:)     ! time points t_i
  REAL(rkind),ALLOCATABLE:: coef_nfixed(:,:)   ! coef data for t_i<= t <t_{i+1}
                                            ! temperature profile
  INTEGER:: ntime_tfixed_max                ! number of time points
  INTEGER:: ndata_tfixed_max                ! number of coef data
  REAL(rkind) :: rho_min_tfixed,rho_max_tfixed  ! range of fixed profile
  REAL(rkind),ALLOCATABLE:: time_tfixed(:)     ! time points t_i
  REAL(rkind),ALLOCATABLE:: coef_tfixed(:,:)   ! coef data for t_i<= t <t_{i+1}

  PUBLIC tr_set_nfixed  ! set coef matrix for n
  PUBLIC tr_set_tfixed  ! set coef matrix for nT
  PUBLIC tr_prof_nfixed ! set fixed density profile
  PUBLIC tr_prof_tfixed ! set fixed temperature profile
  PUBLIC tr_prep_nfixed ! read fixed density pfofile parameters
  PUBLIC tr_prep_tfixed ! read fixed temperature profile parameters
  PUBLIC tr_prep_arfixed ! read Ar density from CSV file
  PUBLIC tr_prep_radfixed ! read qrad from CSV file
  PUBLIC tr_prep_prlfixed ! read PRL from CSV file
  PUBLIC tr_prep_prffixed ! read qrfe from CSV file
  PUBLIC tr_prep_chifixed ! read chi from external file

CONTAINS

  !     ***** Routine for fixed density profile *****
      
  SUBROUTINE tr_set_nfixed(nr,time)

    USE trcomm
    USE trcomx
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nr
    REAL(rkind):: rn_local
    REAL(rkind),INTENT(IN):: time
    INTEGER:: NS,NEQ,NW

    IF(model_nfixed.EQ.0) RETURN
    IF(time.LE.time_nfixed(1)) return
    IF(model_nfixed.EQ.2) THEN
       IF((rm(nr).LT.rho_min_nfixed).OR. &
            (rm(nr).GT.rho_max_nfixed)) RETURN
    END IF
    CALL tr_prof_nfixed(rm(nr),time,rn_local)
    NEQ=NEA(1,1) ! NEQ of electron density equation
    DO NW=1,NEQMAX
       A(NEQ,NW,NR) = 0.D0
       B(NEQ,NW,NR) = 0.D0
       C(NEQ,NW,NR) = 0.D0
    END DO
    D(NEQ,NR)=0.D0
!    B(NEQ,NEQ,NR)=-1.D0/tau_nfixed
!    D(NEQ,NR)=rn_local/tau_nfixed
    RD(NEQ,NR)=1.D0
    DO NS=2,NSMAX
       NEQ=NEA(NS,1) ! NEQ of density equation
       DO NW=1,NEQMAX
          A(NEQ,NW,NR) = 0.D0
          B(NEQ,NW,NR) = 0.D0
          C(NEQ,NW,NR) = 0.D0
       END DO
       D(NEQ,NR)=0.D0
!       B(NEQ,NEQ,NR)=-1.D0/tau_nfixed
!       D(NEQ,NR)=pn(ns)/(pz(ns)*pn(1))*rn_local/tau_nfixed
       RD(NEQ,NR)=1.D0
    END DO
    RETURN
  END SUBROUTINE tr_set_nfixed

  !     ***** Routine for fixed temperature profile *****
      
  SUBROUTINE tr_set_tfixed(nr,time)

    USE trcomm
    USE trcomx
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nr
    REAL(rkind),INTENT(IN):: time
    REAL(rkind):: rt_local
    INTEGER:: NS,NEQ,NW

    IF(model_tfixed.EQ.0) RETURN
    IF(time.LE.time_tfixed(1)) return
    IF(model_tfixed.EQ.2) THEN
       IF((rm(nr).LT.rho_min_tfixed).OR. &
            (rm(nr).GT.rho_max_tfixed)) RETURN
    END IF
    CALL tr_prof_tfixed(rm(nr),time,rt_local)
    NEQ=NEA(1,1) ! NEQ of electron density equation
    DO NW=1,NEQMAX
       A(NEQ,NW,NR) = 0.D0
       B(NEQ,NW,NR) = 0.D0
       C(NEQ,NW,NR) = 0.D0
    END DO
    D(NEQ,NR)=0.D0
!    B(NEQ,NEQ,NR)=-1.D0/tau_tfixed
!    D(NEQ,NR)=rt_local/tau_tfixed
    RD(NEQ,NR)=1.D0
    DO NS=2,NSMAX
       NEQ=NEA(NS,2) ! NEQ of temperature equation
       DO NW=1,NEQMAX
          A(NEQ,NW,NR) = 0.D0
          B(NEQ,NW,NR) = 0.D0
          C(NEQ,NW,NR) = 0.D0
       END DO
       D(NEQ,NR)=0.D0
!       B(NEQ,NEQ,NR)=-1.D0/tau_tfixed
!       D(NEQ,NR)=rt_local/tau_tfixed
       RD(NEQ,NR)=1.D0
    END DO
    RETURN
  END SUBROUTINE tr_set_tfixed
      
  ! *** set fixed density profile ***
  
  SUBROUTINE tr_prof_nfixed(rho,time,rn_local)
  
    USE trcomm
    IMPLICIT NONE    
    REAL(rkind),INTENT(IN):: rho,time
    REAL(rkind),INTENT(OUT):: rn_local
    REAL(rkind):: tr_func_nfixed
    REAL(rkind),ALLOCATABLE:: coef(:)
    REAL(rkind):: factor
    INTEGER:: id,i,ntime

    ! --- find time range ---

    IF(time.LT.time_nfixed(1)) THEN
       RETURN
    ELSE IF (time.GE.time_nfixed(ntime_nfixed_max)) THEN
       id=ntime_nfixed_max
    ELSE
       DO ntime=1,ntime_nfixed_max-1
          IF(time.GE.time_nfixed(ntime).AND. &
               time.LT.time_nfixed(ntime+1)) THEN
             id=ntime
          END IF
       END DO
    END IF

    ! --- set profile coefficients ---
    
    ALLOCATE(coef(0:ndata_nfixed_max))
    IF(id.EQ.ntime_nfixed_max) THEN ! after time_nfixed(ntime_nfixed_max)
       DO i=0,ndata_nfixed_max
          coef(i)=coef_nfixed(i,ntime_nfixed_max)
       END DO
    ELSE ! between time_nfixed(id) and time_nfixed(id+1)
       factor=(time-time_nfixed(id)) &
             /(time_nfixed(id+1)-time_nfixed(id))
       DO i=0,ndata_nfixed_max
          coef(i)=(1.D0-factor)*coef_nfixed(i,id) &
                        +factor*coef_nfixed(i,id+1)
       END DO
    END IF

    ! --- set local density profile ---
    
    rn_local=coef(0) &
         +0.5D0*coef(1) &
         *(tanh((1.D0-coef(2)*coef(3)-rho)/coef(3))+1.D0) &
         +coef(4)*(1.D0-rho*rho)**coef(5) &
         +0.5D0*coef(8)*(1.D0-erf((rho-coef(9))/SQRT(2.D0*coef(10))))
    rn_local=rn_local*1.D-20
    IF(rn_local.LE.0.D0) rn_local=1.D-8
    RETURN
  END SUBROUTINE tr_prof_nfixed

  ! *** set fixed temperature profile ***
  
  SUBROUTINE tr_prof_tfixed(rho,time,rt_local)
  
    USE trcomm
    IMPLICIT NONE    
    REAL(rkind),INTENT(IN):: rho,time
    REAL(rkind),INTENT(OUT):: rt_local
    REAL(rkind):: tr_func_tfixed
    REAL(rkind),ALLOCATABLE:: coef(:)
    REAL(rkind):: factor
    INTEGER:: id,i,ntime

    ! --- find time range ---

    IF(time.LE.time_tfixed(1)) THEN
       RETURN
    ELSE IF (time.GE.time_tfixed(ntime_tfixed_max)) THEN
       id=ntime_tfixed_max
    ELSE
       DO ntime=1,ntime_tfixed_max-1
          IF(time.GE.time_tfixed(ntime).AND. &
               time.LE.time_tfixed(ntime+1)) THEN
             id=ntime
          END IF
       END DO
    END IF

    ! --- set profile coefficients ---
    
    ALLOCATE(coef(0:ndata_tfixed_max))
    IF(id.EQ.0) THEN ! before time_nfixed(1)
       DO i=0,ndata_tfixed_max
          coef(i)=coef_tfixed(i,1)
       END DO
    ELSE IF(id.EQ.ntime_tfixed_max) THEN ! after time_nfixed(ntime_nfixed_max)
       DO i=0,ndata_tfixed_max
          coef(i)=coef_tfixed(i,ntime_tfixed_max)
       END DO
    ELSE ! between time_nfixed(id) and time_nfixed(id+1)
       factor=(time-time_tfixed(id)) &
             /(time_tfixed(id+1)-time_tfixed(id))
       DO i=0,ndata_tfixed_max
          coef(i)=(1.D0-factor)*coef_tfixed(i,id) &
                        +factor*coef_tfixed(i,id+1)
       END DO
    END IF

    ! --- set temperature profile ---
    
    rt_local=coef(0) &
         +0.5D0*coef(1) &
         *(tanh((1.D0-coef(2)*coef(3)-rho)/coef(3))+1.D0) &
         +coef(4)*(1.D0-rho*rho)**coef(5) &
         +0.5D0*coef(8)*(1.D0-erf((rho-coef(9))/SQRT(2.D0*coef(10))))
    rt_local=rt_local*1.D-3
    IF(rt_local.LE.0.D0) rt_local=3.D-5
    RETURN
  END SUBROUTINE tr_prof_tfixed

  ! *** read density profile data from file ***

  SUBROUTINE tr_prep_nfixed
    USE trcomm
    USE libfio
    IMPLICIT NONE
    INTEGER:: nfl,ntime,ndata,ierr

    NFL=12
    CALL fropen(NFL,knam_nfixed,1,0,'fn',ierr)
    READ(NFL,*) ntime_nfixed_max,ndata_nfixed_max,rho_min_nfixed,rho_max_nfixed
    IF(ALLOCATED(time_nfixed)) DEALLOCATE(time_nfixed)
    IF(ALLOCATED(coef_nfixed)) DEALLOCATE(coef_nfixed)
    ALLOCATE(time_nfixed(ntime_nfixed_max))
    ALLOCATE(coef_nfixed(0:ndata_nfixed_max,ntime_nfixed_max))
    DO ntime=1,ntime_nfixed_max
       READ(NFL,*) time_nfixed(ntime)
       READ(NFL,*) (coef_nfixed(ndata,ntime),ndata=0,ndata_nfixed_max)
    END DO
    CLOSE(NFL)
    RETURN
  END SUBROUTINE tr_prep_nfixed
    
  ! *** read temperature profile data from file ***

  SUBROUTINE tr_prep_tfixed
    USE trcomm
    USE libfio
    IMPLICIT NONE
    INTEGER:: nfl,ntime,ndata,ierr

    NFL=12
    CALL fropen(NFL,knam_tfixed,1,0,'ft',ierr)
    READ(NFL,*) ntime_tfixed_max,ndata_tfixed_max,rho_min_tfixed,rho_max_tfixed
    IF(ALLOCATED(time_tfixed)) DEALLOCATE(time_tfixed)
    IF(ALLOCATED(coef_tfixed)) DEALLOCATE(coef_tfixed)
    ALLOCATE(time_tfixed(ntime_tfixed_max))
    ALLOCATE(coef_tfixed(0:ndata_tfixed_max,ntime_tfixed_max))
    DO ntime=1,ntime_tfixed_max
       READ(NFL,*) time_tfixed(ntime)
       READ(NFL,*) (coef_tfixed(ndata,ntime),ndata=0,ndata_tfixed_max)
    END DO
    CLOSE(NFL)
    RETURN
  END SUBROUTINE tr_prep_tfixed

  ! *** read Ar density profile from CSV file ***
  ! CSV format: r/a,ne,nD,nT,nA_thermal,nA_fast,nAr
  ! nAr is in units of 10^20 m^-3

  SUBROUTINE tr_prep_arfixed
    USE trcomm
    USE libfio
    IMPLICIT NONE
    INTEGER:: nfl,nr,ndata_csv,i,j,ierr,ios
    INTEGER,PARAMETER:: NMAX_CSV=300
    REAL(rkind):: rho_csv(NMAX_CSV),nar_csv(NMAX_CSV)
    REAL(rkind):: dummy(5),rho,frac
    CHARACTER(LEN=256):: line

    IF(model_arfixed.EQ.0) RETURN

    NFL=13
    CALL fropen(NFL,knam_arfixed,1,0,'ar',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,'(A)') 'XX tr_prep_arfixed: cannot open file '//TRIM(knam_arfixed)
       RETURN
    END IF

    ! Skip header line
    READ(NFL,'(A)',IOSTAT=ios) line

    ! Read data
    ndata_csv=0
    DO i=1,NMAX_CSV
       READ(NFL,*,IOSTAT=ios) rho_csv(i),dummy(1),dummy(2),dummy(3),dummy(4),dummy(5),nar_csv(i)
       IF(ios.NE.0) EXIT
       ndata_csv=i
    END DO
    CLOSE(NFL)

    IF(ndata_csv.LT.2) THEN
       WRITE(6,'(A)') 'XX tr_prep_arfixed: insufficient data in file'
       RETURN
    END IF

    ! nAr in CSV is already in 10^20/m^-3 (same as ne, nD, nT columns)
    ! No unit conversion needed

    WRITE(6,'(A,I5,A)') '## tr_prep_arfixed: read ',ndata_csv,' points from '//TRIM(knam_arfixed)

    ! Interpolate to TR grid and store in ANAR (simple linear interpolation)
    DO nr=1,NRMAX
       rho = rm(nr)
       ! Find interpolation interval
       IF(rho.LE.rho_csv(1)) THEN
          ANAR(nr) = nar_csv(1)
       ELSE IF(rho.GE.rho_csv(ndata_csv)) THEN
          ANAR(nr) = nar_csv(ndata_csv)
       ELSE
          DO j=1,ndata_csv-1
             IF(rho.GE.rho_csv(j).AND.rho.LT.rho_csv(j+1)) THEN
                frac = (rho - rho_csv(j))/(rho_csv(j+1) - rho_csv(j))
                ANAR(nr) = (1.D0-frac)*nar_csv(j) + frac*nar_csv(j+1)
                EXIT
             END IF
          END DO
       END IF
       ! nAr is already in 10^20 m^-3 units
    END DO

    WRITE(6,'(A,2ES12.4)') '## tr_prep_arfixed: ANAR range = ',MINVAL(ANAR(1:NRMAX)),MAXVAL(ANAR(1:NRMAX))

    RETURN
  END SUBROUTINE tr_prep_arfixed

! ============================================================
!  Read radiation profile (qrad) from external CSV file
!  File format: CSV with header line
!  Columns: r/a, qrad [MW/m^3]
! ============================================================
  SUBROUTINE tr_prep_radfixed
    USE trcomm
    USE libfio
    IMPLICIT NONE
    INTEGER:: nfl,nr,ndata_csv,i,j,ierr,ios
    INTEGER,PARAMETER:: NMAX_CSV=300
    REAL(rkind):: rho_csv(NMAX_CSV),qrad_csv(NMAX_CSV)
    REAL(rkind):: rho,frac
    CHARACTER(LEN=256):: line

    IF(model_radfixed.EQ.0) RETURN

    NFL=14
    CALL fropen(NFL,knam_radfixed,1,0,'qrad',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,'(A)') 'XX tr_prep_radfixed: cannot open file '//TRIM(knam_radfixed)
       model_radfixed=0
       RETURN
    END IF

    ! Skip header line
    READ(NFL,'(A)',IOSTAT=ios) line

    ! Read data: r/a, qrad [MW/m^3]
    ndata_csv=0
    DO i=1,NMAX_CSV
       READ(NFL,*,IOSTAT=ios) rho_csv(i),qrad_csv(i)
       IF(ios.NE.0) EXIT
       ndata_csv=i
    END DO
    CLOSE(NFL)

    IF(ndata_csv.LT.2) THEN
       WRITE(6,'(A)') 'XX tr_prep_radfixed: insufficient data in file'
       model_radfixed=0
       RETURN
    END IF

    WRITE(6,'(A,I5,A)') '## tr_prep_radfixed: read ',ndata_csv,' points from '//TRIM(knam_radfixed)

    ! Interpolate to TR grid and store in PRSUM_ext
    ! CSV is in MW/m^3, internal PRSUM is in W/m^3 → multiply by 1.D6
    DO nr=1,NRMAX
       rho = RG(nr)  ! r/a on TR grid
       IF(rho.LE.rho_csv(1)) THEN
          PRSUM_ext(nr) = qrad_csv(1) * 1.D6
       ELSE IF(rho.GE.rho_csv(ndata_csv)) THEN
          PRSUM_ext(nr) = qrad_csv(ndata_csv) * 1.D6
       ELSE
          DO j=1,ndata_csv-1
             IF(rho.GE.rho_csv(j).AND.rho.LT.rho_csv(j+1)) THEN
                frac = (rho - rho_csv(j))/(rho_csv(j+1) - rho_csv(j))
                PRSUM_ext(nr) = ((1.D0-frac)*qrad_csv(j) + frac*qrad_csv(j+1)) * 1.D6
                EXIT
             END IF
          END DO
       END IF
    END DO

    WRITE(6,'(A,2ES12.4)') '## tr_prep_radfixed: PRSUM_ext range = ', &
         MINVAL(PRSUM_ext(1:NRMAX)),MAXVAL(PRSUM_ext(1:NRMAX))

    RETURN
  END SUBROUTINE tr_prep_radfixed

! ============================================================
!  Read PRL profile (line radiation) from external CSV file
!  File format: CSV with header line
!  Columns: r/a, prl [MW/m^3]
! ============================================================
  SUBROUTINE tr_prep_prlfixed
    USE trcomm
    USE libfio
    IMPLICIT NONE
    INTEGER:: nfl,nr,ndata_csv,i,j,ierr,ios
    INTEGER,PARAMETER:: NMAX_CSV=300
    REAL(rkind):: rho_csv(NMAX_CSV),prl_csv(NMAX_CSV)
    REAL(rkind):: rho,frac
    CHARACTER(LEN=256):: line

    IF(model_prlfixed.EQ.0) RETURN

    NFL=16
    CALL fropen(NFL,knam_prlfixed,1,0,'prl',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,'(A)') 'XX tr_prep_prlfixed: cannot open file '//TRIM(knam_prlfixed)
       model_prlfixed=0
       RETURN
    END IF

    ! Skip header line
    READ(NFL,'(A)',IOSTAT=ios) line

    ! Read data: r/a, prl [MW/m^3]
    ndata_csv=0
    DO i=1,NMAX_CSV
       READ(NFL,*,IOSTAT=ios) rho_csv(i),prl_csv(i)
       IF(ios.NE.0) EXIT
       ndata_csv=i
    END DO
    CLOSE(NFL)

    IF(ndata_csv.LT.2) THEN
       WRITE(6,'(A)') 'XX tr_prep_prlfixed: insufficient data in file'
       model_prlfixed=0
       RETURN
    END IF

    WRITE(6,'(A,I5,A)') '## tr_prep_prlfixed: read ',ndata_csv,' points from '//TRIM(knam_prlfixed)

    ! Interpolate to TR grid and store in PRL_ext
    ! CSV is in MW/m^3, internal PRL is in W/m^3 -> multiply by 1.D6
    DO nr=1,NRMAX
       rho = RG(nr)  ! r/a on TR grid
       IF(rho.LE.rho_csv(1)) THEN
          PRL_ext(nr) = prl_csv(1) * 1.D6
       ELSE IF(rho.GE.rho_csv(ndata_csv)) THEN
          PRL_ext(nr) = prl_csv(ndata_csv) * 1.D6
       ELSE
          DO j=1,ndata_csv-1
             IF(rho.GE.rho_csv(j).AND.rho.LT.rho_csv(j+1)) THEN
                frac = (rho - rho_csv(j))/(rho_csv(j+1) - rho_csv(j))
                PRL_ext(nr) = ((1.D0-frac)*prl_csv(j) + frac*prl_csv(j+1)) * 1.D6
                EXIT
             END IF
          END DO
       END IF
    END DO

    WRITE(6,'(A,2ES12.4)') '## tr_prep_prlfixed: PRL_ext range = ', &
         MINVAL(PRL_ext(1:NRMAX)),MAXVAL(PRL_ext(1:NRMAX))

    RETURN
  END SUBROUTINE tr_prep_prlfixed

! ============================================================
!  Read PRF profile (qrfe) from external CSV file
!  File format: CSV with header line
!  Columns: r/a, qrfe [MW/m^3]
! ============================================================
  SUBROUTINE tr_prep_prffixed
    USE trcomm
    USE libfio
    IMPLICIT NONE
    INTEGER:: nfl,nr,ndata_csv,i,j,ierr,ios
    INTEGER,PARAMETER:: NMAX_CSV=300
    REAL(rkind):: rho_csv(NMAX_CSV),qrfe_csv(NMAX_CSV)
    REAL(rkind):: rho,frac
    CHARACTER(LEN=256):: line

    IF(model_prffixed.EQ.0) RETURN

    NFL=15
    CALL fropen(NFL,knam_prffixed,1,0,'qrfe',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,'(A)') 'XX tr_prep_prffixed: cannot open file '//TRIM(knam_prffixed)
       model_prffixed=0
       RETURN
    END IF

    ! Skip header line
    READ(NFL,'(A)',IOSTAT=ios) line

    ! Read data: r/a, qrfe [MW/m^3]
    ndata_csv=0
    DO i=1,NMAX_CSV
       READ(NFL,*,IOSTAT=ios) rho_csv(i),qrfe_csv(i)
       IF(ios.NE.0) EXIT
       ndata_csv=i
    END DO
    CLOSE(NFL)

    IF(ndata_csv.LT.2) THEN
       WRITE(6,'(A)') 'XX tr_prep_prffixed: insufficient data in file'
       model_prffixed=0
       RETURN
    END IF

    WRITE(6,'(A,I5,A)') '## tr_prep_prffixed: read ',ndata_csv,' points from '//TRIM(knam_prffixed)

    ! Interpolate to TR grid and store in PRF_ext
    ! CSV is in MW/m^3, internal PRF is in W/m^3 → multiply by 1.D6
    DO nr=1,NRMAX
       rho = RG(nr)  ! r/a on TR grid
       IF(rho.LE.rho_csv(1)) THEN
          PRF_ext(nr) = qrfe_csv(1) * 1.D6
       ELSE IF(rho.GE.rho_csv(ndata_csv)) THEN
          PRF_ext(nr) = qrfe_csv(ndata_csv) * 1.D6
       ELSE
          DO j=1,ndata_csv-1
             IF(rho.GE.rho_csv(j).AND.rho.LT.rho_csv(j+1)) THEN
                frac = (rho - rho_csv(j))/(rho_csv(j+1) - rho_csv(j))
                PRF_ext(nr) = ((1.D0-frac)*qrfe_csv(j) + frac*qrfe_csv(j+1)) * 1.D6
                EXIT
             END IF
          END DO
       END IF
    END DO

    WRITE(6,'(A,2ES12.4)') '## tr_prep_prffixed: PRF_ext range = ', &
         MINVAL(PRF_ext(1:NRMAX)),MAXVAL(PRF_ext(1:NRMAX))

    RETURN
  END SUBROUTINE tr_prep_prffixed

  ! *** read chi (thermal diffusivity) from external file ***
  ! File format: space-separated with header line
  ! Columns: r[m]  ne  ni  chi_e  chi_i  Se  Si  ...
  ! chi_e, chi_i are in m^2/s
  ! Note: r is in meters, will be converted to r/a using RA

  SUBROUTINE tr_prep_chifixed
    USE trcomm
    USE libfio
    IMPLICIT NONE
    INTEGER:: nfl,nr,ndata_csv,i,j,ierr,ios
    INTEGER,PARAMETER:: NMAX_CSV=300
    REAL(rkind):: r_csv(NMAX_CSV),rho_csv(NMAX_CSV),chi_e_csv(NMAX_CSV),chi_i_csv(NMAX_CSV)
    REAL(rkind):: dummy_ne,dummy_ni,rho,frac
    CHARACTER(LEN=512):: line

    ! Allow reading chi file when model_chifixed>=1 OR MDLKAI=170:179
    IF(model_chifixed.EQ.0 .AND. (MDLKAI.LT.170 .OR. MDLKAI.GT.179)) RETURN

    NFL=14
    CALL fropen(NFL,knam_chifixed,1,0,'chi',ierr)
    IF(ierr.NE.0) THEN
       WRITE(6,'(A)') 'XX tr_prep_chifixed: cannot open file '//TRIM(knam_chifixed)
       RETURN
    END IF

    ! Skip header line
    READ(NFL,'(A)',IOSTAT=ios) line

    ! Read data: r[m]  ne  ni  chi_e  chi_i  ...
    ndata_csv=0
    DO i=1,NMAX_CSV
       READ(NFL,*,IOSTAT=ios) r_csv(i),dummy_ne,dummy_ni,chi_e_csv(i),chi_i_csv(i)
       IF(ios.NE.0) EXIT
       ndata_csv=i
    END DO
    CLOSE(NFL)

    IF(ndata_csv.LT.2) THEN
       WRITE(6,'(A)') 'XX tr_prep_chifixed: insufficient data in file'
       RETURN
    END IF

    ! Convert r[m] to r/a using RA (minor radius)
    DO i=1,ndata_csv
       rho_csv(i) = r_csv(i) 
    END DO

    WRITE(6,'(A,I5,A)') '## tr_prep_chifixed: read ',ndata_csv,' points from '//TRIM(knam_chifixed)
    WRITE(6,'(A,F8.4,A,F8.4)') '## tr_prep_chifixed: r/a range = ',rho_csv(1),' to ',rho_csv(ndata_csv)

    ! Interpolate to TR grid (simple linear interpolation)
    DO nr=1,NRMAX
       rho = rm(nr)
       ! Find interpolation interval
       IF(rho.LE.rho_csv(1)) THEN
          AKEXT_E(nr) = chi_e_csv(1)
          AKEXT_I(nr) = chi_i_csv(1)
       ELSE IF(rho.GE.rho_csv(ndata_csv)) THEN
          AKEXT_E(nr) = chi_e_csv(ndata_csv)
          AKEXT_I(nr) = chi_i_csv(ndata_csv)
       ELSE
          DO j=1,ndata_csv-1
             IF(rho.GE.rho_csv(j).AND.rho.LT.rho_csv(j+1)) THEN
                frac = (rho - rho_csv(j))/(rho_csv(j+1) - rho_csv(j))
                AKEXT_E(nr) = (1.D0-frac)*chi_e_csv(j) + frac*chi_e_csv(j+1)
                AKEXT_I(nr) = (1.D0-frac)*chi_i_csv(j) + frac*chi_i_csv(j+1)
                EXIT
             END IF
          END DO
       END IF

    END DO

    ! Apply chi modification based on model_chifixed value
    IF(model_chifixed.EQ.1) THEN
       ! model_chifixed=1: Reduce core chi by 0.5x (r/a < 0.4)
       DO nr=1,NRMAX
          rho = rm(nr)
       END DO
       WRITE(6,'(A)') '## tr_prep_chifixed: mode=1, core chi reduced to 0.5x for r/a < 0.4'

    ELSE IF(model_chifixed.EQ.2) THEN
       ! model_chifixed=2: Multiply core chi by chifixed_factor (r/a < 0.4)
       DO nr=1,NRMAX
          rho = rm(nr)
          frac = 0.5D0 * (1.D0 - TANH((rho - 0.4D0) / 0.1D0))  ! 1 at center, 0 at edge
          AKEXT_E(nr) = AKEXT_E(nr) * (1.D0 + (chifixed_factor - 1.D0) * frac)
          AKEXT_I(nr) = AKEXT_I(nr) * (1.D0 + (chifixed_factor - 1.D0) * frac)
       END DO
       WRITE(6,'(A,F8.3,A)') '## tr_prep_chifixed: mode=2, core chi multiplied by ',chifixed_factor,' for r/a < 0.4'

    ELSE IF(model_chifixed.EQ.3) THEN
       ! model_chifixed=3: Multiply edge chi by chifixed_factor (r/a > 0.5)
       DO nr=1,NRMAX
          rho = rm(nr)
          frac = 0.5D0 * (1.D0 + TANH((rho - 0.5D0) / 0.1D0))  ! 0 at center, 1 at edge
          AKEXT_E(nr) = AKEXT_E(nr) * (1.D0 + (chifixed_factor - 1.D0) * frac)
          AKEXT_I(nr) = AKEXT_I(nr) * (1.D0 + (chifixed_factor - 1.D0) * frac)
       END DO
       WRITE(6,'(A,F8.3,A)') '## tr_prep_chifixed: mode=3, edge chi multiplied by ',chifixed_factor,' for r/a > 0.5'
    END IF

    WRITE(6,'(A,2ES12.4)') '## tr_prep_chifixed: chi_e range = ',MINVAL(AKEXT_E(1:NRMAX)),MAXVAL(AKEXT_E(1:NRMAX))
    WRITE(6,'(A,2ES12.4)') '## tr_prep_chifixed: chi_i range = ',MINVAL(AKEXT_I(1:NRMAX)),MAXVAL(AKEXT_I(1:NRMAX))

    RETURN
  END SUBROUTINE tr_prep_chifixed

END MODULE trfixed
