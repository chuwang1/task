! MODULE tiprep

MODULE tiprep

CONTAINS

  SUBROUTINE ti_load_omfit_csv(ierr)
    USE ticomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nid,ios,nr,nsa,ns,nzini,ncol
    REAL(rkind):: rho,ne_in,nar_in,nhe_in,te_in
    REAL(rkind):: skip_cols(10)
    CHARACTER(LEN=256):: line,header

    ierr=0
    IF(LEN_TRIM(omfit_profile_csv).EQ.0) RETURN

    nid=29
    OPEN(nid,FILE=TRIM(omfit_profile_csv),STATUS='OLD',FORM='FORMATTED',IOSTAT=ios)
    IF(ios.NE.0) THEN
       WRITE(6,'(A,1X,A)') 'XX ti_load_omfit_csv: cannot open',TRIM(omfit_profile_csv)
       ierr=1
       RETURN
    END IF

    READ(nid,'(A)',IOSTAT=ios) header
    IF(ios.NE.0) THEN
       WRITE(6,'(A)') 'XX ti_load_omfit_csv: cannot read header'
       CLOSE(nid)
       ierr=2
       RETURN
    END IF

!   --- count columns to skip before rho ---
!       header format: ...,rho,ne,...
!       find position of 'rho' in header to determine skip count
    ncol=0
    CALL ti_csv_find_rho_col(header,ncol)
    WRITE(6,'(A,I3)') '## ti_load_omfit_csv: columns to skip before rho =',ncol

    DO nr=1,nrmax
       READ(nid,'(A)',IOSTAT=ios) line
       IF(ios.NE.0) THEN
          WRITE(6,'(A,I5)') 'XX ti_load_omfit_csv: insufficient rows at NR=',nr
          CLOSE(nid)
          ierr=3
          RETURN
       END IF
       IF(ncol.GT.0) THEN
          READ(line,*,IOSTAT=ios) skip_cols(1:ncol),rho,ne_in,nar_in,nhe_in,te_in
       ELSE
          READ(line,*,IOSTAT=ios) rho,ne_in,nar_in,nhe_in,te_in
       END IF
       IF(ios.NE.0) THEN
          WRITE(6,'(A,I5)') 'XX ti_load_omfit_csv: parse error at NR=',nr
          CLOSE(nid)
          ierr=4
          RETURN
       END IF

       DO nsa=1,nsa_max
          ns=ns_nsa(nsa)
          SELECT CASE(ns)
          CASE(1)
             RNA(nsa,nr)=ne_in
             RTA(nsa,nr)=te_in
             RUA(nsa,nr)=0.D0
          CASE(2)
             RNA(nsa,nr)=ne_in
             RTA(nsa,nr)=te_in
             RUA(nsa,nr)=0.D0
          CASE(3)
             RNA(nsa,nr)=0.D0
             RTA(nsa,nr)=te_in
             RUA(nsa,nr)=0.D0
             nzini=NZINI_NS(ns)
             IF(nzini.GE.0 .AND. NINT(PZA(nsa)).EQ.nzini) RNA(nsa,nr)=nar_in
          END SELECT
       END DO
    END DO

    READ(nid,'(A)',IOSTAT=ios) line
    IF(ios.EQ.0) THEN
       WRITE(6,'(A)') 'XX ti_load_omfit_csv: row count exceeds NRMAX'
       CLOSE(nid)
       ierr=5
       RETURN
    END IF

    CLOSE(nid)
  END SUBROUTINE ti_load_omfit_csv

  SUBROUTINE ti_csv_find_rho_col(header,ncol)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN):: header
    INTEGER,INTENT(OUT):: ncol
    INTEGER:: i,istart,icol
    CHARACTER(LEN=64):: token

    ncol=0
    icol=0
    istart=1
    DO i=1,LEN_TRIM(header)+1
       IF(i.GT.LEN_TRIM(header) .OR. header(i:i).EQ.',') THEN
          token=ADJUSTL(header(istart:i-1))
          IF(TRIM(token).EQ.'rho') THEN
             ncol=icol
             RETURN
          END IF
          icol=icol+1
          istart=i+1
       END IF
    END DO
    ncol=0
  END SUBROUTINE ti_csv_find_rho_col

  SUBROUTINE ti_prep(IERR)

    USE ticomm
    USE plprof
    USE tirecord
    USE ADPOST
    USE ADF11
    USE libmtx
    USE libmpi
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER,SAVE:: init_adpost=0
    INTEGER,SAVE:: init_adas=0
    INTEGER:: NS,NSA,NZ,NEQ,NR,NV,ID_adpost,ID_adas,l1,l2
    REAL(rkind):: RHON,DRX
    TYPE(pl_prf_type),DIMENSION(:),ALLOCATABLE:: PLF

    IERR=0

!    ID_adpost=0
    ID_adpost=1
    ID_adas=0
    DO NS=1,NSMAX
       SELECT CASE(ID_NS(NS))
       CASE(5,6)
          ID_adpost=1
       CASE(10,11,12)
          ID_adas=1
       END SELECT
    END DO

    IF(ID_adpost.EQ.1) THEN
       IF(init_adpost.EQ.0) THEN
          IF(nrank.eq.0) THEN
             l1=LEN_TRIM(adpost_dir)
             L2=LEN_TRIM(adpost_filename)
             CALL read_adpost(adpost_dir(1:l1)//adpost_filename(1:l2),IERR)
             IF(IERR.NE.0) WRITE(6,*) 'XX read_adpost: IERR=',IERR
          END IF
          IF(IERR.NE.0) RETURN
          CALL broadcast_adpost(ierr)
          IF(ierr.NE.0) WRITE(6,*) 'XX broadcast_adpost: IERR=',IERR
          IF(IERR.NE.0) RETURN
          init_adpost=1
       END IF
    END IF
    IF(ID_adas.EQ.1) THEN
       IF(init_adas.EQ.0) THEN
          IF(nrank.EQ.0) THEN
             l1=LEN_TRIM(adas_adf11_dir)
             L2=LEN_TRIM(adas_adf11_filename)
             CALL LOAD_ADF11_bin( &
                  adas_adf11_dir(1:l1)//adas_adf11_filename(1:l2),IERR)
             IF(IERR.NE.0) WRITE(6,*) 'XX load_ADF11_bin: IERR=',IERR
          END IF
          IF(IERR.NE.0) RETURN
          CALL broadcast_ADF11_bin(IERR)
          IF(IERR.NE.0) WRITE(6,*) 'XX broadcast_ADF11_bin: IERR=',IERR
          IF(IERR.NE.0) RETURN
          init_adas=1

          CALL CALC_ADF11(1,1,0.D0,0.D0,DRX,IERR)
          IF(IERR.NE.0) WRITE(6,*) &
               'XX tiprep: CALC_ADF11 error: IERR,nrank=',IERR,NRANK
          IF(IERR.NE.0) RETURN
       END IF
    END IF

    IF(IERR.NE.0) RETURN

!   *** count NSA_MAX: number of active particles ***

    NSA=0
    DO NS=1,NSMAX
       SELECT CASE(ID_NS(NS))
       CASE(0)
          CONTINUE
       CASE(-1,1,2,5,6)
          NSA=NSA+1
       CASE(10,11,12)
          DO NZ=NZMIN_NS(NS),NZMAX_NS(NS)
             NSA=NSA+1
          END DO
       CASE DEFAULT
          WRITE(6,'(A,2I5)') &
               'XX ti_prep: Undefined ID_NS(NS): ID_NS,NS:', &
               ID_NS(NS),NS
          STOP
       END SELECT
    END DO
    nsa_max=NSA

!   *** ALLOCATE array for nsa_max and NRMAX ***

    CALL allocate_ticomm(IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX tiprep: allocate_ticomm ERROR: IERR=',IERR
       STOP
    END IF

!   *** Define NSA variables ***

    NSA=0
    DO NS=1,NSMAX
       SELECT CASE(ID_NS(NS))
       CASE(0)
          ND_adpost(NS)=0
       CASE(-1,1,2)
          NSA=NSA+1
          PMA(NSA)=PM(NS)
          PZA(NSA)=PZ(NS)
          PZ2A(NSA)=PZ(NS)**2
          NS_NSA(NSA)=NS
          NSA_DN(NSA)=0
          NSA_UP(NSA)=0
          ND_adpost(NS)=0
       CASE(5,6)
          NSA=NSA+1
          PMA(NSA)=PM(NS)
          PZA(NSA)=PZ(NS)      ! for ID_NS=5,6, PZA will be ealuated later
          PZ2A(NSA)=PZ(NS)**2  ! for ID_NS=5,6, PZ2A will be ealuated later
          NS_NSA(NSA)=NS
          NSA_DN(NSA)=0
          NSA_UP(NSA)=0
          ND_adpost(NS)=ND_NPA_ADPOST(NPA(NS))
       CASE(10,11,12)
          DO NZ=NZMIN_NS(NS),NZMAX_NS(NS)
             NSA=NSA+1
             PMA(NSA)=PM(NS)
             PZA(NSA)=DBLE(NZ)
             PZ2A(NSA)=DBLE(NZ)**2
             NS_NSA(NSA)=NS
             IF(NZ.EQ.NZMIN_NS(NS)) THEN
                NSA_DN(NSA)=0
             ELSE
                NSA_DN(NSA)=NSA-1
             END IF
             IF(NZ.EQ.NZMAX_NS(NS)) THEN
                NSA_UP(NSA)=0
             ELSE
                NSA_UP(NSA)=NSA+1
             END IF
          END DO
          ND_adpost(NS)=0
       END SELECT
    END DO
    IF(NSA.GT.nsa_max) THEN
       WRITE(6,'(A)') 'XX ti_prep: INCONSISTENT nsa_max'
       STOP
    ELSE IF(NSA.LT.nsa_max) THEN
       IF(nrank.EQ.0) THEN
          WRITE(6,'(A,I5)') &
               '!! calculated nsa_max is lower than assumed nsa_max:',nsa_max
          nsa_max=NSA
          WRITE(6,'(A,I5)') &
               '                              nsa_max is reduced to:',nsa_max
       END IF
    END IF

!   *** Count NEQMAX: Size of equation ***

    NEQ=0
    IF(MODEL_EQB.EQ.1) NEQ=NEQ+1
    DO NS=1,NSMAX
       SELECT CASE(ID_NS(NS))
       CASE(0)
          CONTINUE
       CASE(-1,1,2,5,6)
          IF(MODEL_EQN.EQ.1) NEQ=NEQ+1
          IF(MODEL_EQT.EQ.1) NEQ=NEQ+1
          IF(MODEL_EQU.EQ.1) NEQ=NEQ+1
       CASE(10)
          IF(MODEL_EQN.EQ.1) NEQ=NEQ+NZMAX_NS(NS)-NZMIN_NS(NS)+1
          IF(MODEL_EQT.EQ.1) NEQ=NEQ+1
          IF(MODEL_EQU.EQ.1) NEQ=NEQ+1
       CASE(11)
          IF(MODEL_EQN.EQ.1) NEQ=NEQ+NZMAX_NS(NS)-NZMIN_NS(NS)+1
          IF(MODEL_EQT.EQ.1) NEQ=NEQ+NZMAX_NS(NS)-NZMIN_NS(NS)+1
          IF(MODEL_EQU.EQ.1) NEQ=NEQ+1
       CASE(12)
          IF(MODEL_EQN.EQ.1) NEQ=NEQ+NZMAX_NS(NS)-NZMIN_NS(NS)+1
          IF(MODEL_EQT.EQ.1) NEQ=NEQ+NZMAX_NS(NS)-NZMIN_NS(NS)+1
          IF(MODEL_EQU.EQ.1) NEQ=NEQ+NZMAX_NS(NS)-NZMIN_NS(NS)+1
       END SELECT
    END DO
    NEQMAX=NEQ

!   *** ALLOCATE array for NEQMAX and NRMAX ***

    CALL allocate_neqmax(IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX ti_prep: allocate_neqmax ERROR: IERR=',IERR
       STOP
    END IF

!   *** Define NEQ variables ***

    NEQ=0
    IF(MODEL_EQB.EQ.1) THEN
       NEQ=NEQ+1
       NSA_NEQ(NEQ)=0
       NV_NEQ(NEQ)=1
    END IF
    DO NSA=1,nsa_max
       NS=NS_NSA(NSA)
       IF(MODEL_EQN.EQ.1) THEN
          SELECT CASE(ID_NS(NS))
          CASE(-1,1,2,5,6,10,11,12)
             NEQ=NEQ+1
             NSA_NEQ(NEQ)=NSA
             NV_NEQ(NEQ)=1
          END SELECT
       END IF
       IF(MODEL_EQT.EQ.1) THEN
          SELECT CASE(ID_NS(NS))
          CASE(-1,1,2,5,6,10,11)
             NEQ=NEQ+1
             NSA_NEQ(NEQ)=NSA
             NV_NEQ(NEQ)=2
          CASE(12)
             IF(NINT(PZ(NSA)).EQ.NZMIN_NS(NS)) THEN
                NEQ=NEQ+1
                NSA_NEQ(NEQ)=NSA
                NV_NEQ(NEQ)=2
             END IF
          END SELECT
       END IF
       IF(MODEL_EQU.EQ.1) THEN
          SELECT CASE(ID_NS(NS))
          CASE(-1,1,2,5,6,10)
             NEQ=NEQ+1
             NSA_NEQ(NEQ)=NSA
             NV_NEQ(NEQ)=3
          CASE(11,12)
             IF(NINT(PZ(NSA)).EQ.NZMIN_NS(NS)) THEN
                NEQ=NEQ+1
                NSA_NEQ(NEQ)=NSA
                NV_NEQ(NEQ)=3
             END IF
          END SELECT
       END IF
    END DO
    IF(NEQ.NE.NEQMAX) THEN
       WRITE(6,*) 'XX ti_prep: INCONSISTENT NEQMAX'
       STOP
    END IF

    NEQ_NVNSA(1:3,0:nsa_max)=0.D0
    DO NEQ=1,NEQMAX
       NV=NV_NEQ(NEQ)
       NSA=NSA_NEQ(NEQ)
       NEQ_NVNSA(NV,NSA)=NEQ
    END DO

!   *** Display NEQ variables ***

    IF(nrank.EQ.0) THEN
       WRITE(6,*)
       WRITE(6,'(A)') &
          '       NEQ   NS  NSA  NPA          PM        PZ   NV  ID_NS KID_NS'
       WRITE(6,'(5X,4I5,F12.5,F10.3,I5,I7,5X,A2)') &
               (NEQ,NS_NSA(NSA_NEQ(NEQ)),NSA_NEQ(NEQ), &
                NPA(NS_NSA(NSA_NEQ(NEQ))), &
                PMA(NSA_NEQ(NEQ)), &
                PZA(NSA_NEQ(NEQ)), &
                NV_NEQ(NEQ), &
                ID_NS(NS_NSA(NSA_NEQ(NEQ))), &
                KID_NS(NS_NSA(NSA_NEQ(NEQ))), &
                NEQ=1,NEQMAX)
    END IF

!   *** radial mesh ***

    DR=RA/NRMAX

    RG(0)=0.D0
    DO NR=1,NRMAX
       RM(NR)=(NR-0.5D0)*DR
       RG(NR)=NR*DR
    END DO

    voltot=0.D0
    DO NR=1,NRMAX
       DVRHO(NR)=4.D0*PI*PI*RR*RKAP*RM(NR)
       voltot=voltot+DVRHO(NR)*DR
    END DO
    DVRHOS=4.D0*PI*PI*RR*RKAP*RG(NRMAX)

!   *** initial profile ***

    ALLOCATE(PLF(NSMAX))
     DO NR=1,NRMAX
        RHON=RM(NR)/RA
        CALL pl_prof(RHON,plf)
       DO NSA=1,nsa_max
          NS=NS_NSA(NSA)
          IF(ID_NS(NS).LE.2) THEN
             RNA(NSA,NR)=PLF(NS)%RN
          ELSE
             RNA(NSA,NR)=0.D0
             IF(PZA(NSA).EQ.0.D0.AND.NR.EQ.NRMAX) RNA(NSA,NR)=PLF(NS)%RN
          END IF
          RTA(NSA,NR)=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
          RUA(NSA,NR)=PLF(NS)%RU
       END DO
!       WRITE(6,'(A,I5,1P3E12.4)') &
!            'NR,RHON,RNE,RTE=',NR,RHON,RNA(1,NR),RTA(1,NR)
     END DO
     DEALLOCATE(PLF)

     CALL ti_load_omfit_csv(ierr)
     IF(ierr.NE.0) RETURN

! --- setup boundary conditions ---

    DO NEQ=1,NEQMAX
       NV=NV_NEQ(NEQ)
       NSA=NSA_NEQ(NEQ)
       NS=NS_NSA(NSA)
       MODEL_BNDA(NV,NSA)=1    
       BND_VALUEA(NV,NSA)=0.D0
       IF(NINT(PZA(NSA)).EQ.NZMIN_NS(NS)) THEN
          MODEL_BNDA(NV,NSA)=MODEL_BND(NV,NS)
          BND_VALUEA(NV,NSA)=BND_VALUE(NV,NS)
       END IF
!       write(6,'(A,7I5,1P2E12.4)') &
!            'PREP:',NEQ,NV,NSA,NS,NINT(PZA(NSA)),NZMIN_NS(NS), &
!            MODEL_BND(NV,NS),BND_VALUE(NV,NS),BND_VALUEA(NV,NSA)
    END DO

! --- time initializaion ---

    NT=0
    T=0.D0
    ngt_max=0
    ngr_max=0
    IF(ALLOCATED(gt)) DEALLOCATE(gt)
    IF(ALLOCATED(gvt)) DEALLOCATE(gvt)
    IF(ALLOCATED(gvta)) DEALLOCATE(gvta)
    ngt_allocate_max=0
    IF(ALLOCATED(grt)) DEALLOCATE(grt)
    IF(ALLOCATED(gvrt)) DEALLOCATE(gvrt)
    IF(ALLOCATED(gvrta)) DEALLOCATE(gvrta)
    ngr_allocate_max=0

! --- motrix solver initializaion ---

    imax=NRMAX*NEQMAX
    jwidth=4*NEQMAX-1
    CALL mtx_setup(imax,istart,iend,jwidth)
    NR_START=(istart+NEQMAX-1)/NEQMAX
    NR_END=(iend+NEQMAX-1)/NEQMAX
    CALL mtx_cleanup

    WRITE(6,'(A,I5,6I8)') 'nrank: imax/s/e,nrmax/s/e=', &
                        nrank, imax,istart,iend,nrmax,nr_start,nr_end
    CALL mtx_barrier

! --- record initial values ---

    IF(MOD(NT,NTSTEP ).EQ.0) CALL ti_snap         ! integrate and save
    IF(MOD(NT,NGTSTEP).EQ.0) CALL ti_record_ngt   ! save for time history
    IF(MOD(NT,NGRSTEP).EQ.0) CALL ti_record_ngr   ! save for radial profile

    RETURN
  END SUBROUTINE ti_prep
END MODULE tiprep
