C     Pure EQ driver:
C       1) load gfile with MODELG=5
C       2) rebuild MDLEQF>=5 spline profiles inside EQ
C       3) solve GS and save eqdata
C
      PROGRAM MODELG5_REBUILD
      USE plinit,ONLY: pl_init
      USE plparm,ONLY: pl_parm
      INCLUDE '../eq/eqcomc.inc'
      INCLUDE '../eq/eqcom3.inc'
      INCLUDE '../eq/eqcom4.inc'
      CHARACTER GFILE*256,OUTEQ*256,RIPARG*64,LCFSARG*64,INITARG*64
      INTEGER NARG
      INTEGER IOSIP,IOSLC,IOSIN,IUSEARG,IUSELCFS_ARG,IINITPSI_ARG
      INTEGER IERR_FATAL
      REAL*8 RIPSET,RAXG,ZAXG
      LOGICAL LSETIP
C
      IERR=0
      IERR_FATAL=0
      LSETIP=.FALSE.
      IUSELCFS_ARG=-1
      IINITPSI_ARG=0
      GFILE='in/g260206.20000_teq_0114'
      OUTEQ='eqdata_modelg5_rebuild'
C
      NARG=IARGC()
      IF(NARG.GE.1) CALL GETARG(1,GFILE)
      IF(NARG.GE.2) CALL GETARG(2,OUTEQ)
      IF(NARG.GE.3) THEN
         CALL GETARG(3,RIPARG)
         IF(NARG.GE.4.AND.
     &      (TRIM(RIPARG).EQ.''.OR.TRIM(RIPARG).EQ.'-')) THEN
            CONTINUE
         ELSE
            READ(RIPARG,*,IOSTAT=IOSIP) RIPSET
            IF(IOSIP.EQ.0) THEN
               LSETIP=.TRUE.
            ELSE
               READ(RIPARG,*,IOSTAT=IOSLC) IUSEARG
               IF(NARG.EQ.3.AND.IOSLC.EQ.0.AND.
     &            (IUSEARG.EQ.0.OR.IUSEARG.EQ.1)) THEN
                  IUSELCFS_ARG=IUSEARG
               ELSE
                  WRITE(6,*) 'XX MODELG5_REBUILD: invalid arg3 = ',
     &                       TRIM(RIPARG)
                  WRITE(6,*)
     &             '   arg3 should be RIP(MA) or LCFS mode (0/1)'
                  STOP 1
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF(NARG.GE.4) THEN
         CALL GETARG(4,LCFSARG)
         READ(LCFSARG,*,IOSTAT=IOSLC) IUSEARG
         IF(IOSLC.EQ.0.AND.(IUSEARG.EQ.0.OR.IUSEARG.EQ.1)) THEN
            IUSELCFS_ARG=IUSEARG
         ELSE
            WRITE(6,*) 'XX MODELG5_REBUILD: invalid LCFS mode arg4 = ',
     &                 TRIM(LCFSARG)
            STOP 1
         ENDIF
      ENDIF
      IF(NARG.GE.5) THEN
         CALL GETARG(5,INITARG)
         READ(INITARG,*,IOSTAT=IOSIN) IUSEARG
         IF(IOSIN.EQ.0.AND.(IUSEARG.EQ.0.OR.IUSEARG.EQ.1)) THEN
            IINITPSI_ARG=IUSEARG
         ELSE
            WRITE(6,*)
     &      'XX MODELG5_REBUILD: invalid init-psi mode arg5 = ',
     &      TRIM(INITARG)
            STOP 1
         ENDIF
      ENDIF
C
      WRITE(6,'(A)')   '## PURE EQ MODELG=5 REBUILD'
      WRITE(6,'(A,A)') '## GFILE : ',TRIM(GFILE)
      WRITE(6,'(A,A)') '## OUTEQ : ',TRIM(OUTEQ)
      IF(LSETIP) WRITE(6,'(A,1PE12.4)') '## RIP target(MA): ',RIPSET
C
      CALL GSOPEN
      OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
C
      CALL PL_INIT
      CALL EQINIT
      CALL PL_PARM(1,'plparm',IERR)
      CALL EQPARM(1,'eqparm',IERR)
      CALL EQPARM(2,'NLPMAX=800',IERR)
      CALL EQPARM(2,'NTVMAX=1024',IERR)
      CALL EQPARM(2,'NTGMAX=64',IERR)
      CALL EQPARM(2,'NPRINT=2',IERR)
      IF(IUSELCFS_ARG.GE.0) IUSELCFS=IUSELCFS_ARG
C
      MODELG=5
      KNAMEQ=TRIM(GFILE)
      CALL EQLOAD(5,KNAMEQ,IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX MODELG5_REBUILD: EQLOAD: IERR=',IERR
         IERR_FATAL=1
         GOTO 9000
      ENDIF
C
      WRITE(6,'(A,1P2E12.4)') '## LOADED: PSIPA,PSITA=',PSIPA,PSITA
      WRITE(6,'(A,I3)') '## LCFS boundary mode (0=EQFBND,1=gfile): ',
     &                  IUSELCFS
      WRITE(6,'(A,I3)')
     &   '## Init PSI mode (0=analytic,1=from gfile PSIRZ): ',
     &   IINITPSI_ARG
      RAXG=RAXIS
      ZAXG=ZAXIS
      CALL LOG_AXIS_DELTA('after EQLOAD',RAXG,ZAXG)
      MDLEQF=6
      IF(LSETIP.AND.RIPSET.GT.0.D0) THEN
C        --- MDLEQF=6 with Ip matching ---
         RIP=RIPSET
         WRITE(6,'(A,1PE12.4)')
     &      '## MDLEQF=6: P,F + Ip match, RIP(MA)=',RIP
      ELSE
C        --- MDLEQF=6 with TJ=1 (use F directly, no Ip matching) ---
C        Keep gfile RIP for EQPSIN initial psi estimate, then set to 0
         WRITE(6,'(A)')
     &      '## MDLEQF=6: P,F direct (TJ=1, no Ip match)'
      ENDIF
      WRITE(6,'(A)') '## STEP: EQMESH'
      CALL EQMESH
      WRITE(6,'(A)') '## STEP: EQPSIN'
      CALL EQPSIN
      IF(.NOT.(LSETIP.AND.RIPSET.GT.0.D0)) THEN
C        Set RIP=0 after EQPSIN to trigger TJ=1 in EQRHSV
         RIP=0.D0
      ENDIF
      CALL LOG_AXIS_DELTA('after EQPSIN',RAXG,ZAXG)
      WRITE(6,'(A)') '## STEP: EQDEFB'
      CALL EQDEFB
      CALL LOG_AXIS_DELTA('after EQDEFB',RAXG,ZAXG)
      IF(IINITPSI_ARG.EQ.1) THEN
         WRITE(6,'(A)') '## STEP: EQPSI_FROM_PSIRZ'
         CALL EQPSI_FROM_PSIRZ(IERR)
         IF(IERR.NE.0) THEN
            WRITE(6,*) 'XX MODELG5_REBUILD: EQPSI_FROM_PSIRZ: IERR=',IERR
            IERR_FATAL=1
            GOTO 9000
         ENDIF
      ENDIF
C
      WRITE(6,'(A)') '## STEP: EQSET_MODELG5_TREQ'
      CALL EQSET_MODELG5_TREQ(IERR)
      CALL LOG_AXIS_DELTA('after EQSET_TREQ',RAXG,ZAXG)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX MODELG5_REBUILD: EQSET_MODELG5_TREQ: IERR=',IERR
         IERR_FATAL=1
         GOTO 9000
      ENDIF
C
      WRITE(6,'(A)') '## STEP: EQLOOP'
      CALL EQLOOP(IERR)
      CALL LOG_AXIS_DELTA('after EQLOOP',RAXG,ZAXG)
      IF(IERR.NE.0.AND.IERR.NE.100) THEN
         WRITE(6,*) 'XX MODELG5_REBUILD: EQLOOP: IERR=',IERR
         IERR_FATAL=1
         GOTO 9000
      ENDIF
C
      WRITE(6,'(A)') '## STEP: EQTORZ'
      CALL EQTORZ
      WRITE(6,'(A)') '## STEP: EQCALP'
      CALL EQCALP
      WRITE(6,'(A,1P4E12.4)') '## Pre-EQCALQ: RAXIS,ZAXIS,PSI0=',
     &     RAXIS,ZAXIS,PSI0
      WRITE(6,'(A,1P4E12.4)') '## Pre-EQCALQ: PSIPA,PSITA,RIPX=',
     &     PSIPA,PSITA,RIPX
      WRITE(6,'(A)') '## STEP: EQCALQ'
      CALL EQCALQ(IERR)
      CALL LOG_AXIS_DELTA('after EQCALQ',RAXG,ZAXG)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX MODELG5_REBUILD: EQCALQ: IERR=',IERR
         WRITE(6,*) '!! Skipping EQCALQ error, saving anyway'
         IERR_FATAL=0
      ENDIF
C
      KNAMEQ=TRIM(OUTEQ)
      CALL EQSAVE
      WRITE(6,'(A,1P3E12.4)') '## DONE: RIPX(MA),PSIPA,PSITA=',
     &                        RIPX,PSIPA,PSITA
C
 9000 CONTINUE
      CLOSE(7)
      CALL GSCLOS
      IF(IERR_FATAL.NE.0) THEN
         STOP 1
      ENDIF
      STOP
      END

      SUBROUTINE LOG_AXIS_DELTA(STAGE,RREF,ZREF)
      INCLUDE '../eq/eqcomc.inc'
      CHARACTER*(*) STAGE
      REAL*8 RREF,ZREF,DR,DZ,DD
      DR=RAXIS-RREF
      DZ=ZAXIS-ZREF
      DD=SQRT(DR*DR+DZ*DZ)
      WRITE(6,'(A,A)') '## AXIS ',STAGE
      WRITE(6,'(A,1P2E12.4)') '   current (R,Z)=',RAXIS,ZAXIS
      WRITE(6,'(A,1P2E12.4)') '   gfile   (R,Z)=',RREF,ZREF
      WRITE(6,'(A,1P3E12.4)') '   dR,dZ,|d|    =',DR,DZ,DD
      RETURN
      END
C
C
      SUBROUTINE EQSET_MODELG5_TREQ(IERR)
      USE libspl1d
      INCLUDE '../eq/eqcomq.inc'
      INCLUDE '../eq/eqcom4.inc'
      DIMENSION WORKP(NTRM+2),WORKF(NTRM+2),WORKQ(NTRM+2)
      DIMENSION PSITMP(NTRM+2),QTMP(NTRM+2)
      DIMENSION DERIV(NTRM+2)
      DIMENSION PSING(NPSM),DERIVN(NPSM)
      DIMENSION UPPN(4,NPSM),UFTN(4,NPSM),UQQN(4,NPSM)
      DIMENSION DERIVV(NRVM)
      DIMENSION XPNTR(NTRM+2),UQQR(4,NTRM+2)
      REAL*8 QFLOOR,PSIMAXG,QMINR,QMAXR,SCL,PSIWSCL
      REAL*8 PSIQCUT,QREF,QMINRAW,QMAXRAW
      LOGICAL IS_FINITE_D
      INTEGER NCUT0
C
      IERR=0
      QFLOOR=1.D-8
      PSIQCUT=0.86D0
      IF(NPSMAX.LT.3) THEN
         IERR=9301
         RETURN
      ENDIF
      IF(NRVMAX.LT.3) THEN
         IERR=9302
         RETURN
      ENDIF
      PSIMAXG=PSIPS(NPSMAX)
      PSIWSCL=PSIMAXG
      IF(PSIMAXG.LE.0.D0) THEN
         IERR=9303
         RETURN
      ENDIF
C
C     --- Validate gfile profiles on PSIPS grid ---
      DO NPS=2,NPSMAX
         IF(PSIPS(NPS).LE.PSIPS(NPS-1)) THEN
            IERR=9304
            WRITE(6,*) 'XX EQSET_MODELG5_TREQ: PSIPS not monotonic'
            RETURN
         ENDIF
      ENDDO
      DO NPS=1,NPSMAX
         IF(.NOT.IS_FINITE_D(PPPS(NPS)).OR.
     &      .NOT.IS_FINITE_D(TTPS(NPS)).OR.
     &      .NOT.IS_FINITE_D(QQPS(NPS))) THEN
            IERR=9305
            WRITE(6,*) 'XX EQSET_MODELG5_TREQ: invalid PPPS/TTPS/QQPS'
            RETURN
         ENDIF
      ENDDO
C
      NTRMAX=MIN(NPSMAX,NTRM)
      IF(NTRMAX.LT.2) THEN
         IERR=9306
         RETURN
      ENDIF
C
C     --- Build normalized psi_p splines for P,F,q from gfile ---
      DO NPS=1,NPSMAX
         PSING(NPS)=PSIPS(NPS)/PSIMAXG
      ENDDO
      CALL SPL1D(PSING,PPPS,DERIVN,UPPN,NPSMAX,0,IERR)
      IF(IERR.NE.0) THEN
         IERR=9310
         RETURN
      ENDIF
      CALL SPL1D(PSING,TTPS,DERIVN,UFTN,NPSMAX,0,IERR)
      IF(IERR.NE.0) THEN
         IERR=9311
         RETURN
      ENDIF
      CALL SPL1D(PSING,QQPS,DERIVN,UQQN,NPSMAX,0,IERR)
      IF(IERR.NE.0) THEN
         IERR=9312
         RETURN
      ENDIF
C
C     --- Sample P,F,q on uniform psi_pn grid ---
      QMINRAW=1.D30
      QMAXRAW=-1.D30
      DO NTR=1,NTRMAX
         PSIPNL=DBLE(NTR-1)/DBLE(NTRMAX-1)
         XPNTR(NTR)=PSIPNL
         CALL SPL1DF(PSIPNL,PPSL,PSING,UPPN,NPSMAX,IERR)
         IF(IERR.NE.0) RETURN
         CALL SPL1DF(PSIPNL,FPSL,PSING,UFTN,NPSMAX,IERR)
         IF(IERR.NE.0) RETURN
         CALL SPL1DF(PSIPNL,QPSL,PSING,UQQN,NPSMAX,IERR)
         IF(IERR.NE.0) RETURN
         IF(.NOT.IS_FINITE_D(QPSL).OR.QPSL.LE.QFLOOR) THEN
            IERR=9313
            WRITE(6,'(A,I5,1P2E12.4)') 'XX Q invalid at NTR=',
     &           NTR,PSIPNL,QPSL
            RETURN
         ENDIF
C
         WORKP(NTR)=PPSL*1.D-6
         WORKF(NTR)=FPSL
         WORKQ(NTR)=QPSL
         IF(QPSL.LT.QMINRAW) QMINRAW=QPSL
         IF(QPSL.GT.QMAXRAW) QMAXRAW=QPSL
      ENDDO
C
C     --- Hard q cut: keep q constant for psi_n >= PSIQCUT ---
      NCUT0=0
      DO NTR=1,NTRMAX
         IF(NCUT0.EQ.0.AND.XPNTR(NTR).GE.PSIQCUT) NCUT0=NTR
      ENDDO
      IF(NCUT0.GE.2.AND.NCUT0.LT.NTRMAX) THEN
         QREF=WORKQ(NCUT0)
         DO NTR=NCUT0+1,NTRMAX
            WORKQ(NTR)=QREF
         ENDDO
      ENDIF
C
      QMINR=1.D30
      QMAXR=-1.D30
      DO NTR=1,NTRMAX
         IF(.NOT.IS_FINITE_D(WORKQ(NTR)).OR.WORKQ(NTR).LE.QFLOOR) THEN
            IERR=9313
            WRITE(6,'(A,I5,1P2E12.4)')
     &           'XX regularized q invalid at NTR=',
     &           NTR,XPNTR(NTR),WORKQ(NTR)
            RETURN
         ENDIF
         QTMP(NTR)=WORKQ(NTR)
         IF(WORKQ(NTR).LT.QMINR) QMINR=WORKQ(NTR)
         IF(WORKQ(NTR).GT.QMAXR) QMAXR=WORKQ(NTR)
      ENDDO
C
      CALL SPL1D(XPNTR,WORKQ,DERIV,UQQR,NTRMAX,0,IERR)
      IF(IERR.NE.0) THEN
         IERR=9318
         RETURN
      ENDIF
C
C     --- Rebuild psi_t(psi_p) by integrating q(psi_p) ---
      PSITMP(1)=0.D0
      DO NTR=2,NTRMAX
         PSIPLM=DBLE(NTR-2)/DBLE(NTRMAX-1)*PSIWSCL
         PSIPLP=DBLE(NTR-1)/DBLE(NTRMAX-1)*PSIWSCL
         PSITMP(NTR)=PSITMP(NTR-1)
     &              +2.D0*QTMP(NTR)*QTMP(NTR-1)
     &              /(QTMP(NTR)+QTMP(NTR-1))
     &              *(PSIPLP-PSIPLM)
      ENDDO
      PSITA_NEW=PSITMP(NTRMAX)
      IF(PSITA_NEW.LE.0.D0.OR..NOT.IS_FINITE_D(PSITA_NEW)) THEN
         IERR=9314
         WRITE(6,*) 'XX EQSET_MODELG5_TREQ: invalid PSITA_NEW',PSITA_NEW
         RETURN
      ENDIF
C
      PSITRX(1)=0.D0
      DO NTR=2,NTRMAX
         PSITRX(NTR)=PSITMP(NTR)/PSITA_NEW
         IF(PSITRX(NTR).LE.PSITRX(NTR-1)+1.D-10)
     &      PSITRX(NTR)=PSITRX(NTR-1)+1.D-10
      ENDDO
      SCL=PSITRX(NTRMAX)
      IF(SCL.LE.0.D0) THEN
         IERR=9315
         RETURN
      ENDIF
      DO NTR=2,NTRMAX
         PSITRX(NTR)=PSITRX(NTR)/SCL
      ENDDO
      PSITRX(NTRMAX)=1.D0
C
      CALL SPL1D(PSITRX,WORKP,DERIV,UPPSI,NTRMAX,0,IERR)
      IF(IERR.NE.0) RETURN
      CALL SPL1D(PSITRX,WORKF,DERIV,UFPSI,NTRMAX,0,IERR)
      IF(IERR.NE.0) RETURN
      CALL SPL1D(PSITRX,WORKQ,DERIV,UQPSI,NTRMAX,0,IERR)
      IF(IERR.NE.0) RETURN
C
C     --- Rebuild global PSIPNV/PSITV/QPV/TTV maps used by EQPSITN/EQQPV ---
      QMINV=1.D30
      QMAXV=-1.D30
      DO NRV=1,NRVMAX
         PSIPNV(NRV)=DBLE(NRV-1)/DBLE(NRVMAX-1)
         PSIPV(NRV)=PSIPNV(NRV)*PSIPA
         CALL SPL1DF(PSIPNV(NRV),QPVL,XPNTR,UQQR,NTRMAX,IERR)
         IF(IERR.NE.0) RETURN
         CALL SPL1DF(PSIPNV(NRV),TTVL,PSING,UFTN,NPSMAX,IERR)
         IF(IERR.NE.0) RETURN
         IF(.NOT.IS_FINITE_D(QPVL).OR.QPVL.LE.QFLOOR) THEN
            IERR=9316
            WRITE(6,'(A,I5,1P2E12.4)') 'XX QPV invalid at NRV=',
     &           NRV,PSIPNV(NRV),QPVL
            RETURN
         ENDIF
         QPV(NRV)=QPVL
         TTV(NRV)=TTVL
         IF(QPVL.LT.QMINV) QMINV=QPVL
         IF(QPVL.GT.QMAXV) QMAXV=QPVL
      ENDDO
C
      PSITV(1)=0.D0
      DO NRV=2,NRVMAX
         PSITV(NRV)=PSITV(NRV-1)
     &             +2.D0*QPV(NRV)*QPV(NRV-1)
     &             /(QPV(NRV)+QPV(NRV-1))
     &             *(PSIWSCL/DBLE(NRVMAX-1))
      ENDDO
      PSITA=PSITV(NRVMAX)
      IF(PSITA.LE.0.D0.OR..NOT.IS_FINITE_D(PSITA)) THEN
         IERR=9317
         WRITE(6,*) 'XX EQSET_MODELG5_TREQ: invalid remapped PSITA=',PSITA
         RETURN
      ENDIF
C
      CALL SPL1D(PSIPNV,PSITV,DERIVV,UPSITV,NRVMAX,0,IERR)
      IF(IERR.NE.0) RETURN
      CALL SPL1D(PSIPNV,QPV,DERIVV,UQPV,NRVMAX,0,IERR)
      IF(IERR.NE.0) RETURN
      CALL SPL1D(PSIPNV,TTV,DERIVV,UTTV,NRVMAX,0,IERR)
      IF(IERR.NE.0) RETURN
C
      WRITE(6,'(A,1P2E12.4)') '## q range raw from gfile: ',
     &                        QMINRAW,QMAXRAW
      WRITE(6,'(A,1P2E12.4)') '## q range used for remap: ',
     &                        QMINR,QMAXR
      WRITE(6,'(A,1PE12.4)')  '## q edge cut (constant beyond psi_n): ',
     &                        PSIQCUT
      WRITE(6,'(A,1PE12.4)')  '## PSITA remapped from q(PSIPS): ',
     &                        PSITA
C
      RETURN
      END

      LOGICAL FUNCTION IS_FINITE_D(X)
      REAL*8 X
      IF(X.NE.X) THEN
         IS_FINITE_D=.FALSE.
      ELSEIF(ABS(X).GE.1.D300) THEN
         IS_FINITE_D=.FALSE.
      ELSE
         IS_FINITE_D=.TRUE.
      ENDIF
      RETURN
      END
