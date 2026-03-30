C     $Id$
C
C   ************************************** 
C   **       CALCULATE EQUILIBRIUM      **
C   ************************************** 
C
      SUBROUTINE EQCALC(IERR)
C      
      INCLUDE '../eq/eqcomc.inc'
C
      IERR=0
      CALL EQMESH
      IF(Q0*QA.GT.0.D0) THEN
         CALL EQPSIN
      ELSE
         CALL EQPSIR
      ENDIF
      CALL EQDEFB
      CALL EQLOOP(IERR)
         IF(IERR.NE.0.AND.IERR.NE.100) RETURN
      CALL EQTORZ
      CALL EQCALP
      RETURN
      END
C
C   ************************************************ 
C   **       Mesh Definition (sigma, theta)       **
C   ************************************************ 
C
      SUBROUTINE EQMESH
C      
      INCLUDE '../eq/eqcomc.inc'
C
      DSG=1.D0/NSGMAX
      DTG=2.D0*PI/NTGMAX 
C
      DO NSG=1,NSGMAX
         SIGM(NSG)=DSG*(NSG-0.5D0)
      ENDDO
      DO NSG=1,NSGMAX+1
         SIGG(NSG)=DSG*(NSG-1.D0)
      ENDDO
      DO NTG=1,NTGMAX
         THGM(NTG)=DTG*(NTG-0.5D0)
      ENDDO
      DO NTG=1,NTGMAX+1
         THGG(NTG)=DTG*(NTG-1.D0)
      ENDDO
      RETURN
      END
C
C   ************************************************ 
C   **           Initial Psi for tokamak          **
C   ************************************************
C
      SUBROUTINE EQPSIN
C
      USE libspl1d
      INCLUDE '../eq/eqcomc.inc'
      DIMENSION DERIV(NRVM)
C
      IF(MODELG.NE.5) THEN
         RAXIS=RR
         ZAXIS=0.0D0
      ENDIF
C
C     --- assuming elliptic crosssection, flat current profile ---
C
      IF(MOD(MDLEQF,5).EQ.4) THEN
         PSITA=PI*RKAP*RA**2*BB
         PSIPA=PSITA*2/QA
         PSI0=-PSIPA
      ELSE
         PSI0=-0.5D0*RMU0*RIP*1.D6*RR
         PSIPA=-PSI0
         PSITA=PI*RKAP*RA**2*BB
      ENDIF
C
      DRHO=1.D0/(NRVMAX-1)
      DO NRV=1,NRVMAX
         RHOL=(NRV-1)*DRHO
         PSIPNV(NRV)=RHOL*RHOL
         PSIPV(NRV)=PSIPA*RHOL*RHOL
         PSITV(NRV)=PSITA*RHOL*RHOL
         QPV(NRV)=PSITA/PSIPA
         TTV(NRV)=2.D0*PI*RR*BB
      ENDDO
C
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
          PSI(NTG,NSG)=PSI0*(1.D0-SIGM(NSG)*SIGM(NSG))
          DELPSI(NTG,NSG)=0.D0
          HJT(NTG,NSG)=0.D0
      ENDDO
      ENDDO
C
      CALL SPL1D(PSIPNV,PSITV,DERIV,UPSITV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PSITV: IERR=',IERR
      CALL SPL1D(PSIPNV,QPV,DERIV,UQPV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for QPV: IERR=',IERR
      CALL SPL1D(PSIPNV,TTV,DERIV,UTTV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for TTV: IERR=',IERR
C
C      WRITE(6,'(A,I5,1P4E12.4)')
C     &     ('NRV:',NRV,PSIPNV(NRV),PSITV(NRV),QPV(NRV),TTV(NRV),
C     &      NRV=1,NRVMAX)
C
      RETURN
      END
C
C   ************************************************ 
C   **            Initial Psi for RFP             **
C   ************************************************
C
      SUBROUTINE EQPSIR
C
      USE libspl1d
      INCLUDE '../eq/eqcomc.inc'
      DIMENSION DERIV(NRVM)
C
      IF(MODELG.NE.5) THEN
         RAXIS=RR
         ZAXIS=0.0D0
      ENDIF
C
C     --- assuming elliptic crosssection, flat current profile ---
C
      IF(MOD(MDLEQF,5).EQ.4) THEN
         PSITA=PI*RKAP*RA**2*BB
         PSIPA=PSITA*2/QA
         PSI0=-PSIPA
      ELSE
         PSI0=-0.5D0*RMU0*RIP*1.D6*RR
         PSIPA=-PSI0
         PSITA=PI*RKAP*RA**2*BB
      ENDIF
C
      DRHO=1.D0/(NRVMAX-1)
      DO NRV=1,NRVMAX
         RHOL=(NRV-1)*DRHO
         PSIPNV(NRV)=RHOL*RHOL
         PSIPV(NRV)=PSIPA*RHOL*RHOL
         PSITV(NRV)=PSITA*RHOL*RHOL
         QPV(NRV)=PSITA/PSIPA
         TTV(NRV)=2.D0*PI*RR*BB
      ENDDO
C
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
          PSI(NTG,NSG)=PSI0*(1-SIGM(NSG)*SIGM(NSG))
          DELPSI(NTG,NSG)=0.D0
          HJT(NTG,NSG)=0.D0
      ENDDO
      ENDDO
C
      CALL SPL1D(PSIPNV,PSITV,DERIV,UPSITV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PSITV: IERR=',IERR
      CALL SPL1D(PSIPNV,QPV,DERIV,UQPV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for QPV: IERR=',IERR
      CALL SPL1D(PSIPNV,TTV,DERIV,UTTV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for TTV: IERR=',IERR
C
C      WRITE(6,'(A,I5,1P4E12.4)')
C     &     ('NRV:',NRV,PSIPNV(NRV),PSITV(NRV),QPV(NRV),TTV(NRV),
C     &      NRV=1,NRVMAX)
C
      RETURN
      END
C
C   ************************************************
C   **  Initialize PSI(sigma,theta) from PSIRZ    **
C   ************************************************
C
      SUBROUTINE EQPSI_FROM_PSIRZ(IERR)
C
      INCLUDE '../eq/eqcomc.inc'
C
      IERR=0
C
      PSI0=PSIG(RAXIS,ZAXIS)
      PSIPA=-PSI0

      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         PSI(NTG,NSG)=PSIG(RMM(NTG,NSG),ZMM(NTG,NSG))
         DELPSI(NTG,NSG)=0.D0
         HJT(NTG,NSG)=0.D0
      ENDDO
      ENDDO

      WRITE(6,'(A,1P2E12.4)')
     &   '## EQPSI_FROM_PSIRZ: PSI0,PSIPA=',PSI0,PSIPA

      RETURN
      END
C
C   ************************************************ 
C   **          Boundary Definition               **
C   ************************************************
C
      SUBROUTINE EQDEFB
C
      USE libbrent
      INCLUDE '../eq/eqcomc.inc'
C
      DIMENSION DRHOM(NTGM),DRHOG(NTGMP)
      EXTERNAL EQFBND
C
C     ------ Define criteria for the Brent method ------
C
      EPSZ=1.D-8
C
C     ------ Define boundary radius RHOM/G ------
C
      IERRSU=1
      IF(MODELG.EQ.5.AND.NSUMAX.GE.8.AND.IUSELCFS.NE.0) THEN
         CALL EQSET_RHOB_FROM_RSU(IERRSU)
         IF(IERRSU.EQ.0) THEN
            WRITE(6,*) '## EQDEFB: using gfile LCFS boundary (RSU/ZSU)'
         ELSE
            WRITE(6,*) '## EQDEFB: LCFS boundary map failed, fallback'
            WRITE(6,*) '## EQDEFB: EQSET_RHOB_FROM_RSU IERR=',IERRSU
         ENDIF
      ELSEIF(MODELG.EQ.5.AND.IUSELCFS.EQ.0) THEN
         WRITE(6,*) '## EQDEFB: forcing EQFBND boundary (IUSELCFS=0)'
      ENDIF
C
      IF(IERRSU.NE.0) THEN
         DO NTG=1,NTGMAX
            ZBRF=TAN(THGM(NTG))
            THDASH=FBRENT(EQFBND,THGM(NTG)-1.0D0,THGM(NTG)+1.0D0,EPSZ)
            RHOM(NTG)=RA*SQRT(COS(THDASH+RDLT*SIN(THDASH))**2
     &                       +RKAP**2*SIN(THDASH)**2)
C
            ZBRF=TAN(THGG(NTG))
            THDASH=FBRENT(EQFBND,THGG(NTG)-1.0D0,THGG(NTG)+1.0D0,EPSZ)
            RHOG(NTG)=RA*SQRT(COS(THDASH+RDLT*SIN(THDASH))**2
     &                       +RKAP**2*SIN(THDASH)**2)
         ENDDO
      ENDIF
      RHOG(NTGMAX+1)=RHOG(1)
C
C     ------ Define Delta theta on the boundary ------
C
      DRHOM(1)=0.5D0*(RHOM(2)-RHOM(NTGMAX))/DTG
      DO NTG=2,NTGMAX-1
         DRHOM(NTG)=0.5D0*(RHOM(NTG+1)-RHOM(NTG-1))/DTG
      ENDDO
      DRHOM(NTGMAX)=0.5D0*(RHOM(1)-RHOM(NTGMAX-1))/DTG
C
      DRHOG(1)=0.5D0*(RHOG(2)-RHOG(NTGMAX))/DTG
      DO NTG=2,NTGMAX-1
         DRHOG(NTG)=0.5D0*(RHOG(NTG+1)-RHOG(NTG-1))/DTG
      ENDDO
      DRHOG(NTGMAX)=0.5D0*(RHOG(1)-RHOG(NTGMAX-1))/DTG
      DRHOG(NTGMAX+1)=DRHOG(1)
C
C     ------ Calculate major radius on grid points ------
C
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX+1
         RMG(NSG,NTG)=RR+SIGM(NSG)*RHOG(NTG)*COS(THGG(NTG))
      ENDDO
      ENDDO
C
      DO NSG=1,NSGMAX+1
      DO NTG=1,NTGMAX
         RGM(NSG,NTG)=RR+SIGG(NSG)*RHOM(NTG)*COS(THGM(NTG))       
      ENDDO
      ENDDO
C
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         RMM(NTG,NSG)=RR+SIGM(NSG)*RHOM(NTG)*COS(THGM(NTG))
         ZMM(NTG,NSG)=   SIGM(NSG)*RHOM(NTG)*SIN(THGM(NTG))
      ENDDO
      ENDDO
C
C     ------ Calculate factors for coefficient matrix ------
C
      DO NSG=1,NSGMAX+1
      DO NTG=1,NTGMAX
         AA(NSG,NTG)=SIGG(NSG)/RGM(NSG,NTG)
     &              +(DRHOM(NTG)*DRHOM(NTG)*SIGG(NSG))
     &              /(RGM(NSG,NTG)*RHOM(NTG)*RHOM(NTG))
         AB(NSG,NTG)=-DRHOM(NTG)/(RGM(NSG,NTG)*RHOM(NTG))
      ENDDO
      ENDDO
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX+1
         AC(NSG,NTG)=-DRHOG(NTG)/(RMG(NSG,NTG)*RHOG(NTG))
         AD(NSG,NTG)=1/(RMG(NSG,NTG)*SIGM(NSG))
      ENDDO
      ENDDO
C
C      DO NSG=1,2
C         DO NTG=1,NTGMAX
C            WRITE(6,'(2I5,1P4E12.4)') NSG,NTG,AA(NSG,NTG),AB(NSG,NTG),
C     &                                       AC(NSG,NTG),AD(NSG,NTG)
C         ENDDO
C      ENDDO
C      PAUSE

      RETURN
      END

C   ************************************************
C   **  Boundary radius from gfile LCFS points    **
C   ************************************************
C
      SUBROUTINE EQSET_RHOB_FROM_RSU(IERR)
C
      INCLUDE '../eq/eqcomc.inc'
      DIMENSION THS(NSUM),RHS(NSUM)
      DIMENSION THU(NSUM),RHU(NSUM)
      DIMENSION THP(NSUM+2),RHP(NSUM+2)
C
      IERR=0
      IF(NSUMAX.LT.3) THEN
         IERR=1
         RETURN
      ENDIF
C
      N=0
      DO NSU=1,NSUMAX
         DX=RSU(NSU)-RAXIS
         DZ=ZSU(NSU)-ZAXIS
         RL=SQRT(DX*DX+DZ*DZ)
         IF(RL.LE.1.D-8) CYCLE
         TL=ATAN2(DZ,DX)
         IF(TL.LT.0.D0) TL=TL+2.D0*PI
         N=N+1
         THS(N)=TL
         RHS(N)=RL
      ENDDO
C
      IF(N.LT.3) THEN
         IERR=2
         RETURN
      ENDIF
C
C     --- Insertion sort by theta ---
      DO I=2,N
         TWS=THS(I)
         RWS=RHS(I)
         J=I-1
 100     CONTINUE
         IF(J.GE.1) THEN
            IF(THS(J).GT.TWS) THEN
               THS(J+1)=THS(J)
               RHS(J+1)=RHS(J)
               J=J-1
               GOTO 100
            ENDIF
         ENDIF
         THS(J+1)=TWS
         RHS(J+1)=RWS
      ENDDO
C
C     --- Merge duplicate/near-duplicate angles ---
      NU=1
      THU(1)=THS(1)
      RHU(1)=RHS(1)
      DO I=2,N
         IF(THS(I)-THU(NU).GT.1.D-8) THEN
            NU=NU+1
            THU(NU)=THS(I)
            RHU(NU)=RHS(I)
         ELSE
            RHU(NU)=0.5D0*(RHU(NU)+RHS(I))
         ENDIF
      ENDDO
C
      IF(NU.LT.3) THEN
         IERR=3
         RETURN
      ENDIF
C
      THP(1)=THU(NU)-2.D0*PI
      RHP(1)=RHU(NU)
      DO I=1,NU
         THP(I+1)=THU(I)
         RHP(I+1)=RHU(I)
      ENDDO
      THP(NU+2)=THU(1)+2.D0*PI
      RHP(NU+2)=RHU(1)
C
C     --- Interpolate rho(theta) on EQ mesh angles ---
      DO NTG=1,NTGMAX
         CALL EQRHOINTERP(THGM(NTG),THP,RHP,NU+2,RL,IERRL)
         IF(IERRL.NE.0.OR.RL.LE.1.D-8) THEN
            IERR=4
            RETURN
         ENDIF
         RHOM(NTG)=RL
      ENDDO
C
      DO NTG=1,NTGMAX
         CALL EQRHOINTERP(THGG(NTG),THP,RHP,NU+2,RL,IERRL)
         IF(IERRL.NE.0.OR.RL.LE.1.D-8) THEN
            IERR=5
            RETURN
         ENDIF
         RHOG(NTG)=RL
      ENDDO
C
      RETURN
      END

      SUBROUTINE EQRHOINTERP(THETA,THP,RHP,NP,RHOL,IERR)
       INCLUDE '../eq/eqcomc.inc'
       DIMENSION THP(NP),RHP(NP)
C
      IERR=0
      TL=THETA
      IF(TL.LT.0.D0) TL=TL+2.D0*PI
      IF(TL.GE.2.D0*PI) TL=TL-2.D0*PI
C
      J=1
 200  IF(J.LT.NP.AND.TL.GT.THP(J+1)) THEN
         J=J+1
         GOTO 200
      ENDIF
      IF(J.GE.NP) THEN
         IERR=1
         RHOL=0.D0
         RETURN
      ENDIF
C
      DT=THP(J+1)-THP(J)
      IF(ABS(DT).LE.1.D-12) THEN
         IERR=2
         RHOL=0.D0
         RETURN
      ENDIF
      W=(TL-THP(J))/DT
       RHOL=(1.D0-W)*RHP(J)+W*RHP(J+1)
       RETURN
       END
C
      SUBROUTINE EQGET_RHOB(THETA,RHOL,IERR)
      INCLUDE '../eq/eqcomc.inc'
C
      IERR=0
      TL=THETA
      IF(TL.LT.0.D0) TL=TL+2.D0*PI
      IF(TL.GE.2.D0*PI) TL=TL-2.D0*PI
C
      DO NTG=1,NTGMAX
         TH0=THGG(NTG)
         TH1=THGG(NTG+1)
         IF(TL.GE.TH0.AND.TL.LE.TH1) THEN
            DT=TH1-TH0
            IF(ABS(DT).LE.1.D-12) THEN
               IERR=2
               RHOL=0.D0
               RETURN
            ENDIF
            W=(TL-TH0)/DT
            RHOL=(1.D0-W)*RHOG(NTG)+W*RHOG(NTG+1)
            RETURN
         ENDIF
      ENDDO
C
      IERR=1
      RHOL=0.D0
      RETURN
      END
C
C   ************************************************
C   **         Boundary shape function            **
C   ************************************************
C
      FUNCTION EQFBND(X)
C      
      INCLUDE '../eq/eqcomc.inc'
C
      EQFBND=ZBRF*COS(X+RDLT*SIN(X))-RKAP*SIN(X)
      RETURN
      END
C
C   ************************************************
C   **               Iteration Loop               **
C   ************************************************
C
      SUBROUTINE EQLOOP(IERR)
C
      INCLUDE '../eq/eqcomc.inc'
      REAL*8 PSIB,DELB
C
      IERR=0
C
      DO NLOOP=1,NLPMAX
         CALL EQBAND
         CALL EQRHSV(IERR)
         IF(IERR.NE.0) RETURN
         CALL EQSOLV
C
         SUM0=0.D0
         SUM1=0.D0
         IBAD=0
         DO NSG=1,NSGMAX
             DO NTG=1,NTGMAX
                IF(PSI(NTG,NSG).NE.PSI(NTG,NSG).OR.
     &             DELPSI(NTG,NSG).NE.DELPSI(NTG,NSG)) THEN
                   IBAD=1
                   NSGB=NSG
                   NTGB=NTG
                   PSIB=PSI(NTG,NSG)
                   DELB=DELPSI(NTG,NSG)
                   GOTO 210
                ENDIF
                SUM0=SUM0+PSI(NTG,NSG)**2
                SUM1=SUM1+DELPSI(NTG,NSG)**2
             ENDDO
          ENDDO
  210     CONTINUE
          IF(IBAD.NE.0) THEN
             WRITE(6,'(A,I5,A,I5)')
     &          'XX EQLOOP: NaN in PSI/DELPSI at NSG=',
     &          NSGB,', NTG=',NTGB
             WRITE(6,'(A,1P4E12.4)')
     &          '   PSI,DELPSI,RAXIS,ZAXIS=',PSIB,DELB,RAXIS,ZAXIS
             IERR=901
             RETURN
          ENDIF
          IF(SUM0.LE.1.D-30.OR.SUM0.NE.SUM0.OR.SUM1.NE.SUM1) THEN
             WRITE(6,'(A,1P4E12.4)')
     &          'XX EQLOOP: invalid SUM0/SUM1,RAXIS,ZAXIS=',
     &          SUM0,SUM1,RAXIS,ZAXIS
             IERR=902
             RETURN
          ENDIF
          SUML=SQRT(SUM1/SUM0)
          IF(SUML.NE.SUML) THEN
             WRITE(6,'(A,1P4E12.4)')
     &          'XX EQLOOP: SUML is NaN, SUM0,SUM1,RAXIS=',
     &          SUM0,SUM1,RAXIS,ZAXIS
             IERR=903
             RETURN
          ENDIF
C
          IF(SUML.LT.EPSEQ) THEN
            IF(NPRINT.GE.1) THEN
               WRITE(6,'(A,1P4E14.6,I5)')
     &           'SUML,R/ZAXIS,PSI0=',SUML,RAXIS,ZAXIS,PSI0,NLOOP
            ENDIF
            RETURN
         ELSE
            IF((NPRINT.EQ.1.AND.NLOOP.EQ.1).OR.
     &          NPRINT.GE.2) THEN
               WRITE(6,'(A,1P4E14.6)')
     &           'SUML,R/ZAXIS,PSI0=',SUML,RAXIS,ZAXIS,PSI0
            ENDIF
         ENDIF
      ENDDO
C
      IF(NPRINT.GE.1) THEN
         WRITE(6,'(A,1P4E14.6,I5)')
     &           'SUML,R/ZAXIS,PSI0=',SUML,RAXIS,ZAXIS,PSI0,NLOOP
      ENDIF
      WRITE(6,*) 'XX EQLOOP: NLOOP exceeds NLPMAX'
      IERR=100
C
      RETURN
      END
C
C   ***********************************************
C   **            Matrix calculation             **
C   ***********************************************
C
      SUBROUTINE EQBAND
C
      INCLUDE '../eq/eqcomc.inc'
C
C     ------ Define matrix length MMAX and matrix half width NBND ------
C
      MMAX=NSGMAX*NTGMAX
      NBND=2*NTGMAX
C
C     ------ Initialize band matrix coefficients Q ------
C
      DO N=1,MMAX
      DO M=1,2*NBND-1
          Q(M,N)=0.D0
      ENDDO
      ENDDO
C
C     ------ Calculate band matrix coefficients Q ------
C
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         I=(NSG-1)*NTGMAX+NTG
         IF(NTG.EQ.1)THEN
            Q(NBND         -1,I)= (AB(NSG,NTG)+AC(NSG,NTG))
     &                           /(4.D0*DSG*DTG)
            Q(NBND  +NTGMAX-1,I)= (AB(NSG,NTG)-AB(NSG+1,NTG))
     &                           /(4.D0*DSG*DTG)
     &                           +AD(NSG,NTG)/(DTG*DTG)
            Q(NBND+2*NTGMAX-1,I)=-(AB(NSG+1,NTG)+AC(NSG,NTG))
     &                           /(4.D0*DSG*DTG)
         ELSE
            Q(NBND  -NTGMAX-1,I)= (AB(NSG,NTG)+AC(NSG,NTG))
     &                           /(4.D0*DSG*DTG)
            Q(NBND         -1,I)=(AB(NSG,NTG)-AB(NSG+1,NTG))
     &                           /(4.D0*DSG*DTG)
     &                           +AD(NSG,NTG)/(DTG*DTG)
            Q(NBND  +NTGMAX-1,I)=-(AB(NSG+1,NTG)+AC(NSG,NTG))
     &                           /(4.D0*DSG*DTG)
         ENDIF
         Q(NBND-NTGMAX,I)=  AA(NSG,NTG)/(DSG*DSG)
     &                   -(AC(NSG,NTG+1)-AC(NSG,NTG))/(4.D0*DSG*DTG)
         Q(NBND       ,I)=-(AA(NSG+1,NTG)+AA(NSG,NTG))/(DSG*DSG)
     &                   -(AD(NSG,NTG+1)+AD(NSG,NTG))/(DTG*DTG)
         Q(NBND+NTGMAX,I)=  AA(NSG+1,NTG)/(DSG*DSG)
     &                   +(AC(NSG,NTG+1)-AC(NSG,NTG))/(4.D0*DSG*DTG)
C
C     ------ Set periodic condition ------
C
         IF(NTG.EQ.NTGMAX) THEN
            Q(NBND-2*NTGMAX+1,I)=-(AB(NSG,NTG)+AC(NSG,NTG+1))
     &                            /(4.D0*DSG*DTG)
            Q(NBND-  NTGMAX+1,I)= (AB(NSG+1,NTG)-AB(NSG,NTG))
     &                            /(4.D0*DSG*DTG)
     &                          +  AD(NSG,NTG+1)/(DTG*DTG)
            Q(NBND         +1,I)= (AB(NSG+1,NTG)+AC(NSG,NTG+1))
     &                            /(4.D0*DSG*DTG)
         ELSE
            Q(NBND-NTGMAX+1,I)=-(AB(NSG,NTG)+AC(NSG,NTG+1))
     &                          /(4.D0*DSG*DTG)
            Q(NBND       +1,I)= (AB(NSG+1,NTG)-AB(NSG,NTG))
     &                          /(4.D0*DSG*DTG)
     &                        +  AD(NSG,NTG+1)/(DTG*DTG)
            Q(NBND+NTGMAX+1,I)= (AB(NSG+1,NTG)+AC(NSG,NTG+1))
     &                          /(4.D0*DSG*DTG)
         ENDIF
      ENDDO
      ENDDO
C
C     ------ Set radial boundary condition ------
C
      DO I=1,NTGMAX
      DO J=1,NTGMAX/2
         Q(NBND+J-I,I)=Q(NBND+J-I,I)
     &                +Q(NBND+J-I-NTGMAX/2,I)
         Q(NBND+J-I-NTGMAX/2,I)=0.D0
         Q(NBND+J-I+NTGMAX/2,I)=Q(NBND+J-I+NTGMAX/2,I)
     &                         +Q(NBND+J-I-NTGMAX,I)
         Q(NBND+J-I-NTGMAX,I)=0.D0
      ENDDO
      ENDDO
      DO I=1,NTGMAX
      DO J=1,NTGMAX
         Q(NBND+J-I,I+(NSGMAX-1)*NTGMAX)
     &                   =Q(NBND+J-I,I+(NSGMAX-1)*NTGMAX)
     &                   -Q(NBND+J-I+NTGMAX,I+(NSGMAX-1)*NTGMAX)
         Q(NBND+J-I+NTGMAX,I+(NSGMAX-1)*NTGMAX)=0.D0
      ENDDO
      ENDDO
      RETURN
      END
C
C   ************************************************ 
C   **                RHS vector                  **
C   ************************************************
C
      SUBROUTINE EQRHSV(IERR)
C
      USE libspl1d
      INCLUDE '../eq/eqcomc.inc'
      DIMENSION DERIV(NRVM)
C
      IERR=0
C
C     ----- calculate PSIRZ(R,Z) from PSI(sigma, theta) -----
C
      CALL EQTORZ
C
C     ----- calculate flux average from PSIRZ(R,Z) -----
C
      CALL EQCALV(IERR)
      IF(IERR.NE.0) RETURN
C
C     ----- calculate right hand side vector -----
C
      IMDLEQF=MOD(MDLEQF,5)
C
C     ----- Given pressure and toroidal current profiles -----
C
      IF(IMDLEQF.EQ.0) THEN
         RRC=RAXIS
         FJP=0.D0
         FJT=0.D0
         DO NSG=1,NSGMAX
         DO NTG=1,NTGMAX
            PSIPNL=1.D0-PSI(NTG,NSG)/PSI0
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQJPSI(PSIPNL,HJPSID,HJPSI)
            CALL EQTPSI(PSIPNL,TPSI,DTPSI)
            CALL EQOPSI(PSIPNL,OMGPSI,DOMGPSI)
            PP(NTG,NSG)=PPSI
            HJP(NTG,NSG)=-2.D0*PI*RMM(NTG,NSG)*DPPSI
            HJP1(NTG,NSG)=-2.D0*PI*RMM(NTG,NSG)*DPPSI
            HJT1(NTG,NSG)=-HJPSI
            HJP2A=EXP(RMM(NTG,NSG)**2*OMGPSI**2*AMP/(2.D0*TPSI))
            HJP2B=EXP(RRC**2         *OMGPSI**2*AMP/(2.D0*TPSI))
            HJP2C=HJP2A-(RRC**2/RMM(NTG,NSG)**2)*HJP2B
            HJP2D=HJP2A-(RRC**4/RMM(NTG,NSG)**4)*HJP2B
            HJP2E=-0.5D0*2.D0*PI*PPSI*RMM(NTG,NSG)**3
            HJP2F=AMP*(2.D0*OMGPSI*DOMGPSI/TPSI
     &           -DTPSI*OMGPSI**2/TPSI**2)
            HJP2(NTG,NSG)=HJP2C*HJP1(NTG,NSG)+HJP2D*HJP2E*HJP2F
            HJT2(NTG,NSG)=(RRC/RMM(NTG,NSG))*HJT1(NTG,NSG)
            DVOL=SIGM(NSG)*RHOM(NTG)*RHOM(NTG)*DSG*DTG
C
            FJP=FJP+HJP2(NTG,NSG)*DVOL
            FJT=FJT+HJT2(NTG,NSG)*DVOL
         ENDDO
         ENDDO
         TJ=(RIP*1.D6-FJP)/FJT
C
         DO NSG=1,NSGMAX
         DO NTG=1,NTGMAX
            HJT(NTG,NSG)=HJP2(NTG,NSG)+TJ*HJT2(NTG,NSG)
            PSIPNL=1.D0-PSI(NTG,NSG)/PSI0
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQJPSI(PSIPNL,HJPSID,HJPSI)
            CALL EQTPSI(PSIPNL,TPSI,DTPSI)
            CALL EQOPSI(PSIPNL,OMGPSI,DOMGPSI)
            TT(NTG,NSG)=SQRT((2.D0*PI*BB*RR)**2
     &                 +2.D0*RMU0*RRC
     &                 *(2.D0*PI*TJ*HJPSID-RRC*PPSI
     &                 *EXP(RRC**2*OMGPSI**2*AMP/(2.D0*TPSI))))
            RHO(NTG,NSG)=(PPSI*AMP/TPSI)
     &                  *EXP(RMM(NTG,NSG)**2*OMGPSI**2*AMP/(2.D0*TPSI))
         ENDDO
         ENDDO
         RIPX=RIP
C
C         DO NSG=1,3
C            DO NTG=1,NTGMAX
C               WRITE(6,'(2I5,1P3E12.4)')
C     &              NSG,NTG,RMM(NTG,NSG),ZMM(NTG,NSG),PSI(NTG,NSG)
C               WRITE(6,'(2I5,1P3E12.4)')
C     &              NSG,NTG,HJP1(NTG,NSG),HJT1(NTG,NSG),PP(NTG,NSG)
C               WRITE(6,'(2I5,1P3E12.4)')
C     &              NSG,NTG,HJP2(NTG,NSG),HJT2(NTG,NSG),TT(NTG,NSG)
C            ENDDO
C         ENDDO
C         PAUSE
C
C     ----- Given pressure and poloidal current profiles -----
C
      ELSEIF(IMDLEQF.EQ.1) THEN
         FJP=0.D0
         FJT1=0.D0
         FJT2=0.D0
         DO NSG=1,NSGMAX
         DO NTG=1,NTGMAX
            PSIPNL=1.D0-PSI(NTG,NSG)/PSI0
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQFPSI(PSIPNL,FPSI,DFPSI)
            PP(NTG,NSG)=PPSI
            HJP(NTG,NSG)=-2.D0*PI*RMM(NTG,NSG)*DPPSI
            HJP1(NTG,NSG)=-2.D0*PI*RMM(NTG,NSG)*DPPSI
            HJT1(NTG,NSG)=-2.D0*PI*BB*RR       *DFPSI
     &                                /(2.D0*PI*RMU0*RMM(NTG,NSG))
            HJT2(NTG,NSG)=-(FPSI-2.D0*PI*BB*RR)*DFPSI
     &                                /(2.D0*PI*RMU0*RMM(NTG,NSG))
            HJP2(NTG,NSG)=HJP1(NTG,NSG)
            DVOL=SIGM(NSG)*RHOM(NTG)*RHOM(NTG)*DSG*DTG
            FJP =FJP +HJP2(NTG,NSG)*DVOL
            FJT1=FJT1+HJT1(NTG,NSG)*DVOL
            FJT2=FJT2+HJT2(NTG,NSG)*DVOL
         ENDDO
         ENDDO
         IF(RIP.LE.0.D0) THEN
C           --- No Ip matching: use F directly (TJ=1) ---
            TJ=1.D0
         ELSEIF(FJT1.GT.0.D0) THEN
            TJ=(-FJT1+SQRT(FJT1**2+4.D0*FJT2*(RIP*1.D6-FJP)))
     &         /(2.D0*FJT2)
         ELSE
            TJ=(-FJT1-SQRT(FJT1**2+4.D0*FJT2*(RIP*1.D6-FJP)))
     &         /(2.D0*FJT2)
         ENDIF
C
         DO NSG=1,NSGMAX
         DO NTG=1,NTGMAX
            HJT(NTG,NSG)=HJP2(NTG,NSG)+TJ*HJT1(NTG,NSG)
     &                                +TJ*TJ*HJT2(NTG,NSG)
            PSIPNL=1.D0-PSI(NTG,NSG)/PSI0
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQFPSI(PSIPNL,FPSI,DFPSI)
            CALL EQTPSI(PSIPNL,TPSI,DTPSI)
            TT(NTG,NSG)=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
            RHO(NTG,NSG)=0.D0
         ENDDO
         ENDDO
         IF(RIP.LE.0.D0) THEN
            RIPX=(FJP+TJ*FJT1+TJ*TJ*FJT2)*1.D-6
         ELSE
            RIPX=RIP
         ENDIF
C
C     ----- Given pressure and parallel current profiles -----
C
      ELSEIF(IMDLEQF.EQ.2) THEN
         CALL EQIPJP
C         DO NR=1,11
C            PSIPNL=0.002D0*(NR-1)
C            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
C            WRITE(6,'(A,I5,1P3E12.4)') 'NR,PSIPNL,PPSI,DPPSI=',
C     &           NR,PSIPNL,PPSI,DPPSI
C         ENDDO
         FJP=0.D0
         FJT1=0.D0
         FJT2=0.D0
         DO NSG=1,NSGMAX
         DO NTG=1,NTGMAX
            PSIPNL=1.D0-PSI(NTG,NSG)/PSI0
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQFIPV(PSIPNL,FPSI,DFPSI)
            PP(NTG,NSG)=PPSI
            HJP(NTG,NSG)=-2.D0*PI*RMM(NTG,NSG)*DPPSI
            HJP1(NTG,NSG)=-2.D0*PI*RMM(NTG,NSG)*DPPSI
            HJT1(NTG,NSG)=-2.D0*PI*BB*RR       *DFPSI
     &                              /(2.D0*PI*RMU0*RMM(NTG,NSG))
            HJT2(NTG,NSG)=-(FPSI-2.D0*PI*BB*RR)*DFPSI
     &                              /(2.D0*PI*RMU0*RMM(NTG,NSG))
            HJP2(NTG,NSG)=HJP1(NTG,NSG)
            DVOL=SIGM(NSG)*RHOM(NTG)*RHOM(NTG)*DSG*DTG
            FJP=FJP+HJP2(NTG,NSG)*DVOL
            FJT1=FJT1+HJT1(NTG,NSG)*DVOL
            FJT2=FJT2+HJT2(NTG,NSG)*DVOL
         ENDDO
         ENDDO
         IF(FJT1.GT.0.D0) THEN
            TJ=(-FJT1+SQRT(FJT1**2+4.D0*FJT2*(RIP*1.D6-FJP)))
     &         /(2.D0*FJT2)
         ELSE
            TJ=(-FJT1-SQRT(FJT1**2+4.D0*FJT2*(RIP*1.D6-FJP)))
     &         /(2.D0*FJT2)
         ENDIF
C
         DO NSG=1,NSGMAX
         DO NTG=1,NTGMAX
            HJT(NTG,NSG)=HJP2(NTG,NSG)+TJ*HJT1(NTG,NSG)
     &                                +TJ*TJ*HJT2(NTG,NSG)
            PSIPNL=1.D0-PSI(NTG,NSG)/PSI0
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQFPSI(PSIPNL,FPSI,DFPSI)
            TT(NTG,NSG)=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
            RHO(NTG,NSG)=0.D0
         ENDDO
         ENDDO 
         RIPX=RIP
C
C     ----- Given pressure and parallel current profiles -----
C
      ELSEIF(IMDLEQF.EQ.3) THEN
         CALL EQIPJP
C         DO NR=1,11
C            PSIPNL=0.002D0*(NR-1)
C            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
C            CALL EQFIPV(PSIPNL,FPSI,DFPSI)
C            WRITE(6,'(A,I5,1P5E12.4)') 'PSI,PP,FFF=',
C     &           NR,PSIPNL,PPSI,DPPSI,FPSI,DFPSI
C         ENDDO
         FJP=0.D0
         FJT=0.D0
         DO NSG=1,NSGMAX
         DO NTG=1,NTGMAX
            PSIPNL=1.D0-PSI(NTG,NSG)/PSI0
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQFIPV(PSIPNL,FPSI,DFPSI)
            PP(NTG,NSG)=PPSI
            HJP(NTG,NSG)=-2.D0*PI*RMM(NTG,NSG)*DPPSI
C            IF(NTG.EQ.1.OR.NTG.EQ.NTGMAX/2+1) THEN
C               WRITE(6,'(A,I5,1P4E12.4)') 
C     &              'NSG,PSIPNL,RMM,DPPSI,HJP=',
C     &               NSG,PSIPNL,RMM(NTG,NSG),DPPSI,HJP(NTG,NSG)
C            ENDIF
            HJP1(NTG,NSG)=-2.D0*PI*RMM(NTG,NSG)*DPPSI
            HJT1(NTG,NSG)=-FPSI*DFPSI/(2.D0*PI*RMU0*RMM(NTG,NSG))
            HJP2(NTG,NSG)=HJP1(NTG,NSG)
            DVOL=SIGM(NSG)*RHOM(NTG)*RHOM(NTG)*DSG*DTG
            FJP=FJP+HJP2(NTG,NSG)*DVOL
            FJT=FJT+HJT1(NTG,NSG)*DVOL
C
            HJT(NTG,NSG)=HJP2(NTG,NSG)+HJT1(NTG,NSG)
            TT(NTG,NSG)=FPSI
            RHO(NTG,NSG)=0.D0
         ENDDO
         ENDDO 
         RIPX=(FJP+FJT)*1.D-6
C         WRITE(6,'(A,1P3E12.4)') 'RIPX=',RIPX,FJP*1.D-6,FJT*1.D-6
C
C     ----- Given pressure and safety factor profiles ------
C
      ELSEIF(IMDLEQF.EQ.4) THEN
         CALL EQIPQP
         FJP=0.D0
         FJT=0.D0
         DO NSG=1,NSGMAX
         DO NTG=1,NTGMAX
            PSIPNL=1.D0-PSI(NTG,NSG)/PSI0
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQFIPV(PSIPNL,FPSI,DFPSI)
            PP(NTG,NSG)=PPSI
            HJP(NTG,NSG)=-2.D0*PI*RMM(NTG,NSG)*DPPSI
            HJP1(NTG,NSG)=-2.D0*PI*RMM(NTG,NSG)*DPPSI
            HJT1(NTG,NSG)=-FPSI*DFPSI/(2.D0*PI*RMU0*RMM(NTG,NSG))
            HJP2(NTG,NSG)=HJP1(NTG,NSG)
            DVOL=SIGM(NSG)*RHOM(NTG)*RHOM(NTG)*DSG*DTG
            FJP=FJP+HJP2(NTG,NSG)*DVOL
            FJT=FJT+HJT1(NTG,NSG)*DVOL
C
            HJT(NTG,NSG)=HJP2(NTG,NSG)+HJT1(NTG,NSG)
            TT(NTG,NSG)=FPSI
            RHO(NTG,NSG)=0.D0
         ENDDO
         ENDDO 
         RIPX=(FJP+FJT)*1.D-6
C         WRITE(6,'(A,1P3E12.4)') 'RIPX=',RIPX,FJP*1.D-6,FJT*1.D-6
      ENDIF
C
C     ----- CALCULATE POLOIDAL CURRENT -----
C
      IMDLEQF=MOD(MDLEQF,5)
      DO NRV=1,NRVMAX
         PSIPNL=PSIPNV(NRV)
         IF (IMDLEQF.EQ.0) THEN
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQJPSI(PSIPNL,HJPSID,HJPSI)
            CALL EQTPSI(PSIPNL,TPSI,DTPSI)
            CALL EQOPSI(PSIPNL,OMGPSI,DOMGPSI)
            TTVL=SQRT((2.D0*PI*BB*RR)**2
     &                  +2.D0*RMU0*RRC
     &                  *(2.D0*PI*TJ*HJPSID-RRC*PPSI
     &                  *EXP(RRC**2*OMGPSI**2*AMP/(2.D0*TPSI))))
         ELSEIF (IMDLEQF.EQ.1) THEN
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQFPSI(PSIPNL,FPSI,DFPSI)
            TTVL=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
         ELSEIF (IMDLEQF.EQ.2) THEN
            CALL EQFIPV(PSIPNL,FPSI,DFPSI)
            TTVL=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
         ELSEIF (IMDLEQF.EQ.3) THEN
            CALL EQFIPV(PSIPNL,FPSI,DFPSI)
            TTVL=FPSI
         ELSEIF (IMDLEQF.EQ.4) THEN
            CALL EQFIPV(PSIPNL,FPSI,DFPSI)
            TTVL=FPSI
         ENDIF
         TTV(NRV)=TTVL
      ENDDO
C
      CALL SPL1D(PSIPNV,PSITV,DERIV,UPSITV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PSITV: IERR=',IERR
      CALL SPL1D(PSIPNV,QPV,DERIV,UQPV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for QPV: IERR=',IERR
      CALL SPL1D(PSIPNV,TTV,DERIV,UTTV,NRVMAX,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for TTV: IERR=',IERR
C
C      WRITE(6,'(A,I5,1P4E12.4)')
C     &     ('NRV:',NRV,PSIPNV(NRV),PSITV(NRV),QPV(NRV),TTV(NRV),
C     &      NRV=1,NRVMAX)
C
      RETURN
      END
C
C   ************************************************
C   **              Matrix Solver                 **
C   ************************************************
C
      SUBROUTINE EQSOLV
C
      USE libbnd
      INCLUDE '../eq/eqcomc.inc'
C
      REAL(rkind),ALLOCATABLE:: FJT(:),PSIOLD(:,:)

      ALLOCATE(FJT(MLM),PSIOLD(NTGM,NSGM))
C
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         PSIOLD(NTG,NSG)=PSI(NTG,NSG)
      ENDDO
      ENDDO
C
      DO NSG=1,NSGMAX
         I=(NSG-1)*NTGMAX
         DO NTG=1,NTGMAX
            FJT(I+NTG)=2.D0*PI*RMU0*HJT(NTG,NSG)
     &                 *SIGM(NSG)*RHOM(NTG)*RHOM(NTG)
         ENDDO
      ENDDO
C
C      DO NSG=1,2
C         DO NTG=1,NTGMAX
C            WRITE(6,'(2I5,1PE12.4)') NSG,NTG,HJT(NSG,NTG)
C         ENDDO
C      ENDDO
C      PAUSE
C
C      DO I=1,3
C         WRITE(6,'(1p5E12.4)') FJT(I),(Q(J,I),J=1,4*NTGMAX-1)
C      ENDDO
C
      CALL BANDRD(Q,FJT,NTGMAX*NSGMAX,4*NTGMAX-1,MWM,IERR)
         IF(IERR.NE.0) THEN
            WRITE(6,*) 'XX EQSOLV: BANDRD ERROR: IERR = ',IERR
         ENDIF
C
C      WRITE(6,'(1p5E12.4)') (FJT(I),I=1,2*NTGMAX)
C
      DO NSG=1,NSGMAX
         I=(NSG-1)*NTGMAX
         DO NTG=1,NTGMAX
            PSI(NTG,NSG)=FJT(I+NTG)
         ENDDO
      ENDDO
C      WRITE(6,'(1p5E12.4)') (PSI(I,1),I=1,NTGMAX)
C      WRITE(6,'(1p5E12.4)') (PSI(I,2),I=1,NTGMAX)
C      WRITE(6,'(1p5E12.4)') (PSI(I,3),I=1,NTGMAX)
C      WRITE(6,'(1p5E12.4)') (PSI(I,4),I=1,NTGMAX)
C     
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         DELPSI(NTG,NSG)=PSI(NTG,NSG)-PSIOLD(NTG,NSG)
      ENDDO
      ENDDO
      RETURN
      END
C
C   ************************************************
C   **          sigma,theta to R,Z                **
C   ************************************************
C
      SUBROUTINE EQTORZ
C
       USE libbrent
       INCLUDE '../eq/eqcomc.inc'
       EXTERNAL EQFBND
C
      RMIN= RR-RB
      RMAX= RR+RB
      ZMIN=-RKAP*RB
      ZMAX= RKAP*RB
C
      SIG1=1.D0
      SIG2=0.95D0
C
      DTRG=(RMAX-RMIN)/(NRGMAX-1)
      DTZG=(ZMAX-ZMIN)/(NZGMAX-1)
      EPSZ=1.D-8
C
      DO NRG=1,NRGMAX
         RG(NRG)=RMIN+DTRG*(NRG-1)
      ENDDO
      DO NZG=1,NZGMAX
         ZG(NZG)=ZMIN+DTZG*(NZG-1)
      ENDDO
C
      CALL EQSETF
C
      DO NRG=1,NRGMAX
      DO NZG=1,NZGMAX
         IF((RG(NRG)-RR.EQ.0.D0).AND.
     &      (ZG(NZG).EQ.0.D0)) THEN
            RHOL=0.D0
            SIGL=0.D0
         ELSE
            THL=ATAN2(ZG(NZG),RG(NRG)-RR)
            IF(THL.LT.0.D0) THL=THL+2.D0*PI
            IF(MODELG.EQ.5.AND.NSUMAX.GE.8.AND.IUSELCFS.NE.0) THEN
               CALL EQGET_RHOB(THL,RHOL,IERRL)
               IF(IERRL.NE.0.OR.RHOL.LE.1.D-8) THEN
                  ZBRF=TAN(THL)
                  THDASH=FBRENT(EQFBND,THL-1.0D0,THL+1.0D0,EPSZ)
                  RHOL=RA*SQRT(COS(THDASH+RDLT*SIN(THDASH))**2
     &                        +RKAP**2*SIN(THDASH)**2)
               ENDIF
            ELSE
               ZBRF=TAN(THL)
               THDASH=FBRENT(EQFBND,THL-1.0D0,THL+1.0D0,EPSZ)
               RHOL=RA*SQRT(COS(THDASH+RDLT*SIN(THDASH))**2
     &                     +RKAP**2*SIN(THDASH)**2)
            ENDIF
            SIGL=SQRT((RG(NRG)-RR)**2+ZG(NZG)**2)/RHOL
         ENDIF
         IF(SIGL.LT.1.D0) THEN
            PSIRZ(NRG,NZG)=PSIF(SIGL,THL)
            HJTRZ(NRG,NZG)=HJTF(SIGL,THL)
         ELSE
            PSI1=PSIF(SIG1,THL)
            PSI2=PSIF(SIG2,THL)
            PSIRZ(NRG,NZG)=PSI2+(PSI1-PSI2)*(SIGL-SIG2)/(SIG1-SIG2)
            HJTRZ(NRG,NZG)=0.D0
         ENDIF
      ENDDO
      ENDDO
C
      RETURN
      END
C
C   ******************************************
C   ***** Calculate Spline Coeff for Psi *****
C   ******************************************
C
      SUBROUTINE EQSETF
C
      USE libspl2d
      INCLUDE '../eq/eqcomc.inc'
C
      REAL(rkind),ALLOCATABLE:: PSISX(:,:),PSITX(:,:)
      REAL(rkind),ALLOCATABLE:: PSISTX(:,:)

      ALLOCATE(PSISX(NTGPM,NSGPM),PSITX(NTGPM,NSGPM))
      ALLOCATE(PSISTX(NTGPM,NSGPM))
C
C     ----- mesh extended in sigma (radius) and theta (periodic) -----
C
      SIGMX(1)=0.D0
      DO NSG=1,NSGMAX
         SIGMX(NSG+1)=SIGM(NSG)
      ENDDO
      SIGMX(NSGMAX+2)=1.D0
C
      THGMX(1)=0.D0
      DO NTG=1,NTGMAX
         THGMX(NTG+1)=THGM(NTG)
      ENDDO
      THGMX(NTGMAX+2)=2.D0*PI
C
      SUMPSI=0.D0
      SUMHJT=0.D0
      DO NTG=1,NTGMAX
         SUMPSI=SUMPSI+(9.D0*PSI(NTG,1)-PSI(NTG,2))/8.D0
         SUMHJT=SUMHJT+(9.D0*HJT(NTG,1)-HJT(NTG,2))/8.D0
      ENDDO
      PSIL=SUMPSI/NTGMAX
      HJTL=SUMHJT/NTGMAX
      DO NTG=1,NTGMAX+2
         PSIST(NTG,1)=PSIL
         HJTST(NTG,1)=HJTL
      ENDDO
C
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         PSIST(NTG+1,NSG+1)=PSI(NTG,NSG)
         HJTST(NTG+1,NSG+1)=HJT(NTG,NSG)
      ENDDO
      ENDDO
C
      DO NSG=1,NSGMAX
         PSIST(       1,NSG+1)=(9.D0*PSI(     1,NSG)-PSI(       2,NSG))
     &                         /16.D0
     &                        +(9.D0*PSI(NTGMAX,NSG)-PSI(NTGMAX-1,NSG))
     &                         /16.D0
         PSIST(NTGMAX+2,NSG+1)=PSIST(1,NSG+1)
         HJTST(       1,NSG+1)=(9.D0*HJT(     1,NSG)-HJT(       2,NSG))
     &                         /16.D0
     &                        +(9.D0*HJT(NTGMAX,NSG)-HJT(NTGMAX-1,NSG))
     &                         /16.D0
         HJTST(NTGMAX+2,NSG+1)=HJTST(1,NSG+1)
      ENDDO
C
      DO NTG=1,NTGMAX+2
         PSIST(NTG,NSGMAX+2)=0.D0
         HJTST(NTG,NSGMAX+2)=0.D0
      ENDDO
C
      CALL SPL2D(THGMX,SIGMX,PSIST,PSITX,PSISX,PSISTX,UPSIST,
     &           NTGPM,NTGMAX+2,NSGMAX+2,4,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL2D for PSIX: IERR=',IERR
      CALL SPL2D(THGMX,SIGMX,HJTST,PSITX,PSISX,PSISTX,UHJTST,
     &           NTGPM,NTGMAX+2,NSGMAX+2,4,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX SPL2D for PSIX: IERR=',IERR
      RETURN
      END
C
C   *******************************************
C   ***** Calculate Psi at (sigma, theta) *****
C   *******************************************
C
      FUNCTION PSIF(RSIG,RTHG)
C
      USE libspl2d
      INCLUDE '../eq/eqcomc.inc'
C
      CALL SPL2DF(RTHG,RSIG,PSIL,THGMX,SIGMX,UPSIST,
     &            NTGPM,NTGMAX+2,NSGMAX+2,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX PSIF: SPL2DF ERROR : IERR=',IERR
      PSIF=PSIL
      RETURN
      END
C
C   *******************************************
C   ***** Calculate Hjt at (sigma, theta) *****
C   *******************************************
C
      FUNCTION HJTF(RSIG,RTHG)
C
      USE libspl2d
      INCLUDE '../eq/eqcomc.inc'
C
      CALL SPL2DF(RTHG,RSIG,HJTL,THGMX,SIGMX,UHJTST,
     &            NTGPM,NTGMAX+2,NSGMAX+2,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX HJTF: SPL2DF ERROR : IERR=',IERR
      HJTF=HJTL
      RETURN
      END
C
C   ************************************************
C   **      CALCULATE pp,tt,temp,omega,rho        **
C   ************************************************
C
      SUBROUTINE EQCALP
C
      INCLUDE '../eq/eqcomc.inc'
C
      IMDLEQF=MOD(MDLEQF,5)
      DPS=PSIPA/(NPSMAX-1)
      DO NPS=1,NPSMAX
         PSIPS(NPS)=DPS*(NPS-1)
         PSIPNL=PSIPS(NPS)/PSIPA
C
         IF (IMDLEQF.EQ.0) THEN
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQJPSI(PSIPNL,HJPSID,HJPSI)
            CALL EQTPSI(PSIPNL,TPSI,DTPSI)
            CALL EQOPSI(PSIPNL,OMGPSI,DOMGPSI)
            PPPS(NPS)=PPSI
            OMPS(NPS)=OMGPSI
            TTPS(NPS)=SQRT((2.D0*PI*BB*RR)**2
     &                  +2.D0*RMU0*RRC
     &                  *(2.D0*PI*TJ*HJPSID-RRC*PPSI
     &                  *EXP(RRC**2*OMGPSI**2*AMP/(2.D0*TPSI))))
            TEPS(NPS)=TPSI/(AEE*1.D3)
         ELSEIF (IMDLEQF.EQ.1) THEN
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQFPSI(PSIPNL,FPSI,DFPSI)
            PPPS(NPS)=PPSI
            TTPS(NPS)=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
            OMPS(NPS)=0.D0
            TEPS(NPS)=0.D0
         ELSEIF (IMDLEQF.EQ.2) THEN
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQFIPV(PSIPNL,FPSI,DFPSI)
            PPPS(NPS)=PPSI
            TTPS(NPS)=2.D0*PI*BB*RR+TJ*(FPSI-2.D0*PI*BB*RR)
            OMPS(NPS)=0.D0
            TEPS(NPS)=0.D0
         ELSEIF (IMDLEQF.EQ.3) THEN
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQFIPV(PSIPNL,FPSI,DFPSI)
            PPPS(NPS)=PPSI
            TTPS(NPS)=FPSI
            OMPS(NPS)=0.D0
            TEPS(NPS)=0.D0
         ELSEIF (IMDLEQF.EQ.4) THEN
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQFIPV(PSIPNL,FPSI,DFPSI)
            PPPS(NPS)=PPSI
            TTPS(NPS)=FPSI
            OMPS(NPS)=0.D0
            TEPS(NPS)=0.D0
         ENDIF
      ENDDO
      RETURN
      END
