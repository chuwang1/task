! trpnf.f90

MODULE trpnf

  PRIVATE
  PUBLIC tr_prep_pnf
  PUBLIC tr_pnf
  PRIVATE tr_nf_dt
  PRIVATE tr_nf_dd1
  PRIVATE tr_nf_dd2
  PRIVATE tr_nf_dd3
  PRIVATE tr_nf_dhe31
  PRIVATE tr_nf_dhe32
  PRIVATE tr_nf_tt
  PRIVATE tr_nf_the31
  PRIVATE tr_nf_the32
  PRIVATE tr_nf_the33
  PRIVATE tr_nf_the34
  PRIVATE tr_nf_the35
  PRIVATE tr_nf_the36

CONTAINS

  SUBROUTINE tr_prep_pnf

    USE trcomm
    USE libnf
    IMPLICIT NONE
    INTEGER:: nnf,id_nf

    ! --- initialize fusion cross section and reaction rate ---
    
    IF(model_pnf.EQ.0) RETURN

    DO nnf=1,nnfmax
       id_nf=id_nf_nnf(nnf)
       ns1_nnf(nnf)=ns1_idnf(id_nf)
       ns2_nnf(nnf)=ns2_idnf(id_nf)
       nsp_nnf(nnf)=nsp_idnf(id_nf)
       wgt_nnf(nnf)=wgt_idnf(id_nf)
       eng_nnf(nnf)=eng_idnf(id_nf)
       enn_nnf(nnf)=enn_idnf(id_nf)
    END DO
    
    WRITE(6,*) 'nnf: id_nf,ns1,ns2,nsp,wgt,eng,enn'
    DO nnf=1,nnfmax
       WRITE(6,'(5I4,3ES12.4)') &
            nnf,id_nf_nnf(nnf),ns1_nnf(nnf),ns2_nnf(nnf),nsp_nnf(nnf), &
            wgt_nnf(nnf),eng_nnf(nnf),enn_nnf(nnf)
    END DO

  END SUBROUTINE tr_prep_pnf

  ! *** calculate fusion power ***

  SUBROUTINE tr_pnf

    USE trcomm
    USE libnf
    USE trlib
    IMPLICIT NONE
    REAL(rkind):: ANE,TE,P1,VC3,VCR,WF,VF,TAUS,HYF
    REAL(rkind):: PN1,PN2,PT1,RATE_NF,SNF
    REAL(rkind):: wgt,eng,enn
    INTEGER:: nnf,nr,id_nf,ns1,ns2,nsp,ns

    SNF_NSNNFNR(1:NSMAX,1:NNFMAX,1:NRMAX)=0.D0   ! particle source
    PNF_NSNNFNR(1:NSMAX,1:NNFMAX,1:NRMAX)=0.D0   ! fusion power source
    PNFIN_NNFNR(1:NNFMAX,1:NRMAX)=0.D0           ! fast ion creation
    PNFCL_NSNNFNR(1:NSMAX,1:NNFMAX,1:NRMAX)=0.D0 ! collisional transfer in
    SNFNN_NNFNR(1:NNFMAX,1:NRMAX)=0.D0  ! neutron number
    PNFNN_NNFNR(1:NNFMAX,1:NRMAX)=0.D0  ! neutron power
    
    IF (model_pnf.EQ.0 .AND. ANY(model_nnf > 0)) THEN
        DO nnf=1, NNFM
           SELECT CASE(model_nnf(nnf))
           CASE(0)
              TAUF(nnf,1:NRMAX)=1.D0
           CASE(1:4)
              CALL TRNFDT(nnf)
           CASE(11:14)
              CALL TRNFDD(nnf)
           END SELECT
        END DO
        GOTO 100
    END IF
    DO nnf=1,nnfmax
       id_nf=id_nf_nnf(nnf)
       ns1=ns1_nnf(nnf)
       ns2=ns2_nnf(nnf)
       wgt=wgt_nnf(nnf)
       nsp=nsp_nnf(nnf)
       eng=eng_nnf(nnf)
       enn=enn_nnf(nnf)
       DO NR=1,NRMAX
          PN1=RN(NR,ns1)
          PN2=RN(NR,ns2)
          PT1=RT(NR,ns1)
          RATE_NF=sigmav_nf(id_nf,PT1)
          SNF=wgt*PN1*PN2*1.D20*RATE_NF
          SNF_NSNNFNR(ns1,nnf,nr)=SNF_NSNNFNR(ns1,nnf,nr)-SNF
          SNF_NSNNFNR(ns2,nnf,nr)=SNF_NSNNFNR(ns2,nnf,nr)-SNF
          SNF_NSNNFNR(nsp,nnf,nr)=SNF_NSNNFNR(nsp,nnf,nr)+SNF
          PNF_NSNNFNR(nsp,nnf,nr)=PNF_NSNNFNR(nsp,nnf,nr)+eng*SNF*1.D20
          IF(enn.GT.0.D0) THEN
             SNFNN_NNFNR(nnf,nr)=SNFNN_NNFNR(nnf,nr)+SNF
             PNFNN_NNFNR(nnf,nr)=PNFNN_NNFNR(nnf,nr)+enn*SNF*1.D20
          END IF
       END DO
    END DO

    DO NR=1,NRMAX
       ANE= RN(NR,NS_e)
       TE = RT(NR,NS_e)
       P1   = 3.D0*SQRT(0.5D0*PI)*AME/ANE *(ABS(TE)*RKEV/AME)**1.5D0
       VC3=0.D0
       DO NS=1,NSMAX
          IF(PZ(NS).GT.0.D0) &    ! sum over ions
               VC3=VC3+P1*RN(NR,NS)*PZ(NS)**2/(PA(NS)*AMP)
       END DO
       VCR  = VC3**(1.D0/3.D0)
       DO nnf=1,nnfmax
          nsp=nsp_nnf(nnf)
          WF = RW(NR,NNBMAX+NNF)
          VF =SQRT(2.D0*eng_nnf(nnf)*RKEV/(PA(ns)*AMP))
          HYF=HY(VF/VCR)
          TAUS = 0.2D0*PA(ns)*ABS(TE)**1.5D0 &
               /(PZ(ns)**2*ANE*COULOG(1,ns,ANE,TE))
          TAUF(NNF,NR)= 0.5D0*TAUS*(1.D0-HYF)
       END DO
    END DO
          
100 CONTINUE
          
    ! --- following variables are used in trcalc at every step ---
    
    DO NR=1,NRMAX
       DO NS=1,NSMAX
          SNF_NSNR(NS,NR)=SUM(SNF_NSNNFNR(NS,1:NNFMAX,NR))
          PNF_NSNR(NS,NR)=SUM(PNF_NSNNFNR(NS,1:NNFMAX,NR))
          PNFCL_NSNR(NS,NR)=SUM(PNFCL_NSNNFNR(NS,1:NNFMAX,NR))
       END DO
       PNFIN_NR(NR)=SUM(PNFIN_NNFNR(1:NNFMAX,NR))
    END DO

    RETURN
  END SUBROUTINE tr_pnf

!     ***********************************************************
!           Nuclear reaction (DT)
!     ***********************************************************
  SUBROUTINE TRNFDT(nnf)
    USE TRCOMM
    USE libsigma
    USE trlib
    IMPLICIT NONE
    INTEGER,INTENT(IN):: nnf
    REAL(rkind)   :: &
         ANE, EC, HYF, P1, PTNT, SS, SSB, TAUS, &
         TD, TE, TT, VC3, VCA3, VCD3, VCR, VCT3, VF, WF, ZEFFM, PB
    INTEGER:: NR,NNB,NS_beam

    VF =SQRT(2.D0*3.5D3 *RKEV/AMA) ! alpha velocity

      DO NR=1,NRMAX
         SNF_NSNNFNR(1:NSMAX,NNF,NR)=0.D0
         PNF_NSNNFNR(1:NSMAX,NNF,NR)=0.D0
         ANE= RN(NR,NS_e)
         TE = RT(NR,NS_e)
         TD = RT(NR,NS_D)
         TT = RT(NR,NS_T)
         SS = SIGMAM(TD,TT)
         IF(model_nnf(nnf).GE.3) THEN
            ZEFFM = (PZ(NS_D  )*PZ(NS_D  )*RN(NR,NS_D  )/PA(NS_D  ) &
                    +PZ(NS_T)  *PZ(NS_T  )*RN(NR,NS_T  )/PA(NS_T  ) &
                    +PZ(NS_He4)*PZ(NS_He4)*RN(NR,NS_He4)/PA(NS_He4) &
                    +PZC(NR)*PZC(NR) *ANC(NR) /12.D0 &
                    +PZFE(NR)*PZFE(NR)*ANFE(NR)/52.D0)/ANE
            SSB=0.D0
            PB=0.D0
            DO NNB=1,NNBMAX
               IF(pnbin(nnb).GT.0.D0) THEN
               NS_beam=ns_nnb(nnb)
               IF (NS_beam.EQ.NS_D) THEN
               EC  = 14.8D0*TE*PA(NS_beam)*ZEFFM**(2.D0/3.D0)
               TAUS= 0.2D0*PA(NS_beam)*ABS(TE)**1.5D0 &
                    /(PZ(NS_beam)**2*ANE*COULOG(1,NS_beam,ANE,TE))
               PTNT= PNB_NNBNR(nnb,NR)*TAUS &
                    /(RN(NR,NS_beam)*1.D20*PNBENG(NNB)*RKEV)
               SSB = SSB+SIGMAB(PNBENG(NNB),EC,TT,PTNT)
               PB  = PB +PNB_NNBNR(NNB,NR)
               END IF
               END IF
            END DO
         ELSE
            SSB=0.D0
         ENDIF
         
         SNF_NSNNFNR(NS_He4,NNF,NR) = (SS+SSB)*RN(NR,NS_D)*RN(NR,NS_T)*1.D20 
         PNF_NSNNFNR(NS_He4,NNF,NR) = SNF_NSNNFNR(NS_He4,NNF,NR)*3.5D3*RKEV*1.D20
         IF(MOD(model_nnf(nnf),2).EQ.1) &
              SNF_NSNNFNR(NS_He4,NNF,NR)=0.D0
         SNF_NSNNFNR(NS_D,NNF,NR) =-SNF_NSNNFNR(NS_He4,NNF,NR)
         SNF_NSNNFNR(NS_T,NNF,NR) =-SNF_NSNNFNR(NS_He4,NNF,NR)
         PNF_NSNNFNR(NS_D,NNF,NR) &
              =SNF_NSNNFNR(NS_He4,NNF,NR)*RT(NR,NS_D)*RKEV*1.D20
         PNF_NSNNFNR(NS_T,NNF,NR) &
              =SNF_NSNNFNR(NS_He4,NNF,NR)*RT(NR,NS_T)*RKEV*1.D20
      ENDDO

      DO NR=1,NRMAX
         ANE= RN(NR,NS_e)
         TE = RT(NR,NS_e)
         WF = RW(NR,NNBMAX+NNF)
         P1   = 3.D0*SQRT(0.5D0*PI)*AME/ANE *(ABS(TE)*RKEV/AME)**1.5D0
         VCD3 = P1*RN(NR,NS_D)*PZ(NS_D)**2/AMD
         VCT3 = P1*RN(NR,NS_T)*PZ(NS_T)**2/AMT
         VCA3 = P1*RN(NR,NS_He4)*PZ(NS_He4)**2/AMA
         VC3  = VCD3+VCT3+VCA3
         VCR  = VC3**(1.D0/3.D0)
         HYF=HY(VF/VCR)
         TAUS = 0.2D0*PA(NS_He4)*ABS(TE)**1.5D0 &
              /(PZ(NS_He4)**2*ANE*COULOG(1,2,ANE,TE))
         TAUF(NNF,NR)= 0.5D0*TAUS*(1.D0-HYF)
         RNF(NR,NNBMAX+NNF) &
              = 2.D0*LOG(1.D0+(VF/VCR)**3)*WF /(3.D0*(1.D0-HYF)*3.5D3)
         IF(RNF(NR,NNBMAX+NNF).GT.0.D0) THEN
            RTF(NR,NNBMAX+NNF)= WF/RNF(NR,NNBMAX+NNF)
         ELSE
            RTF(NR,NNBMAX+NNF)= 0.D0
         ENDIF
         PNFIN_NNFNR(NNF,NR) = WF*RKEV*1.D20/TAUF(nnf,NR)
         PNFCL_NSNNFNR(NS_e,  NNF,NR)=    (1.D0-HYF)*PNFIN_NNFNR(NNF,NR)
         PNFCL_NSNNFNR(NS_D,  NNF,NR)=(VCD3/VC3)*HYF*PNFIN_NNFNR(NNF,NR)
         PNFCL_NSNNFNR(NS_T,  NNF,NR)=(VCT3/VC3)*HYF*PNFIN_NNFNR(NNF,NR)
         PNFCL_NSNNFNR(NS_He4,NNF,NR)=(VCA3/VC3)*HYF*PNFIN_NNFNR(NNF,NR)
      ENDDO
      RETURN
  END SUBROUTINE TRNFDT

!     ***********************************************************
!           Nuclear reaction (DD)
!     ***********************************************************
  SUBROUTINE TRNFDD(nnf)
         USE TRCOMM
         USE trlib
         IMPLICIT NONE
         INTEGER,INTENT(IN):: nnf
         REAL(rkind)   :: &
              ANE, EC, HYF, P1, PTNT, SS, SSB, TAUS, &
              TD, TE, TT, VC3, VCA3, VCD3, VCR, VCT3, VF, WF, ZEFFM, PB
         INTEGER:: NR,NNB,NS_beam
     
         VF =SQRT(2.D0*3.5D3 *RKEV/AMA) ! alpha velocity
   
         DO NR=1,NRMAX
             SNF_NSNR(1:NSMAX,NR)=0.D0
             PNF_NSNR(1:NSMAX,NR)=0.D0
             ANE= RN(NR,NS_e)
             TE = RT(NR,NS_e)
             TD = RT(NR,NS_D)
             TT = RT(NR,NS_T)
             SS = SIGMAM_DD(TD)
             IF(model_nnf(nnf).GE.13) THEN
               ZEFFM = (PZ(NS_D  )*PZ(NS_D  )*RN(NR,NS_D  )/PA(NS_D  ) &
                       +PZ(NS_T  )*PZ(NS_T  )*RN(NR,NS_T  )/PA(NS_T  ) &
                       +PZ(NS_He4)*PZ(NS_He4)*RN(NR,NS_He4)/PA(NS_He4) &
                       +PZC(NR) *PZC(NR) *ANC(NR) /12.D0 &
                       +PZFE(NR)*PZFE(NR)*ANFE(NR)/52.D0)/ANE
                SSB=0.D0
                PB=0.D0
                DO NNB=1,NNBMAX
                   NS_beam=ns_nnb(nnb)
                   IF(pnbin(nnb).GT.0.D0) THEN
                   EC  = 14.8D0*TE*PA(NS_beam)*ZEFFM**(2.D0/3.D0)
                   TAUS= 0.2D0*PA(NS_beam)*ABS(TE)**1.5D0 &
                        /(PZ(NS_beam)**2*ANE*COULOG(1,NS_beam,ANE,TE))
                   PTNT= PNB_NNBNR(nnb,NR)*TAUS &
                        /(RN(NR,NS_beam)*1.D20*PNBENG(NNB)*RKEV)
                   SSB = SSB+PNB_NNBNR(NNB,NR)*SIGMAB(PNBENG(NNB),EC,TD,PTNT)
                   PB  = PB +PNB_NNBNR(NNB,NR)
                   END IF
                END DO
                IF(PB.NE.0.D0) THEN
                   SSB=SSB/PB
                ELSE
                   SSB=0.D0
                END IF
             ELSE
                SSB=0.D0
             ENDIF
             SNF_NSNNFNR(NS_T,NNF,NR) = 0.5*0.5*(SS+SSB)*RN(NR,2)*RN(NR,2)*1.D20
             PNF_NSNNFNR(NS_T,NNF,NR) = 2*SNF_NSNNFNR(NS_T,NNF,NR)*((1.01D3+3.03D3+0.82D3)/2)*RKEV*1.D20
             IF(MOD(model_nnf(nnf),2).EQ.1) SNF_NSNNFNR(NS_He4,NNF,NR)=0.D0
             SNF_NSNNFNR(NS_D,NNF,NR) =-SNF_NSNNFNR(NS_T,NNF,NR)
          ENDDO
 
         DO NR=1,NRMAX
            ANE= RN(NR,1)
            TE = RT(NR,1)
            WF = RW(NR,NNBMAX+NNF)
            P1   = 3.D0*SQRT(0.5D0*PI)*AME/ANE *(ABS(TE)*RKEV/AME)**1.5D0
            VCD3 = P1*RN(NR,2)*PZ(NS_D)**2/AMD
            VCT3 = P1*RN(NR,3)*PZ(NS_T)**2/AMT
            VCA3 = P1*RN(NR,4)*PZ(NS_He4)**2/AMA
            VC3  = VCD3+VCT3+VCA3
            VCR  = VC3**(1.D0/3.D0)
            HYF=HY(VF/VCR)
            TAUS = 0.2D0*PA(NS_He4)*ABS(TE)**1.5D0 /(PZ(NS_He4)**2*ANE*COULOG(1,2,ANE,TE))
            TAUF(nnf,NR) = 0.5D0*TAUS*(1.D0-HYF)
            RNF(NR,NNBMAX+NNF)= 2.D0*LOG(1.D0+(VF/VCR)**3)*WF /(3.D0*(1.D0-HYF)*3.5D3)
            IF(RNF(NR,NNBMAX+NNF).GT.0.D0) THEN
               RTF(NR,NNBMAX+NNF)= WF/RNF(NR,NNBMAX+NNF)
            ELSE
               RTF(NR,NNBMAX+NNF)= 0.D0
            ENDIF
            PNFIN_NNFNR(NNF,NR) = WF*RKEV*1.D20/TAUF(nnf,NR)
            PNFCL_NSNNFNR(NS_e,  NNF,NR)=    (1.D0-HYF)*PNFIN_NNFNR(NNF,NR)
            PNFCL_NSNNFNR(NS_D,  NNF,NR)=(VCD3/VC3)*HYF*PNFIN_NNFNR(NNF,NR)
            PNFCL_NSNNFNR(NS_T,  NNF,NR)=(VCT3/VC3)*HYF*PNFIN_NNFNR(NNF,NR)
            PNFCL_NSNNFNR(NS_He4,NNF,NR)=(VCA3/VC3)*HYF*PNFIN_NNFNR(NNF,NR)
         ENDDO
         RETURN
  END SUBROUTINE TRNFDD

      FUNCTION SIGMAM(TD,TT)
         USE bpsd_kinds
         USE libsigma
         IMPLICIT NONE
         REAL(rkind),INTENT(IN)::  TD,TT
         REAL(rkind):: TI,SIGMAM
         TI = (3.D0*ABS(TD)+2.D0*ABS(TT))/5.D0
         SIGMAM=sigmavm_dt(TI)*(1E-6)
         RETURN
      END FUNCTION SIGMAM

      FUNCTION SIGMAM_DD(TD)
          USE trcomm,ONLY: rkind
          IMPLICIT NONE
          REAL(rkind) TD,SIGMAM_DD
          REAL(rkind) TI
          TI=TD
          IF(TI.GE.4) THEN 
             SIGMAM_DD = ((2.1069E-20)+ (-3.7748E-20)*TI + (1.1242E-20)*TI**2 + (-1.9831E-22)*TI**3 + (1.3421E-24)*TI**4)*1E-6
          ELSE
             SIGMAM_DD = (4.428*1E-20)*1E-6
          ENDIF    
          RETURN
      END FUNCTION SIGMAM_DD

      FUNCTION SIGMAB(EB,EC,TI,PTNT)
      USE trcomm,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind) EB,EC,TI,PTNT,SIGMAB
      REAL(rkind) XB,XC,AG1,AG2,AG3,AL1,AL2,AL3,X1,X2,X3,X4,SA
      XB=SQRT(EB/127.D0)
      XC=SQRT(EC/127.D0)
      AG1= 1.06D0-0.058D0*SQRT(ABS(TI))
      AG2= 1.06D0-0.058D0*SQRT(ABS(TI))
      AG3= 0.33D0
      AL1= 1.D0/(0.40D0+0.032D0*TI)
      AL2=-1.D0/(0.91D0+0.016D0*SQRT(ABS(TI)))
      AL3=-0.11D0
      X1=0.97D0-AG1/AL1
      X2=0.97D0
      X3=0.97D0+(AG2-AG3)/(AL3-AL2)
      X4=0.97D0+3.D0
      IF(XB.LT.X1) THEN
         SA=0.D0
      ELSE
         SA=-SIGMBS(X1,AG1,AL1,XC)
         IF(XB.LT.X2) THEN
            SA=SA+SIGMBS(XB,AG1,AL1,XC)
         ELSE
            SA=SA+SIGMBS(X2,AG1,AL1,XC)-SIGMBS(X2,AG2,AL2,XC)
            IF(XB.LT.X3) THEN
               SA=SA+SIGMBS(XB,AG2,AL2,XC)
            ELSE
               SA=SA+SIGMBS(X3,AG2,AL2,XC)-SIGMBS(X3,AG3,AL3,XC)
               IF(XB.LT.X4) THEN
                  SA=SA+SIGMBS(XB,AG3,AL3,XC)
               ELSE
                  SA=SA+SIGMBS(X4,AG3,AL3,XC)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      SIGMAB=PTNT*1.67D-21*SA
      RETURN
      END FUNCTION SIGMAB

      FUNCTION SIGMBS(XX,RGG,RGL,XC)
      USE trcomm,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind) :: XX,RGG,RGL,XC,SIGMBS
      REAL(rkind) :: X
      X=XX/XC
      SIGMBS=((RGG-0.97D0*RGL)/3.D0+XC*RGL/6.D0)*LOG(X*X*X+1.D0) &
           +XC*RGL*(X-LOG(X+1.D0)/2.D0 &
                   -ATAN((2.D0*X-1.D0)/SQRT(3.D0))/SQRT(3.D0))
      RETURN
    END FUNCTION SIGMBS

END MODULE trpnf
