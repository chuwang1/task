C     $Id$
C
C     ***** EXPORT EQGS1D DATA TO CSV *****
C
      SUBROUTINE EQGS1D_EXPORT
C
C     This subroutine exports all 23 plots data from EQGS1D to CSV files
C
      USE libspl1d
      INCLUDE 'eqcomq.inc'
C
      DIMENSION GX(NRM),GY(NRM,6)
      CHARACTER*80 FNAME
      INTEGER IPLOT,IERRL
C
      WRITE(6,*) '# EXPORTING EQGS1D DATA TO CSV FILES...'
C
C     Initialize X-axis (PSIP)
      DO NR=1,NRPMAX
         GX(NR)=GUCLIP(PSIP(NR))
      ENDDO
C
C     ===== PLOT 1: PPS (Pressure) =====
      IPLOT = 1
      FNAME = 'eqgs1d_01_PPS.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,PPS_MPa'
      DO NR=1,NRPMAX
         GY(NR,1)=GUCLIP(PPS(NR))*1.E-6
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GY(NR,1)
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] PPS saved'
C
C     ===== PLOT 2: TTS (Toroidal Flux) =====
      IPLOT = 2
      FNAME = 'eqgs1d_02_TTS.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,TTS_Wb'
      DO NR=1,NRPMAX
         GY(NR,1)=GUCLIP(TTS(NR))
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GY(NR,1)
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] TTS saved'
C
C     ===== PLOT 3: RLEN (Magnetic field line length) =====
      IPLOT = 3
      FNAME = 'eqgs1d_03_RLEN.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,RLEN_m'
      DO NR=1,NRPMAX
         GY(NR,1)=GUCLIP(RLEN(NR))
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GY(NR,1)
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] RLEN saved'
C
C     ===== PLOT 4: QPS (Safety factor) =====
      IPLOT = 4
      FNAME = 'eqgs1d_04_QPS.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,QPS'
      DO NR=1,NRPMAX
         GY(NR,1)=GUCLIP(QPS(NR))
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GY(NR,1)
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] QPS saved'
C
C     ===== PLOT 5: VPS (Volume) =====
      IPLOT = 5
      FNAME = 'eqgs1d_05_VPS.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,VPS_m3'
      DO NR=1,NRPMAX
         GY(NR,1)=GUCLIP(VPS(NR))
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GY(NR,1)
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] VPS saved'
C
C     ===== PLOT 6: SPS (Surface area) =====
      IPLOT = 6
      FNAME = 'eqgs1d_06_SPS.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,SPS_m2'
      DO NR=1,NRPMAX
         GY(NR,1)=GUCLIP(SPS(NR))
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GY(NR,1)
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] SPS saved'
C
C     ===== PLOT 7: RRMIN/MAX, ZZMIN/MAX, BBMIN/MAX =====
      IPLOT = 7
      FNAME = 'eqgs1d_07_MINMAX.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)')
     &  'PSIP,RRMIN_m,RRMAX_m,ZZMIN_m,ZZMAX_m,BBMIN_T,BBMAX_T'
      DO NR=1,NRPMAX
         WRITE(31,'(E16.8,6(A,E16.8))') GX(NR),',',
     &      GUCLIP(RRMIN(NR)),',',GUCLIP(RRMAX(NR)),',',
     &      GUCLIP(ZZMIN(NR)),',',GUCLIP(ZZMAX(NR)),',',
     &      GUCLIP(BBMIN(NR)),',',GUCLIP(BBMAX(NR))
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] MIN/MAX saved'
C
C     ===== PLOT 8: BBMIN/MAX =====
      IPLOT = 8
      FNAME = 'eqgs1d_08_BB.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,BBMIN_T,BBMAX_T'
      DO NR=1,NRPMAX
         WRITE(31,'(E16.8,2(A,E16.8))') GX(NR),',',
     &      GUCLIP(BBMIN(NR)),',',GUCLIP(BBMAX(NR))
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] BB saved'
C
C     ===== PLOT 9-16: Average values =====
      IPLOT = 9
      FNAME = 'eqgs1d_09_AVERR2.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,AVERR2_m2'
      DO NR=1,NRPMAX
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GUCLIP(AVERR2(NR))
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] AVERR2 saved'
C
      IPLOT = 10
      FNAME = 'eqgs1d_10_AVEIR2.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,AVEIR2_per_m2'
      DO NR=1,NRPMAX
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GUCLIP(AVEIR2(NR))
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] AVEIR2 saved'
C
      IPLOT = 11
      FNAME = 'eqgs1d_11_AVEBB2.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,AVEBB2_T2'
      DO NR=1,NRPMAX
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GUCLIP(AVEBB2(NR))
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] AVEBB2 saved'
C
      IPLOT = 12
      FNAME = 'eqgs1d_12_AVEIB2.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,AVEIB2_per_T2'
      DO NR=1,NRPMAX
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GUCLIP(AVEIB2(NR))
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] AVEIB2 saved'
C
      IPLOT = 13
      FNAME = 'eqgs1d_13_AVEGV2.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,AVEGV2_m4'
      DO NR=1,NRPMAX
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GUCLIP(AVEGV2(NR))
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] AVEGV2 saved'
C
      IPLOT = 14
      FNAME = 'eqgs1d_14_AVEGVR2.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,AVEGVR2_m2'
      DO NR=1,NRPMAX
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GUCLIP(AVEGVR2(NR))
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] AVEGVR2 saved'
C
      IPLOT = 15
      FNAME = 'eqgs1d_15_AVEGP2.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,AVEGP2'
      DO NR=1,NRPMAX
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GUCLIP(AVEGP2(NR))
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] AVEGP2 saved'
C
      IPLOT = 16
      FNAME = 'eqgs1d_16_AVEJ.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,AVEJPR_MAm2,AVEJTR_MAm2'
      DO NR=1,NRPMAX
         WRITE(31,'(E16.8,2(A,E16.8))') GX(NR),',',
     &      GUCLIP(AVEJPR(NR))*1.E-6,',',GUCLIP(AVEJTR(NR))*1.E-6
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] AVEJ saved'
C
C     ===== PLOT 17-20: Geometry parameters =====
      IPLOT = 17
      FNAME = 'eqgs1d_17_RRPSI.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,RRPSI_m'
      DO NR=1,NRPMAX
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GUCLIP(RRPSI(NR))
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] RRPSI saved'
C
      IPLOT = 18
      FNAME = 'eqgs1d_18_RSPSI.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,RSPSI_m'
      DO NR=1,NRPMAX
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GUCLIP(RSPSI(NR))
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] RSPSI saved'
C
      IPLOT = 19
      FNAME = 'eqgs1d_19_ELIPPSI.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,ELIPPSI'
      DO NR=1,NRPMAX
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GUCLIP(ELIPPSI(NR))
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] ELIPPSI saved'
C
      IPLOT = 20
      FNAME = 'eqgs1d_20_TRIGPSI.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,TRIGPSI'
      DO NR=1,NRPMAX
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GUCLIP(TRIGPSI(NR))
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] TRIGPSI saved'
C
C     ===== PLOT 21-23: Derivatives =====
      IPLOT = 21
      FNAME = 'eqgs1d_21_DVDPSIP.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,DVDPSIP_m3perWb'
      DO NR=1,NRPMAX
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GUCLIP(DVDPSIP(NR))
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] DVDPSIP saved'
C
      IPLOT = 22
      FNAME = 'eqgs1d_22_DVDPSIT.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,DVDPSIT_m3perWb'
      DO NR=1,NRPMAX
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GUCLIP(DVDPSIT(NR))
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] DVDPSIT saved'
C
      IPLOT = 23
      FNAME = 'eqgs1d_23_AVEGV.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIP,AVEGV_m2'
      DO NR=1,NRPMAX
         WRITE(31,'(E16.8,A,E16.8)') GX(NR),',',GUCLIP(AVEGV(NR))
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] AVEGV saved'
C
C     ===== PLOT 24: Exact EQIPQP inputs on NRVMAX grid =====
      IPLOT = 24
      IF(IEQSNAP.GE.2) THEN
         WRITE(6,*) '  [24] use stored EQIPQP snapshot from eqdata'
      ELSEIF(IEQSNAP.EQ.1) THEN
         WRITE(6,*) '  [24] use stored EQCALV snapshot from eqdata'
      ELSE
         IERRL = 0
         CALL EQCALV(IERRL)
         IF(IERRL.NE.0) THEN
            WRITE(6,'(A,I5)')
     &         'XX EQGS1D_EXPORT: EQCALV failed: IERR=',IERRL
         ENDIF
         CALL EQIPQP
         IEQSNAP=2
         WRITE(6,*) '  [24] snapshot missing; recomputed from PSIRZ'
      ENDIF
C
      FNAME = 'eqgs1d_24_EQIPQP_INPUTS.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'PSIPV,PSIPNV,DPPSI,QPSI,RSV,AVIR2,VPV,AVRR2,'
     &   //'FIPV'
      DO NRV=1,NRVMAX
         PSIPNL=PSIPNV(NRV)
         IF(IEQSNAP.GE.2) THEN
            DPPSI=DPIPV(NRV)
            QPSI=QIPV(NRV)
         ELSE
            CALL EQPPSI(PSIPNL,PPSI,DPPSI)
            CALL EQQPSI(PSIPNL,QPSI)
         ENDIF
         WRITE(31,'(1P,E24.16,8(A1,E24.16))') PSIPV(NRV),',',
     &      PSIPNV(NRV),',',DPPSI,',',QPSI,',',RSV(NRV),',',
     &      AVIR2(NRV),',',VPV(NRV),',',AVRR2(NRV),',',FIPV(NRV)
      ENDDO
      CLOSE(31)
      WRITE(6,'(A,I2,A)') '  [',IPLOT,'] EQIPQP inputs saved'
C
      WRITE(6,*) '# ALL 24 PLOTS EXPORTED TO CSV FILES.'
C
      RETURN
      END
