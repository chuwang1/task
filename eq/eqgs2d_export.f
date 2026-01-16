C     $Id$
C
C     ***** EXPORT EQGS2D DATA TO CSV *****
C
      SUBROUTINE EQGS2D_EXPORT
C
C     This subroutine exports 2D equilibrium data to CSV files:
C     1. PSIRZ(R,Z)    - Poloidal flux on R-Z grid
C     2. DPSIDR(R,Z)   - Radial derivative of poloidal flux
C     3. DPSIDZ(R,Z)   - Vertical derivative of poloidal flux
C     4. RPS(theta,psi) - R coordinate of flux surfaces
C     5. ZPS(theta,psi) - Z coordinate of flux surfaces
C     6. Separatrix boundary points (RSU, ZSU, RSW, ZSW)
C     7. Key parameters (axis position, etc.)
C
      INCLUDE 'eqcomq.inc'
C
      DIMENSION GPSIRZ(NRGM,NZGM),GDPSIDR(NRGM,NZGM),GDPSIDZ(NRGM,NZGM)
      CHARACTER*80 FNAME
      REAL*8 DPSIDR,DPSIDZ
C
      WRITE(6,*) '# EXPORTING EQGS2D DATA TO CSV FILES...'
C
C     Setup PSIG interpolation
      CALL setup_psig
C
C     === 1. Export PSIRZ(R,Z) grid data ===
      FNAME = 'eqgs2d_01_PSIRZ_grid.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')

C     Write header with R grid values
      WRITE(31,'(A)',ADVANCE='NO') 'Z_m / R_m'
      DO NR=1,NRGMAX
         WRITE(31,'(A,E16.8)',ADVANCE='NO') ',',GUCLIP(RG(NR))
      ENDDO
      WRITE(31,*)  ! New line

C     Write data: each row is Z, followed by PSIRZ values
      DO NZ=1,NZGMAX
         WRITE(31,'(E16.8)',ADVANCE='NO') GUCLIP(ZG(NZ))
         DO NR=1,NRGMAX
            GPSIRZ(NR,NZ)=GUCLIP(PSIRZ(NR,NZ))
            WRITE(31,'(A,E16.8)',ADVANCE='NO') ',',GPSIRZ(NR,NZ)
         ENDDO
         WRITE(31,*)  ! New line
      ENDDO
      CLOSE(31)
      WRITE(6,'(A)') '  [1] PSIRZ(R,Z) grid saved'
C
C     === 2. Export DPSIDR(R,Z) grid data ===
      FNAME = 'eqgs2d_02_DPSIDR_grid.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')

      WRITE(31,'(A)',ADVANCE='NO') 'Z_m / R_m'
      DO NR=1,NRGMAX
         WRITE(31,'(A,E16.8)',ADVANCE='NO') ',',GUCLIP(RG(NR))
      ENDDO
      WRITE(31,*)

      DO NZ=1,NZGMAX
         WRITE(31,'(E16.8)',ADVANCE='NO') GUCLIP(ZG(NZ))
         DO NR=1,NRGMAX
            CALL PSIGD(RG(NR),ZG(NZ),DPSIDR,DPSIDZ)
            GDPSIDR(NR,NZ)=DPSIDR
            WRITE(31,'(A,E16.8)',ADVANCE='NO') ',',GDPSIDR(NR,NZ)
         ENDDO
         WRITE(31,*)
      ENDDO
      CLOSE(31)
      WRITE(6,'(A)') '  [2] DPSIDR(R,Z) grid saved'
C
C     === 3. Export DPSIDZ(R,Z) grid data ===
      FNAME = 'eqgs2d_03_DPSIDZ_grid.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')

      WRITE(31,'(A)',ADVANCE='NO') 'Z_m / R_m'
      DO NR=1,NRGMAX
         WRITE(31,'(A,E16.8)',ADVANCE='NO') ',',GUCLIP(RG(NR))
      ENDDO
      WRITE(31,*)

      DO NZ=1,NZGMAX
         WRITE(31,'(E16.8)',ADVANCE='NO') GUCLIP(ZG(NZ))
         DO NR=1,NRGMAX
            CALL PSIGD(RG(NR),ZG(NZ),DPSIDR,DPSIDZ)
            GDPSIDZ(NR,NZ)=DPSIDZ
            WRITE(31,'(A,E16.8)',ADVANCE='NO') ',',GDPSIDZ(NR,NZ)
         ENDDO
         WRITE(31,*)
      ENDDO
      CLOSE(31)
      WRITE(6,'(A)') '  [3] DPSIDZ(R,Z) grid saved'
C
C     === 4. Export flux surfaces R(theta,psi) ===
      FNAME = 'eqgs2d_04_RPS_flux_surfaces.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')

C     Header: theta / psi index
      WRITE(31,'(A)',ADVANCE='NO') 'theta_rad / psi_index'
      DO NR=1,NRMAX
         WRITE(31,'(A,I4)',ADVANCE='NO') ',',NR
      ENDDO
      WRITE(31,*)

C     Data: each row is theta, followed by R values at different psi
      DO NTH=1,NTHMAX
         THETA = 2.0*PI*DBLE(NTH-1)/DBLE(NTHMAX)
         WRITE(31,'(E16.8)',ADVANCE='NO') THETA
         DO NR=1,NRMAX
            WRITE(31,'(A,E16.8)',ADVANCE='NO') ',',GUCLIP(RPS(NTH,NR))
         ENDDO
         WRITE(31,*)
      ENDDO
      CLOSE(31)
      WRITE(6,'(A)') '  [4] RPS(theta,psi) flux surfaces saved'
C
C     === 5. Export flux surfaces Z(theta,psi) ===
      FNAME = 'eqgs2d_05_ZPS_flux_surfaces.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')

      WRITE(31,'(A)',ADVANCE='NO') 'theta_rad / psi_index'
      DO NR=1,NRMAX
         WRITE(31,'(A,I4)',ADVANCE='NO') ',',NR
      ENDDO
      WRITE(31,*)

      DO NTH=1,NTHMAX
         THETA = 2.0*PI*DBLE(NTH-1)/DBLE(NTHMAX)
         WRITE(31,'(E16.8)',ADVANCE='NO') THETA
         DO NR=1,NRMAX
            WRITE(31,'(A,E16.8)',ADVANCE='NO') ',',GUCLIP(ZPS(NTH,NR))
         ENDDO
         WRITE(31,*)
      ENDDO
      CLOSE(31)
      WRITE(6,'(A)') '  [5] ZPS(theta,psi) flux surfaces saved'
C
C     === 6. Export individual flux surface contours ===
      FNAME = 'eqgs2d_06_flux_surface_contours.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'psi_index,theta_index,theta_rad,R_m,Z_m,PSIP'

      DO NR=1,NRMAX
         DO NTH=1,NTHMAX
            THETA = 2.0*PI*DBLE(NTH-1)/DBLE(NTHMAX)
            WRITE(31,'(I4,A,I4,4(A,E16.8))')
     &         NR,',',NTH,',',THETA,',',
     &         GUCLIP(RPS(NTH,NR)),',',GUCLIP(ZPS(NTH,NR)),',',
     &         GUCLIP(PSIP(NR))
         ENDDO
      ENDDO
      CLOSE(31)
      WRITE(6,'(A)') '  [6] Flux surface contours saved'
C
C     === 7. Export separatrix boundary (if exists) ===
      IF(NSUMAX.GT.0) THEN
         FNAME = 'eqgs2d_07_separatrix.csv'
         OPEN(31,FILE=FNAME,STATUS='REPLACE')
         WRITE(31,'(A)') 'point_index,RSU_m,ZSU_m,RSW_m,ZSW_m'

         DO NSU=1,NSUMAX
            WRITE(31,'(I4,4(A,E16.8))') NSU,',',
     &         GUCLIP(RSU(NSU)),',',GUCLIP(ZSU(NSU)),',',
     &         GUCLIP(RSW(NSU)),',',GUCLIP(ZSW(NSU))
         ENDDO
         CLOSE(31)
         WRITE(6,'(A)') '  [7] Separatrix boundary saved'
      ELSE
         WRITE(6,'(A)') '  [7] No separatrix data (NSUMAX=0)'
      ENDIF
C
C     === 8. Export key equilibrium parameters ===
      FNAME = 'eqgs2d_08_parameters.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,*) 'parameter,value,unit,description'

      WRITE(31,'(A,E16.8,A)') 'RAXIS,',RAXIS,',m,axis R'
      WRITE(31,'(A,E16.8,A)') 'ZAXIS,',ZAXIS,',m,axis Z'
      WRITE(31,'(A,E16.8,A)') 'PSITA,',PSITA,',Wb,tor flux'
      WRITE(31,'(A,E16.8,A)') 'PSIPA,',PSIPA,',Wb,pol flux'
      WRITE(31,'(A,E16.8,A)') 'PSI0,',PSI0,',Wb,boundary'
      WRITE(31,'(A,I6,A)') 'NRGMAX,',NRGMAX,',-,R pts'
      WRITE(31,'(A,I6,A)') 'NZGMAX,',NZGMAX,',-,Z pts'
      WRITE(31,'(A,I6,A)') 'NTHMAX,',NTHMAX,',-,theta pts'
      WRITE(31,'(A,I6,A)') 'NRMAX,',NRMAX,',-,radial surfs'
      WRITE(31,'(A,I6,A)') 'NRPMAX,',NRPMAX,',-,inside sep'
      WRITE(31,'(A,I6,A)') 'NSUMAX,',NSUMAX,',-,sep pts'
      CLOSE(31)
      WRITE(6,'(A)') '  [8] Key parameters saved'
C
C     === 9. Export R and Z grid vectors ===
      FNAME = 'eqgs2d_09_RZ_grids.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'index,R_m,Z_m'

      DO NR=1,MAX(NRGMAX,NZGMAX)
         IF(NR.LE.NRGMAX .AND. NR.LE.NZGMAX) THEN
            WRITE(31,'(I4,2(A,E16.8))') NR,',',
     &         GUCLIP(RG(NR)),',',GUCLIP(ZG(NR))
         ELSEIF(NR.LE.NRGMAX) THEN
            WRITE(31,'(I4,A,E16.8,A)') NR,',',GUCLIP(RG(NR)),','
         ELSEIF(NR.LE.NZGMAX) THEN
            WRITE(31,'(I4,A,A,E16.8)') NR,',,',GUCLIP(ZG(NR))
         ENDIF
      ENDDO
      CLOSE(31)
      WRITE(6,'(A)') '  [9] R and Z grid vectors saved'
C
C     === 10. Export PSIP profile (for reference) ===
      FNAME = 'eqgs2d_10_PSIP_profile.csv'
      OPEN(31,FILE=FNAME,STATUS='REPLACE')
      WRITE(31,'(A)') 'psi_index,PSIP_normalized'

      DO NR=1,NRMAX
         WRITE(31,'(I4,A,E16.8)') NR,',',GUCLIP(PSIP(NR))
      ENDDO
      CLOSE(31)
      WRITE(6,'(A)') '  [10] PSIP profile saved'
C
      WRITE(6,*) '# ALL 2D EQUILIBRIUM DATA EXPORTED TO CSV FILES.'
C
      RETURN
      END
