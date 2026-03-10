! trfixed.f90
! Module for reading external transport coefficient profiles

MODULE trfixed

  USE trcomm,ONLY: rkind

  PUBLIC tr_prep_chifixed
  PUBLIC tr_prep_prlfixed

CONTAINS

  ! *** read chi (thermal diffusivity) from external file ***
  ! File format: space-separated with header line
  ! Columns: r[m]  ne  ni  chi_e  chi_i  Se  Si  ...
  ! chi_e, chi_i are in m^2/s

  SUBROUTINE tr_prep_chifixed
    USE trcomm
    USE libfio
    IMPLICIT NONE
    INTEGER:: nfl,nr,ndata_csv,i,j,ierr,ios
    INTEGER,PARAMETER:: NMAX_CSV=300
    REAL(rkind):: r_csv(NMAX_CSV),rho_csv(NMAX_CSV),chi_e_csv(NMAX_CSV),chi_i_csv(NMAX_CSV)
    REAL(rkind):: dummy_ne,dummy_ni,rho,frac
    CHARACTER(LEN=512):: line

    ! Allow reading chi file when model_chifixed>=1 OR MDLKAI=170:189
    IF(model_chifixed.EQ.0 .AND. (MDLKAI.LT.170 .OR. MDLKAI.GT.189)) RETURN

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
          frac = 0.5D0 * (1.D0 - TANH((rho - 0.4D0) / 0.1D0))
          AKEXT_E(nr) = AKEXT_E(nr) * (1.D0 - 0.5D0 * frac)
          AKEXT_I(nr) = AKEXT_I(nr) * (1.D0 - 0.5D0 * frac)
       END DO
       WRITE(6,'(A)') '## tr_prep_chifixed: mode=1, core chi reduced to 0.5x for r/a < 0.4'

    ELSE IF(model_chifixed.EQ.2) THEN
       ! model_chifixed=2: Multiply core chi by chifixed_factor (r/a < 0.4)
       DO nr=1,NRMAX
          rho = rm(nr)
          frac = 0.5D0 * (1.D0 - TANH((rho - 0.4D0) / 0.1D0))
          AKEXT_E(nr) = AKEXT_E(nr) * (1.D0 + (chifixed_factor - 1.D0) * frac)
          AKEXT_I(nr) = AKEXT_I(nr) * (1.D0 + (chifixed_factor - 1.D0) * frac)
       END DO
       WRITE(6,'(A,F8.3,A)') '## tr_prep_chifixed: mode=2, core chi multiplied by ',chifixed_factor,' for r/a < 0.4'

    ELSE IF(model_chifixed.EQ.3) THEN
       ! model_chifixed=3: Multiply edge chi by chifixed_factor (r/a > 0.5)
       DO nr=1,NRMAX
          rho = rm(nr)
          frac = 0.5D0 * (1.D0 + TANH((rho - 0.5D0) / 0.1D0))
          AKEXT_E(nr) = AKEXT_E(nr) * (1.D0 + (chifixed_factor - 1.D0) * frac)
          AKEXT_I(nr) = AKEXT_I(nr) * (1.D0 + (chifixed_factor - 1.D0) * frac)
       END DO
       WRITE(6,'(A,F8.3,A)') '## tr_prep_chifixed: mode=3, edge chi multiplied by ',chifixed_factor,' for r/a > 0.5'
    END IF

    WRITE(6,'(A,2ES12.4)') '## tr_prep_chifixed: chi_e range = ',MINVAL(AKEXT_E(1:NRMAX)),MAXVAL(AKEXT_E(1:NRMAX))
    WRITE(6,'(A,2ES12.4)') '## tr_prep_chifixed: chi_i range = ',MINVAL(AKEXT_I(1:NRMAX)),MAXVAL(AKEXT_I(1:NRMAX))

    RETURN
  END SUBROUTINE tr_prep_chifixed

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
                PRL_ext(nr) = (1.D0-frac)*prl_csv(j) + frac*prl_csv(j+1)
                PRL_ext(nr) = PRL_ext(nr) * 1.D6
                EXIT
             END IF
          END DO
       END IF
    END DO

    RETURN
  END SUBROUTINE tr_prep_prlfixed

END MODULE trfixed
