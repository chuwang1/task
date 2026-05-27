PROGRAM validate_onetwo_prc_driver

  USE bpsd_kinds, ONLY: rkind
  USE tr_onetwo_prc_kernel, ONLY: onetwo_prc_point
  IMPLICIT NONE

  INTEGER, PARAMETER :: input_unit = 20
  INTEGER, PARAMETER :: output_unit = 21

  CHARACTER(LEN=512) :: input_path
  CHARACTER(LEN=512) :: output_path
  CHARACTER(LEN=512) :: arg
  CHARACTER(LEN=1024) :: line
  INTEGER :: ios
  INTEGER :: narg
  INTEGER :: nrow
  REAL(KIND=rkind) :: bt_t
  REAL(KIND=rkind) :: a_m
  REAL(KIND=rkind) :: r0_m
  REAL(KIND=rkind) :: wall_reflection
  REAL(KIND=rkind) :: rho
  REAL(KIND=rkind) :: ne_1e20_m3
  REAL(KIND=rkind) :: te_kev
  REAL(KIND=rkind) :: prc_w_m3
  REAL(KIND=rkind) :: phi_bar

  narg = COMMAND_ARGUMENT_COUNT()
  IF (narg /= 6) THEN
     WRITE(*,'(A)') 'Usage: onetwo_prc_driver input.csv output.csv bt_T a_m R0_m wall_reflection'
     STOP 2
  END IF

  CALL GET_COMMAND_ARGUMENT(1, input_path)
  CALL GET_COMMAND_ARGUMENT(2, output_path)

  CALL GET_COMMAND_ARGUMENT(3, arg)
  READ(arg,*,IOSTAT=ios) bt_t
  IF (ios /= 0) STOP 'Failed to parse bt_T'

  CALL GET_COMMAND_ARGUMENT(4, arg)
  READ(arg,*,IOSTAT=ios) a_m
  IF (ios /= 0) STOP 'Failed to parse a_m'

  CALL GET_COMMAND_ARGUMENT(5, arg)
  READ(arg,*,IOSTAT=ios) r0_m
  IF (ios /= 0) STOP 'Failed to parse R0_m'

  CALL GET_COMMAND_ARGUMENT(6, arg)
  READ(arg,*,IOSTAT=ios) wall_reflection
  IF (ios /= 0) STOP 'Failed to parse wall_reflection'

  OPEN(UNIT=input_unit, FILE=TRIM(input_path), STATUS='OLD', ACTION='READ', &
       IOSTAT=ios)
  IF (ios /= 0) STOP 'Failed to open input CSV'

  OPEN(UNIT=output_unit, FILE=TRIM(output_path), STATUS='REPLACE', &
       ACTION='WRITE', IOSTAT=ios)
  IF (ios /= 0) STOP 'Failed to open output CSV'

  READ(input_unit,'(A)',IOSTAT=ios) line
  IF (ios /= 0) STOP 'Failed to read input CSV header'

  WRITE(output_unit,'(A)') 'rho,ne_1e20_m3,te_kev,prc_w_m3,prc_mw_m3,phi_bar'

  nrow = 0
  DO
     READ(input_unit,'(A)',IOSTAT=ios) line
     IF (ios < 0) EXIT
     IF (ios > 0) STOP 'Failed while reading input CSV'
     IF (LEN_TRIM(line) == 0) CYCLE

     CALL replace_commas_with_spaces(line)
     READ(line,*,IOSTAT=ios) rho, ne_1e20_m3, te_kev
     IF (ios /= 0) STOP 'Failed to parse input CSV row'

     CALL onetwo_prc_point(ne_1e20_m3, te_kev, bt_t, a_m, r0_m, &
          wall_reflection, prc_w_m3, phi_bar)

     WRITE(output_unit,'(*(ES24.16,:,","))') rho, ne_1e20_m3, te_kev, &
          prc_w_m3, prc_w_m3*1.0E-6_rkind, phi_bar
     nrow = nrow + 1
  END DO

  CLOSE(input_unit)
  CLOSE(output_unit)

  WRITE(*,'(A,I0,A)') 'Wrote ', nrow, ' ONETWO PRC rows.'

CONTAINS

  SUBROUTINE replace_commas_with_spaces(text)
    CHARACTER(LEN=*), INTENT(INOUT) :: text
    INTEGER :: i

    DO i = 1, LEN(text)
       IF (text(i:i) == ',') text(i:i) = ' '
    END DO
  END SUBROUTINE replace_commas_with_spaces

END PROGRAM validate_onetwo_prc_driver
