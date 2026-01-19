      SUBROUTINE TR_WRITE_CSV(GX,GY,NXM,NXMAX,NGMAX,STR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NXM, NXMAX, NGMAX
      REAL, INTENT(IN) :: GX(NXMAX), GY(NXM,NGMAX)
      CHARACTER(LEN=*), INTENT(IN) :: STR

      INTEGER, SAVE :: N_CSV_COUNT = 0
      CHARACTER(LEN=256) :: FNAME
      INTEGER :: I, J, UNIT_CSV

      N_CSV_COUNT = N_CSV_COUNT + 1
      UNIT_CSV = 80 + MOD(N_CSV_COUNT, 20) ! Cycle unit numbers if needed, avoiding standard units

      WRITE(FNAME, '(A,I3.3,A)') 'tr_data_', N_CSV_COUNT, '.csv'
      OPEN(UNIT=UNIT_CSV, FILE=TRIM(FNAME), STATUS='REPLACE', ACTION='WRITE')

      WRITE(UNIT_CSV, *) 'Title: ', TRIM(STR)
      WRITE(UNIT_CSV, '(A)', ADVANCE='NO') 'X'
      DO J = 1, NGMAX
          WRITE(UNIT_CSV, '(A,I0)', ADVANCE='NO') ',Y', J
      END DO
      WRITE(UNIT_CSV, *)

      DO I = 1, NXMAX
          WRITE(UNIT_CSV, '(E14.7)', ADVANCE='NO') GX(I)
          DO J = 1, NGMAX
              WRITE(UNIT_CSV, '(A,E14.7)', ADVANCE='NO') ',', GY(I,J)
          END DO
          WRITE(UNIT_CSV, *)
      END DO

      CLOSE(UNIT_CSV)
      WRITE(6, *) ' ## CSV EXPORTED: ', TRIM(FNAME)

      RETURN
      END SUBROUTINE TR_WRITE_CSV
