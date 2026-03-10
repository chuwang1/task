      SUBROUTINE TR_WRITE_CSV(GX,GY,NXM,NXMAX,NGMAX,STR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NXM, NXMAX, NGMAX
      REAL, INTENT(IN) :: GX(NXMAX), GY(NXM,NGMAX)
      CHARACTER(LEN=*), INTENT(IN) :: STR

      INTEGER, SAVE :: N_CSV_COUNT = 0
      CHARACTER(LEN=256) :: FNAME, TITLE_STR, VAR_PART, VAR_NAME, NEW_NAME
      CHARACTER(LEN=64) :: HEADERS(100)
      INTEGER :: I, J, UNIT_CSV, K, POS_LBRACKET, POS_VS, POS_COMMA, NUM_VARS
      INTEGER :: START_IDX, END_IDX, POS_NS
      CHARACTER(LEN=2) :: SUFFIX_ARR(10)

      DATA SUFFIX_ARR /'E ', 'D ', 'T ', 'A ', 'S5', 'S6', 'S7', 'S8', 'S9', 'SX'/

      N_CSV_COUNT = N_CSV_COUNT + 1
      UNIT_CSV = 80 + MOD(N_CSV_COUNT, 20)

      WRITE(FNAME, '(A,I3.3,A)') 'tr_data_', N_CSV_COUNT, '.csv'
      OPEN(UNIT=UNIT_CSV, FILE=TRIM(FNAME), STATUS='REPLACE', ACTION='WRITE')

      WRITE(UNIT_CSV, *) 'Title: ', TRIM(STR)

      TITLE_STR = STR
      IF (LEN_TRIM(TITLE_STR) > 0 .AND. TITLE_STR(1:1) == '@') THEN
          TITLE_STR = TITLE_STR(2:)
      END IF
      K = LEN_TRIM(TITLE_STR)
      IF (K > 0 .AND. TITLE_STR(K:K) == '@') THEN
          TITLE_STR(K:K) = ' '
      END IF

      POS_LBRACKET = INDEX(TITLE_STR, '[')
      POS_VS = INDEX(TITLE_STR, ' vs ')

      END_IDX = LEN_TRIM(TITLE_STR)
      IF (POS_LBRACKET > 0) END_IDX = MIN(END_IDX, POS_LBRACKET - 1)
      IF (POS_VS > 0) END_IDX = MIN(END_IDX, POS_VS - 1)

      VAR_PART = ' '
      IF (END_IDX > 0) VAR_PART = TITLE_STR(1:END_IDX)

      HEADERS = ''
      START_IDX = 1
      NUM_VARS = 0

      DO J = 1, MIN(NGMAX, 100)
          POS_COMMA = INDEX(VAR_PART(START_IDX:), ',')
          IF (POS_COMMA > 0) THEN
              HEADERS(J) = ADJUSTL(VAR_PART(START_IDX:START_IDX+POS_COMMA-2))
              START_IDX = START_IDX + POS_COMMA
              NUM_VARS = NUM_VARS + 1
          ELSE
              IF (LEN_TRIM(VAR_PART(START_IDX:)) > 0) THEN
                  HEADERS(J) = ADJUSTL(VAR_PART(START_IDX:))
                  NUM_VARS = NUM_VARS + 1
              END IF
              EXIT
          END IF
      END DO

      WRITE(UNIT_CSV, '(A)', ADVANCE='NO') 'X'

      POS_NS = 0
      IF (NUM_VARS == 1) THEN
          VAR_NAME = HEADERS(1)
          POS_NS = INDEX(VAR_NAME, '(NS)')
          IF (POS_NS == 0) POS_NS = INDEX(VAR_NAME, '(ns)')
      END IF

      DO J = 1, NGMAX
          WRITE(UNIT_CSV, '(A)', ADVANCE='NO') ','

          IF (NUM_VARS == 1 .AND. NGMAX > 1) THEN
              VAR_NAME = HEADERS(1)
              IF (POS_NS > 0) THEN
                  VAR_NAME = VAR_NAME(1:POS_NS-1)
                  IF (J <= 10) THEN
                      NEW_NAME = TRIM(VAR_NAME) // TRIM(SUFFIX_ARR(J))
                  ELSE
                      WRITE(NEW_NAME, '(A,A,I0)') TRIM(VAR_NAME), '_S', J
                  END IF
              ELSE
                  WRITE(NEW_NAME, '(A,A,I0)') TRIM(VAR_NAME), '_', J
              END IF
              WRITE(UNIT_CSV, '(A)', ADVANCE='NO') TRIM(NEW_NAME)

          ELSEIF (J <= NUM_VARS .AND. J <= 100) THEN
              WRITE(UNIT_CSV, '(A)', ADVANCE='NO') TRIM(HEADERS(J))
          ELSE
              WRITE(UNIT_CSV, '(A,I0)', ADVANCE='NO') 'Y', J
          END IF
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

      SUBROUTINE TR_WRITE_CSV_2D(GT,GX,GF,NTM,NTMAX,NXM,NXMAX,STR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NTM, NTMAX, NXM, NXMAX
      REAL, INTENT(IN) :: GT(NTM), GX(NXM), GF(NTM,NXMAX)
      CHARACTER(LEN=*), INTENT(IN) :: STR

      INTEGER, SAVE :: N_CSV_COUNT = 0
      CHARACTER(LEN=256) :: FNAME
      INTEGER :: I, J, UNIT_CSV

      N_CSV_COUNT = N_CSV_COUNT + 1
      UNIT_CSV = 80 + MOD(N_CSV_COUNT, 20)

      WRITE(FNAME, '(A,I3.3,A)') 'tr_data_', N_CSV_COUNT, '.csv'
      OPEN(UNIT=UNIT_CSV, FILE=TRIM(FNAME), STATUS='REPLACE', ACTION='WRITE')

      WRITE(UNIT_CSV, *) 'Title: ', TRIM(STR)
      WRITE(UNIT_CSV, *) 'Type: 2D_CONTOUR'
      WRITE(UNIT_CSV, *) 'NTMAX: ', NTMAX, ' NXMAX: ', NXMAX

      WRITE(UNIT_CSV, '(A)', ADVANCE='NO') 'Time'
      DO J = 1, NXMAX
          WRITE(UNIT_CSV, '(A,E14.7)', ADVANCE='NO') ',', GX(J)
      END DO
      WRITE(UNIT_CSV, *)

      DO I = 1, NTMAX
          WRITE(UNIT_CSV, '(E14.7)', ADVANCE='NO') GT(I)
          DO J = 1, NXMAX
              WRITE(UNIT_CSV, '(A,E14.7)', ADVANCE='NO') ',', GF(I,J)
          END DO
          WRITE(UNIT_CSV, *)
      END DO

      CLOSE(UNIT_CSV)
      WRITE(6, *) ' ## CSV EXPORTED (2D): ', TRIM(FNAME)

      RETURN
      END SUBROUTINE TR_WRITE_CSV_2D
