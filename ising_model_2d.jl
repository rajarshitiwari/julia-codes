
PI = 4.0*atan(1.0)
Int :: indices(:,:)
Int :: minusplus(:,:)

function total_energy_ising_2d(spins_lattice, jey1, jey2, magnetic_field)::Float64
    # INTEGER, DIMENSION(:,:), INTENT(IN) :: SPINS_LATTICE
    # REAL(8), INTENT(IN) :: JEY1,JEY2
    # REAL(8), INTENT(IN) :: MAGNETIC_FIELD

    co = zeros(Int, 2)
    xi = zeros(Int, 2)
    
    lsize = size(spins_lattice)[1]
    total_energy_ising_2d = 0.0
    energy_nn = 0
    energy_nnn = 0
    for i1 = 1:lsize, i2 = 1:lsize

        xi .= [i1, i2]    #NEAREST NEIGHBOR PART
        for m1 = 1:2, m2 = 1:1
            co .= xi
            co[m1] = minusplus[xi[M1], M2]
            energy_nn = energy_nn + spins_lattice[i1, i2] * spins_lattice[co[1], co[2]]
        end

        #NEXT NEAREST NEIGHBOR PART
        for m1 = 1:2, m2 = 1:1
            co[1] = minusplus[i1, m1]
            co[2] = minusplus[i2, m2]
            energy_nnn = energy_nnn + spins_lattice[i1, i2] * spins_lattice[co[1], co[2]]
        end
    end

    total_energy_ising_2d = jey1 * energy_nn + jey2 * energy_nnn
    total_s = sum(spins_lattice)
    total_energy_ising_2d = total_energy_ising_2d - magnetic_field * total_s
    return
end # function total_energy_ising_2d

SUBROUTINE CALCULATE_AVERAGE_SSSQ(SPINS,TOTAL_S,TOTAL_SSQ)
    USE MYTYPE, ONLY : INDICES
    IMPLICIT NONE
    INTEGER, DIMENSION(:,:), INTENT(IN) :: SPINS
    REAL(8), INTENT(OUT) :: TOTAL_S,TOTAL_SSQ

    TOTAL_S = SUM(SPINS(:,:))
    TOTAL_SSQ = TOTAL_S**2

    RETURN
  END SUBROUTINE CALCULATE_AVERAGE_SSSQ

  SUBROUTINE CREATE_INDICES(LSIZE)
    USE MYTYPE, ONLY :INDICES
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: LSIZE
    INTEGER :: I1,I2

    ALLOCATE(INDICES(LSIZE,LSIZE))

    FORALL (I1 = 1:LSIZE,I2 = 1:LSIZE)
       INDICES(I1,I2) = LSIZE * (I1 - 1) + I2
    END FORALL
    RETURN
  END SUBROUTINE CREATE_INDICES

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ###
  SUBROUTINE GENERATE_RANDOM_SPIN(SPINCONF)
    USE MYTYPE
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: SPINCONF
    REAL(8) :: TEMP
    CALL RANDOM_NUMBER(TEMP)
    IF (TEMP >= 0.5D0)THEN
       SPINCONF = 1
    ELSE
       SPINCONF = -1
    END IF
    RETURN
  END SUBROUTINE GENERATE_RANDOM_SPIN

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ###
### TAKE THE SPIN CONFIGURATION AND CALCULATE SI.SJ      ###
  SUBROUTINE CALCULATE_SI_SJ(SPINS,SI_SJ)
    USE MYTYPE, ONLY : INDICES
    IMPLICIT NONE
    INTEGER, DIMENSION(:,:), INTENT(IN) :: SPINS
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: SI_SJ
    INTEGER :: I1,I2,J1,J2,LSIZE

    LSIZE = SIZE(SPINS,DIM = 1)

    SI_SJ = 0

    #ALL THE SITES HAVE NON ZERO SPINS
    DO I1 = 1, LSIZE
       DO I2 = 1, LSIZE

          DO J1 = 1, LSIZE
             DO J2 = 1, LSIZE
                IF( INDICES(I1,I2) < INDICES(J1,J2)) CYCLE
                SI_SJ(INDICES(I1,I2),INDICES(J1,J2)) = SPINS(I1,I2) * SPINS(J1,J2)
                SI_SJ(INDICES(J1,J2),INDICES(I1,I2)) = SI_SJ(INDICES(I1,I2),INDICES(J1,J2))
             END DO
          END DO

       END DO
    END DO

  END SUBROUTINE CALCULATE_SI_SJ

  SUBROUTINE SELECTED_SQ_FROM_SI_SJ(SI_SJ,QUE,LSIZE,STRUCTURE_FACTOR)

    USE MYTYPE, ONLY: PI,INDICES
    IMPLICIT NONE
    INTEGER, DIMENSION(:,:), INTENT(IN) :: SI_SJ
    INTEGER, DIMENSION(:), INTENT(IN) :: QUE
    INTEGER, INTENT(IN) :: LSIZE
    REAL(8), INTENT(OUT) :: STRUCTURE_FACTOR

    REAL(8) :: TMPDP,COEFF
    INTEGER :: Q(2),RIJ(2),NORMALIZE
    INTEGER :: I1,I2,J1,J2,ITMP

    COEFF = 2.0D0 * PI/DBLE(LSIZE)
    NORMALIZE = LSIZE**4
    TMPDP = 0.0D0
    DO ITMP = 1,LSIZE**2
       TMPDP = TMPDP + DBLE(SI_SJ(ITMP,ITMP))
    END DO
    STRUCTURE_FACTOR = TMPDP

    Q = QUE
    DO I1 = 1,LSIZE
       DO J1 = 1,LSIZE
          RIJ(1) = I1 - J1
          DO I2 = 1,LSIZE
             DO J2 = 1,LSIZE
                RIJ(2) = I2 - J2
                IF(INDICES(J1,J2) >= INDICES(I1,I2))CYCLE
                ITMP = DOT_PRODUCT(Q,RIJ)
                TMPDP = COEFF * DBLE(ITMP)
                TMPDP = DBLE(SI_SJ(INDICES(I1,I2),INDICES(J1,J2))) * 2.0D0 * COS(TMPDP)
                STRUCTURE_FACTOR = STRUCTURE_FACTOR + TMPDP
             END DO
          END DO
       END DO
    END DO

    STRUCTURE_FACTOR = STRUCTURE_FACTOR / DBLE(NORMALIZE)
    RETURN

  END SUBROUTINE SELECTED_SQ_FROM_SI_SJ

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ###
### GIVEN SI.SJ FOR THE SPIN CONFIGURATION, CALCULATE    ###
### THE STRUCTURE FACTOR S(Q) FROM THIS                  ###
  SUBROUTINE STRUCTURE_FACTOR_FROM_SI_SJ(SI_SJ,STRUCTURE_FACTOR)

    USE MYTYPE, ONLY: PI,INDICES
    IMPLICIT NONE
    INTEGER, DIMENSION(:,:), INTENT(IN) :: SI_SJ
    REAL(8), DIMENSION(:,:), INTENT(OUT) :: STRUCTURE_FACTOR

    REAL(8) :: TMPDP,COEFF
    INTEGER :: LSIZE
    INTEGER :: Q1,Q2,Q(2),RIJ(2),NORMALIZE
    INTEGER :: I1,I2,J1,J2,ITMP

    LSIZE = SIZE(STRUCTURE_FACTOR,DIM = 1)
    COEFF = 2.0D0 * PI/DBLE(LSIZE)

    NORMALIZE = LSIZE**4

    TMPDP = 0.0D0
    DO ITMP = 1,LSIZE**2
       TMPDP = TMPDP + DBLE(SI_SJ(ITMP,ITMP))
    END DO
    STRUCTURE_FACTOR = TMPDP

    LOOP_Q1:DO Q1 = 1,LSIZE
       Q(1) = Q1 - 1
       LOOP_Q2:DO Q2 = 1,LSIZE
          Q(2) = Q2 - 1

          DO I1 = 1,LSIZE
             DO J1 = 1,LSIZE
                RIJ(1) = I1 - J1
                DO I2 = 1,LSIZE
                   DO J2 = 1,LSIZE
                      RIJ(2) = I2 - J2
                      IF(INDICES(J1,J2) >= INDICES(I1,I2))CYCLE
                      ITMP = DOT_PRODUCT(Q,RIJ)
                      TMPDP = COEFF * DBLE(ITMP)
                      TMPDP = DBLE(SI_SJ(INDICES(I1,I2),INDICES(J1,J2))) * 2.0D0 * COS(TMPDP)
                      STRUCTURE_FACTOR(Q1,Q2) = STRUCTURE_FACTOR(Q1,Q2) + TMPDP
                   END DO
                END DO
             END DO
          END DO

       END DO LOOP_Q2
    END DO LOOP_Q1

    STRUCTURE_FACTOR = STRUCTURE_FACTOR / DBLE(NORMALIZE)
    RETURN

  END SUBROUTINE STRUCTURE_FACTOR_FROM_SI_SJ

# PROGRAM COOL_ISING_SPINS_3D

#   USE MYTYPE

#   IMPLICIT NONE

#   INTEGER, PARAMETER :: LSIZE = 20
#   INTEGER, PARAMETER :: NSIZE = LSIZE**2

#   REAL(8) :: JEY1,JEY2
#   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: SI_SJ,AVG_SISJ
#   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: SPINS_LATTICE
#   INTEGER :: TMP_SPIN
#   INTEGER :: MCSWEEP     ### INSTANTANEOUS NUMBER OF MC SWEEPS   ###
#   INTEGER :: MAXMCSWEEP  ### MAXIMUM NUMBER OF MC SWEEPS         ###
#   INTEGER :: NEQUIL      ### NUMBER OF MC SWEEPS FOR EQUILIBRIUM ###
#   INTEGER :: I,I1,I2
#   INTEGER :: M1,M2,CO(2),XI(2)
#   INTEGER :: COUNT,IMC,IPASS,IFAIL
#   INTEGER :: OLDCOUNT
#   INTEGER :: OUTCOUNT
#   INTEGER :: NDEL
#   INTEGER :: UNIT_ES,UNIT_MX,UNIT_SQ,UNIT_SP
#   INTEGER :: UNIT_SISJ,UNIT_TMP,UNIT_EPS,UNIT_ACC
#   INTEGER :: NUMBER_OF_TEMPERATURE,COUNT_TEMP
#   REAL(8) :: TEMPERATURE
#   REAL(8) :: MAX_TEMPERATURE
#   REAL(8) :: DELTA_T
#   REAL(8) :: MAGNETIC_FIELD,TOTAL_S,PREVIOUS_S
#   REAL(8) :: CURRENT_S,CHAI,AVERAGE_S,AVERAGE_SSQ,S_SQUARE
#   REAL(8) :: PREVIOUS_ENERGY
#   REAL(8) :: FINAL_ENERGY,MAGNETIZATION,SQAF
#   REAL(8) :: ENERGY_DIFFERENCE
#   REAL(8) :: DIVIDE,MAX_SQ
#   REAL(8) :: AVGE,AVGESQ,C_V,TOTAL_ENERGY
#   REAL(8) :: METRO_POLICE,COMPARE
#   REAL(8) :: TIME1,TIME2,TSF1,TSF2#,X05BAF,

#   #REAL(8) :: TMP(10)
#   CHARACTER(LEN=12) :: TMPSAVE,QSTR1,QSTR2
#   #CHARACTER(LEN=32) :: ARG
#   REAL(8), ALLOCATABLE, DIMENSION(:,:) :: SF_ISING

#   UNIT_ES   = 100
#   UNIT_MX   = 101
#   UNIT_SP   = 102
#   UNIT_SQ   = 103
#   UNIT_SISJ = 104
#   UNIT_EPS  = 105
#   UNIT_ACC  = 106
#   UNIT_TMP  = 111

#   #DEFINE THE TIME OF SWEEP AND NUMBER OF STEPS TO EQUILIBRIATE.
#   MAXMCSWEEP = 10000#0
#   NEQUIL = 5000#0
#   NDEL = 20
#   NUMBER_OF_TEMPERATURE = 60

#   #======READ THE INPUT VARIABLES FIRSTLY FROM ARGUMENT
#   #IF THAT NOT GIVEN THEN FROM STANDARD INPUT=======#
#   JEY1 = -1.0D0
#   PRINT*,"ENTER JEY2,MAX TEMPERATURE"
#   READ(*,*)JEY2,MAX_TEMPERATURE #JEY2 = 0.5D0

#   MAGNETIC_FIELD = 0.0D0
#   #MAX_TEMPERATURE = 3.0D0 * ABS(JEY1)
#   DELTA_T = MAX_TEMPERATURE / DBLE(NUMBER_OF_TEMPERATURE - 1)

#   ALLOCATE(SI_SJ(NSIZE,NSIZE),AVG_SISJ(NSIZE,NSIZE))
#   ALLOCATE(SF_ISING(0:LSIZE-1,0:LSIZE-1))
#   ALLOCATE(SPINS_LATTICE(LSIZE,LSIZE))
#   ALLOCATE(MINUSPLUS(LSIZE,2),ONETON(LSIZE))

#   CALL CREATE_INDICES(LSIZE)

#   FORALL(I = 1:LSIZE)ONETON(I) = I
#   MINUSPLUS(:,1) = CSHIFT(ONETON,-1)
#   MINUSPLUS(:,2) = CSHIFT(ONETON,1)

#   OPEN(UNIT_TMP, FILE = 'ising_2d.inp')
#   WRITE(UNIT_TMP,'("JEY1            = ",F15.5)')JEY1
#   WRITE(UNIT_TMP,'("JEY2            = ",F15.5)')JEY2
#   WRITE(UNIT_TMP,'("LSIZE           = ",I15)')LSIZE
#   WRITE(UNIT_TMP,'("NSIZE           = ",I15)')NSIZE
#   WRITE(UNIT_TMP,'("SWEEP LENGTH    = ",I15)')MAXMCSWEEP
#   WRITE(UNIT_TMP,'("EQUIL LENGTH    = ",I15)')NEQUIL
#   WRITE(UNIT_TMP,'("NO. OF T        = ",I15)')NUMBER_OF_TEMPERATURE
#   WRITE(UNIT_TMP,'("MAXIMUM T       = ",F15.5)')MAX_TEMPERATURE
#   WRITE(UNIT_TMP,'("DELTA T         = ",F15.5)')DELTA_T
#   WRITE(UNIT_TMP,'("MAGNETIC FIELD  = ",3(F15.5,1X))')MAGNETIC_FIELD

#   CLOSE(UNIT_TMP)

#   WRITE(*,'("JEY1            = ",F15.5)')JEY1
#   WRITE(*,'("JEY2            = ",F15.5)')JEY2
#   WRITE(*,'("LSIZE           = ",I15)')LSIZE
#   WRITE(*,'("NSIZE           = ",I15)')NSIZE
#   WRITE(*,'("SWEEP LENGTH    = ",I15)')MAXMCSWEEP
#   WRITE(*,'("EQUIL LENGTH    = ",I15)')NEQUIL
#   WRITE(*,'("NO. OF T        = ",I15)')NUMBER_OF_TEMPERATURE
#   WRITE(*,'("MAXIMUM T       = ",F15.5)')MAX_TEMPERATURE
#   WRITE(*,'("DELTA T         = ",F15.5)')DELTA_T
#   WRITE(*,'("MAGNETIC FIELD  = ",3(F15.5,1X))')MAGNETIC_FIELD

#   #GENERATE THE SNAPSHOT OF SPINS (RANDOM CONFIGURATION)#
#   OPEN(UNIT_TMP,FILE = 'starting_spins.dat')
#   DO I1 = 1,LSIZE
#      DO I2 = 1,LSIZE
#         #SPINS_LATTICE(I1,I2) = 1
#         CALL GENERATE_RANDOM_SPIN(SPINS_LATTICE(I1,I2))
#         WRITE(UNIT_TMP,"(3(I2,1X))")I1,I2,SPINS_LATTICE(I1,I2)
#      END DO
#   END DO
#   CLOSE(UNIT_TMP)

#   #=====================================#
#   #=======OPEN THE FILES TO WRITE=======#
#   OPEN(UNIT_ES, FILE = 'energy_susceptibility.dat')
#   OPEN(UNIT_MX, FILE = 'magnetization_chai.dat')
#   OPEN(UNIT_SP, FILE = 'spin_configurations.dat')
#   OPEN(UNIT_SQ, FILE = 'structure_factors.dat')
#   OPEN(UNIT_SISJ, FILE = 'avgsisj.dat')
#   OPEN(UNIT_ACC, FILE = 'acceptence.dat')

#   WRITE(UNIT_ES,'(3(A15,1X))')'TEMPERATURE','AVERAGE_ENERGY','SPECIFIC_HEAT'
#   WRITE(UNIT_MX,'(6(A15,1X))')'TEMPERATURE','SQ_FERRO','SQ_G_TYPE_AF',&
#        &'CHAI','S','SSQ'
#   #======================================#

#   TOTAL_ENERGY = TOTAL_ENERGY_ISING_2D(SPINS_LATTICE,JEY1,JEY2,MAGNETIC_FIELD)
#   PRINT'("TOTAL ENERGY IN STARTING = ",F20.10)',TOTAL_ENERGY/NSIZE
#   PRINT"('STARTING THE T LOOP')"

#   TOTAL_S = SUM(SPINS_LATTICE(:,:))

#   PREVIOUS_ENERGY = TOTAL_ENERGY
#   TEMPERATURE = MAX_TEMPERATURE
#   COUNT_TEMP = 0

#   TEMPERATURE_LOOP : DO WHILE(TEMPERATURE > 0.0D0)

#      MCSWEEP = 1
#      COUNT = 0
#      OLDCOUNT = 0
#      OUTCOUNT = 0
#      AVG_SISJ = 0.0D0
#      AVGE = 0.0D0
#      AVGESQ = 0.0D0
#      MAGNETIZATION = 0.0D0
#      SQAF = 0.0D0
#      AVERAGE_S = 0.0D0
#      AVERAGE_SSQ = 0.0D0
#      CHAI = 0.0D0
#      IMC = 0
#      IFAIL = 0
#      IPASS = 0
#      CALL CPU_TIME(TIME1) #TIME1 = X05BAF()
#      MC_LOOP : DO WHILE(MCSWEEP <= MAXMCSWEEP)

#         SWEEP_I1 : DO I1 = 1, LSIZE
#            SWEEP_I2 : DO I2 = 1, LSIZE

#               TMP_SPIN = - SPINS_LATTICE(I1,I2) #FLIP THE TEMP SPIN

#               #@######################################
#               #@ CALCULATE THE ENERGY DIFFERENCE NOW #
#               #@######################################

#               ENERGY_DIFFERENCE = 0.0D0
#               XI = [I1,I2] #NEAREST NEIGHBOR PART
#               DO M1 = 1,2
#                  DO M2 = 1,2
#                     CO = XI
#                     CO(M1) = MINUSPLUS(XI(M1),M2)
#                     ENERGY_DIFFERENCE = ENERGY_DIFFERENCE + 2.0D0 * JEY1 * DBLE(TMP_SPIN * SPINS_LATTICE(CO(1),CO(2)))
#                  END DO
#               END DO

#               XI = [I1,I2]
#               DO M1 = 1,2
#                  DO M2 = 1,2
#                     CO(1) = MINUSPLUS(I1,M1)
#                     CO(2) = MINUSPLUS(I2,M2)
#                     ENERGY_DIFFERENCE = ENERGY_DIFFERENCE + 2.0D0 * JEY2 * DBLE(TMP_SPIN * SPINS_LATTICE(CO(1),CO(2)))
#                  END DO
#               END DO

#               PREVIOUS_S = SPINS_LATTICE(I1,I2)
#               CURRENT_S = TMP_SPIN

#               ENERGY_DIFFERENCE = ENERGY_DIFFERENCE - MAGNETIC_FIELD*CURRENT_S + MAGNETIC_FIELD*PREVIOUS_S

#               IF(ENERGY_DIFFERENCE < 0.0D0)THEN
#                  SPINS_LATTICE(I1,I2) = TMP_SPIN #UPDATE THE TEMP SPIN
#                  FINAL_ENERGY = PREVIOUS_ENERGY + ENERGY_DIFFERENCE
#                  IPASS = IPASS + 1
#               ELSE
#                  #CALL RAN3(METRO_POLICE)
#                  CALL RANDOM_NUMBER(METRO_POLICE)
#                  COMPARE = ENERGY_DIFFERENCE
#                  COMPARE = COMPARE/(TEMPERATURE + 1.0E-10_8)
#                  COMPARE = EXP(-COMPARE)
#                  IF(METRO_POLICE < COMPARE)THEN
#                     SPINS_LATTICE(I1,I2) = TMP_SPIN #UPDATE THE TEMP SPIN
#                     FINAL_ENERGY = PREVIOUS_ENERGY + ENERGY_DIFFERENCE
#                     IPASS = IPASS + 1
#                  ELSE
#                     FINAL_ENERGY = PREVIOUS_ENERGY
#                     IFAIL = IFAIL + 1
#                  END IF
#               END IF
#               IMC = IMC + 1
#            END DO SWEEP_I2
#         END DO SWEEP_I1
        
#         IF(MCSWEEP > NEQUIL)THEN
#            COUNT = COUNT + 1
#            TAKING_AVG :IF(COUNT - OLDCOUNT == NDEL)THEN
#               OUTCOUNT = OUTCOUNT + 1
              
#               CALL CALCULATE_SI_SJ(SPINS_LATTICE,SI_SJ)
#               CALL CALCULATE_AVERAGE_SSSQ(SPINS_LATTICE,TOTAL_S,S_SQUARE)
#               #TOTAL_S = SUM(SPINS_LATTICE(:,:))
#               #S_SQUARE = TOTAL_S**2
              
#               TOTAL_ENERGY = TOTAL_ENERGY_ISING_2D(SPINS_LATTICE,JEY1,JEY2,MAGNETIC_FIELD)

#               AVERAGE_S = AVERAGE_S + TOTAL_S
#               AVERAGE_SSQ = AVERAGE_SSQ + S_SQUARE
#               AVG_SISJ = AVG_SISJ + SI_SJ
#               AVGESQ = AVGESQ + TOTAL_ENERGY**2
#               AVGE = AVGE + TOTAL_ENERGY

#               OLDCOUNT = COUNT
#            END IF TAKING_AVG
#         END IF

#         MCSWEEP = MCSWEEP + 1
#      END DO MC_LOOP

#      DIVIDE = DBLE(OUTCOUNT)
#      AVG_SISJ = AVG_SISJ/DIVIDE
#      AVGE = AVGE/DIVIDE
#      AVGESQ = AVGESQ/DIVIDE
#      AVERAGE_S = AVERAGE_S/DIVIDE
#      AVERAGE_SSQ = AVERAGE_SSQ/DIVIDE

#      WRITE(UNIT_ACC,'(3(F15.10,1X))')TEMPERATURE,DBLE(IPASS)/DBLE(IMC),DBLE(IFAIL)/DBLE(IMC)
     
# ###########################################################
# ###########################################################
#      CALL CPU_TIME(TIME2) #TIME2 = X05BAF()#
#      WRITE(*,"(A,2X,F15.10)",ADVANCE='NO')'RUNNING PROGRAM AT T =',TEMPERATURE
#      WRITE(*,'(3X,"TIME TAKEN =",F12.2," SECONDS")',ADVANCE='NO')(TIME2-TIME1)

#      CALL CPU_TIME(TSF1)        #TSF1 = X05BAF()
#      CALL STRUCTURE_FACTOR_FROM_SI_SJ(AVG_SISJ,SF_ISING)
#      CALL CPU_TIME(TSF2)        #TSF2 = X05BAF()
#      MAX_SQ = MAXVAL(SF_ISING)
#      WRITE(*,'(3X,"TIME FOR S(Q) =",1X,F12.2,1X,"SECONDS. MAXIMUM S(Q) =",1X,F20.10)')TSF2-TSF1,MAX_SQ

#      C_V = (AVGESQ - (AVGE)**2)/DBLE(NSIZE)
#      C_V = C_V / (TEMPERATURE**2)
#      AVGE = AVGE/NSIZE

#      CHAI = AVERAGE_SSQ/DBLE(NSIZE) - (AVERAGE_S)**2 /DBLE(NSIZE)
#      CHAI = CHAI / TEMPERATURE
#      SQAF = SF_ISING(LSIZE/2,LSIZE/2)
#      MAGNETIZATION = SF_ISING(0,0)

#      WRITE(UNIT_ES,'(2(F20.10,1X),E20.10)')TEMPERATURE,AVGE,C_V
#      WRITE(UNIT_MX,'(6(F20.10,1X))')TEMPERATURE,MAGNETIZATION,SQAF,CHAI,AVERAGE_S/DBLE(NSIZE),AVERAGE_SSQ/DBLE(NSIZE)

#      WRITE(UNIT_SP,*)'T=',TEMPERATURE
#      DO I1 = 1, LSIZE
#         DO I2 = 1, LSIZE
#            WRITE(UNIT_SP,'(3(I2,1X))')I1,I2,SPINS_LATTICE(I1,I2)
#            CALL NUMTOSTR(I1-1,QSTR1)
#            CALL NUMTOSTR(I2-1,QSTR2)
#            TMPSAVE = TRIM(QSTR1)//'_'//TRIM(QSTR2)
#            WRITE(UNIT_SISJ,'(3(I2,1X))')I1,I2,AVG_SISJ(1,INDICES(I1,I2))
#            WRITE(UNIT_SQ,'(A9,1X,1(F20.10,1X),E20.10)')TMPSAVE,TEMPERATURE,SF_ISING(I1-1,I2-1)
#         END DO
#      END DO
#      WRITE(UNIT_SP,*)
#      WRITE(UNIT_SQ,*)

#      TEMPERATURE = TEMPERATURE - DELTA_T
#      COUNT_TEMP = COUNT_TEMP + 1

#   END DO TEMPERATURE_LOOP

#   CLOSE(UNIT_ES)
#   CLOSE(UNIT_MX)
#   CLOSE(UNIT_SP)
#   CLOSE(UNIT_SQ)
#   CLOSE(UNIT_SISJ)
#   CLOSE(UNIT_EPS)
#   CLOSE(UNIT_TMP)

#   OPEN(101,FILE = 'spin_gs.dat')
#   OPEN(102,FILE = 'sisj_gs.dat')
#   DO I1 = 1, LSIZE
#      DO I2 = 1, LSIZE
#         WRITE(101,'(3(I2,1X))')I1,I2,SPINS_LATTICE(I1,I2)
#         WRITE(102,'(3(I2,1X))')I1,I2,SI_SJ(1,INDICES(I1,I2))
#      END DO;WRITE(101,*);WRITE(102,*)
#   END DO
#   CLOSE(101)
#   CLOSE(102)
#   PRINT*,'GROUND STATE SPIN SAVED IN FILE=','spin_gs.dat'
#   PRINT*,'LOWEST E',TOTAL_ENERGY_ISING_2D(SPINS_LATTICE,JEY1,JEY2,MAGNETIC_FIELD)/DBLE(NSIZE)
#   PRINT*,'NUMBER OF Ts =',COUNT_TEMP

# CONTAINS
# END PROGRAM COOL_ISING_SPINS_3D
