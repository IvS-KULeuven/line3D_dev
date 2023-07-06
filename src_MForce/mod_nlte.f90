MODULE mod_nlte
    USE prog_type 
    IMPLICIT NONE
    ! Module for comutation of an approximat NLTE occupation number densities
    ! usin the atomic data files from FASTWIND.
    !
    ! The Module was developed by Luka Poniatowski in 2021
    ! It contains routines ratially based on original LTE_Line_module.F90 
    ! and NLTE corrections based on Lucy & Abbott 1993; Puls et al. (2000).
    !
    ! This module contains the definitions of variables and variable types
    ! used to allocate the atomic, partition functions, occupation numbers, etc.
    !
    ! version 1.0
    ! ----
    !
    ! Identifiers used in the description for definitions of dimensions and variable types:
    ! UI - Unsingen INTEGER
    ! DP - DOUBLE PRECISION
    ! 1D/2D/3D - 1/2/3 dimentional array
    ! AT - ATOMIC_DATA_TYPE
    ! OT - NLTE_NUMBER_TYPE
    ! ----
    !
    ! Within this module variable names are reserd as follows:
    ! Z - Nuclear charge of the atom :: UI
    ! I - Ionisation state   :: UI
    ! L - Energy level  :: UI
    ! Ll / Lu - Lower / Upper Energy level :: UI
    ! Xg - Energy of the gorund level  :: DP
    ! Xl - Excitation energy of a level  :: DP
    ! XI - Ionisaation energy from the ground level of a ion :: DP
    ! gg - Degeneracy of the groung level  :: DP
    ! gL - Degenerecy of a level :: DP
    ! ATOMIC_DATA - Structure containing atomic data :: AT
    ! Number_denc - Structure containing number dencities
    ! ---

    TYPE ATOMIC_DATA_TYPE
            ! This variable type is used to allocate the atomic data

            ! The variables are as follows:

            ! Maximum Nuclear charge Z, Ionisation state I,
            !  and Energy level L available in the atomic table
            INTEGER(I4B) :: MaxZ, MaxI, MaxL

            ! List of maximum ionisation states of each Z
            ! The list has to be 1D array of dimension (MaxZ)
            INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: List_I

            ! List of maximum energy level for ion I of each Z
            ! The list has to be 2D array of dimensions (MaxI, MaxZ)
            INTEGER(I4B), DIMENSION(:, :), ALLOCATABLE :: List_L

            ! Lists of degeneracies, excitation energies and flag for Meta-stable/other
            !  for a sfecific ionisation saste I of an atom Z
            ! Lists must be 3D array of dimensions (MaxL,MaxI,MaxZ)
            REAL(DP), DIMENSION(:, :, :), ALLOCATABLE :: Degeneracy
            REAL(DP), DIMENSION(:, :, :), ALLOCATABLE :: Lev_energy
            INTEGER(I2B) , DIMENSION(:, :, :), ALLOCATABLE :: Lev_is_mstbl

            ! List of ionisation eneries of a specific ion I of atom Z
            ! List must be 2D array of dimensions (MaxI,MaxZ)
            REAL(DP), DIMENSION(:, :), ALLOCATABLE :: Ion_energy

            ! Fraction of recombination processes leading to the ground ion I of atom Z
            ! List must be 2D array of dimensions (MaxI,MaxZ)
            REAL(DP), DIMENSION(:, :), ALLOCATABLE :: Zeta


            ! Atomic weights from fastwind
            REAL(DP) :: Atomic_weight(30) = &
                        [1.008d0, 4.003d0, 6.941d0, 9.012d0, 10.811d0, 12.011d0, &
                            14.007d0, 16.000d0, 18.998d0, 20.180d0, 22.990d0, 24.305d0, &
                            26.982d0, 28.085d0, 30.974d0, 32.066d0, 35.453d0, 39.948d0, &
                            39.098d0, 40.078d0, 44.956d0, 47.88d0, 50.941d0, 51.996d0, &
                            54.938d0, 55.847d0, 58.933d0, 58.69d0, 63.546d0, 65.39d0]

            CHARACTER(len=3) :: Names(30)

        CONTAINS
            PROCEDURE, PASS :: Initialise => Read_Atom_Data
            PROCEDURE, PASS :: Info => Print_ATOMIC_DATA_TYPE_info
    END TYPE ATOMIC_DATA_TYPE
    !========================================================================

    TYPE NLTE_NUMBER_TYPE
            ! This variable type is used to allocate the number densities

            ! The variables are as follows:

            ! Number aboundences of elements with nuclear charge Z realtive to He
            REAL(DP), DIMENSION(:), ALLOCATABLE :: Aboundance

            ! Number dencity of nucleis with nuclear charge Z = rho*Aboundance*Y/(He_mass)
            REAL(DP), DIMENSION(:), ALLOCATABLE :: Nuclei

            ! Partition function of ionisation stage I of element Z (I,Z)
            REAL(DP), DIMENSION(:, :), ALLOCATABLE :: Partition

            ! Number dencity of ions I of element Z (I,Z)
            REAL(DP), DIMENSION(:, :), ALLOCATABLE :: Ions

            ! Occupation number dencities of level L of ion I of element Z (L,I,Z)
            REAL(DP), DIMENSION(:, :, :), ALLOCATABLE :: Occupation

            ! Number dencity of electrones
            REAL(DP) :: Electron

            ! Mass fractions of H and He
            REAL(DP) :: X_frac = 0.0d0
            REAL(DP) :: Y_frac = 0.0d0

            ! Radiation+Equilibrium Temperature assigned to the type
            REAL(DP) :: T = 0.0d0
            ! Electorn Temperature assigned to the type
            REAL(DP) :: Te = 0.0d0
            ! Dilution factor assigned to the type
            REAL(DP) :: W = 0.0d0
            ! Mass dencity assigned to the type
            REAL(DP) :: rho = 0.0d0

            ! h_plank*cc/kb as addopted from fastwind
            REAL(DP) :: Bolz_norm = 1.43883549673335d0

            ! Saha normalization factor (2pi m_e k_B/h^2)^(3/2) addopted from FW
            REAL(DP) :: Saha_norm = 2.414535115067086e+15

            ! Mass of hydrogen
            REAL(DP) :: mH = 1.67372360d-24

            ! Pointer to the Atomic data the occupation numbers are initiated with
            TYPE(ATOMIC_DATA_TYPE), POINTER :: ATOMIC_DATA

        CONTAINS
            PROCEDURE, PASS :: Initialise => Read_Composition
            PROCEDURE, PASS :: Set => Comput_NLTE_Occupation_Numbers
            PROCEDURE, PASS :: Clear => Clear_NLTE_NUMBER_TYPE
            PROCEDURE, PASS :: Info => Print_NLTE_NUMBER_TYPE_info
    END TYPE NLTE_NUMBER_TYPE
    !========================================================================

    PUBLIC :: LoGo 
    !========================================================================

    PRIVATE :: Read_Atom_Data, Print_ATOMIC_DATA_TYPE_info, &
               Read_Composition, Comput_NLTE_Occupation_Numbers, &
               Comput_Partition_Functions, Comput_Ionisation_structure, &
               Comput_Excitation_structure, Print_NLTE_NUMBER_TYPE_info

    !========================================================================
    !========================================================================
    !========================================================================
    !========================================================================

    CONTAINS

        SUBROUTINE Read_Atom_Data(OBJ, verbose)
            ! In-Out and optional variables
            CLASS(ATOMIC_DATA_TYPE), INTENT(INOUT) ::  OBJ
            LOGICAL, OPTIONAL :: verbose

            ! Tepmroray variables used in reading and sortind tables
            LOGICAL :: ver
            INTEGER, PARAMETER  :: Lines_genr_info = 149
            INTEGER :: Core_cahrge(Lines_genr_info) = 0.d0 ! list of nuclear charges Z available
            INTEGER :: Ionis_stage(Lines_genr_info) = 0.d0 ! list of ionisation stages I for each Z
            INTEGER :: Numbr_level(Lines_genr_info) = 0.d0 ! list of number of energy levels E for each I in each Z
            INTEGER :: MaxL = 0, MaxI = 0, MaxZ = 0, Z = 0, I = 0, L = 0, N_L = 0, ind = 0, IO_status = 0, iduma(50) = 0
            REAL(DP) :: Ionisation_cs = 0.d0
            CHARACTER(LEN=10) :: Lev_flag = '         '

            ! Dummy variables for variout uses
            DOUBLE PRECISION :: rdum = 0.d0
            INTEGER :: idum = 0
            CHARACTER(LEN=10) :: cdum = '          '

            ! String containing the dectination to the project directiory
            CHARACTER(256) DATA_DIR

            ! Pars the input optianal variables
            ! If verbose = .true. specified give datailed status
            ver = .FALSE.
            IF (PRESENT(verbose)) ver = verbose

            IF (ver) WRITE(stdout,*) ' '
            IF (ver) WRITE(stdout,*) 'Initialising the Atomic Data Type:'

            ! Get the destination ot the data location (check if the enviroment was set)
            CALL get_environment_variable("MFORCE_DIR", DATA_DIR)
            IF (LEN_TRIM(DATA_DIR) .EQ. 0) STOP '---- Error: MFORCE_DIR not set -----'
            DATA_DIR = TRIM(DATA_DIR)//'/DATA'

            ! Chech if atomic data is already Initialised
            IF (ALLOCATED(OBJ%List_I)) THEN
                IF (ver) WRITE(stdout,*) '  - Already Initialised, Nothing to do here'
                RETURN
            END IF

            ! create the temproraly storage for the data from the generalinfo file containing
            ! the  complit list of available energy level (E) for aveilable inoisation stage (I) of
            ! given element with nuclear charge (Z). Organisation of file as follows:
            !
            ! Nu  Nu  Nu  Nu  Nu  Nu  Nu  Nu  Nu  Nu  Nu
            ! v   |   |   x   x   x   |   x   x   |   x
            ! Line counter C.         |           |
            !     |   |               |           |
            !     v   |               |           v
            !     Nuclear chagre Z.   v           Number of available energy levels E.
            !         v               Number of available lines for I in Z.
            !         Ionisation stage I of Z.

            ! open and read the generalinfo
            OPEN (UNIT=100, FILE=TRIM(DATA_DIR)//'/generalinfo', STATUS='OLD', IOSTAT=IO_status)
            IF (IO_status .NE. 0) THEN
                WRITE(stdout,*) '---- Error Opening the '//TRIM(DATA_DIR)//'/generalinfo file STATUS ', &
                    IO_status, ' -----'
                STOP
            END IF

            DO ind = 1, Lines_genr_info
                READ (100, *) idum, Core_cahrge(ind), Ionis_stage(ind), &
                    idum, idum, idum, idum, idum, idum, Numbr_level(ind), idum
            END DO
            CLOSE (100)
            IF (ver) WRITE(stdout,*) '  - General information read'
            ! generalinfo file read and closed

            ! Create output storage for maximum available E and I for each Z
            OBJ%MaxZ = MAXVAL(Core_cahrge)
            OBJ%MaxI = MAXVAL(Ionis_stage) + 1
            MaxZ = OBJ%MaxZ
            MaxI = OBJ%MaxI

            ALLOCATE (OBJ%List_I(MaxZ))
            OBJ%List_I = 0
            OBJ%MaxL = MAXVAL(Numbr_level)
            MaxL = OBJ%MaxL

            ALLOCATE (OBJ%List_L(MaxI, MaxZ))
            OBJ%List_L = 0

            IF (ver) WRITE(stdout,*) '    > Maximum nuclear charge:  ', OBJ%MaxZ
            IF (ver) WRITE(stdout,*) '    > Maximum ionisation stage:', OBJ%MaxI
            IF (ver) WRITE(stdout,*) '    > Maximum energy level:    ', OBJ%MaxL

            ! Create output storage for the data from the atomnl2_meta_v2.0 file containing the
            ! the lists ionisation energies and degenerecies of every energy level E
            ! for every ionisation state I of spicies with nuclear charge Z.
            ! Orgatisation of file as follows:
            !
            ! non indented lines
            !
            ! Nu  Nu  Nu  Nu     Nu  Nu  Nu
            ! v   |   |   |      |   x   x
            ! Counter corresponding to C.
            !     v   |   |      |
            !     Nuclear charge Z.
            !         v   |      |
            !         Ionisation stage I  of Z.
            !             v      |
            !             Number of energy levels E available for I of Z.
            !                    v
            !                    Ionisation energy of the  gorund lavel Xg (cm^-1).
            !
            ! indented lines
            ! Nu  Nu  Nu  Nu  Nu  St  Nu  Nu  Nu  Nu  Nu  Nu  Nu  Nu
            ! v   x   x   |   |   |   x   |   x   x   |   x   x   x
            ! Counter of the energy level E.          v
            !             v   |   |       |           Ionisation Cross section
            !             Gorund to level E excitation energy Xl (cm^-1).
            !              Ionisation energy of the level E is then -hc*(Xg - Xl).
            !                 v   |       |
            !                 Degeneracy of the level E.
            !                     v       |
            !                     String identifire of
            !                     matastabil 'm' and subordinant's'  levels.
            !                             v
            !                             Transition to ground energy (cm^-1).

            ALLOCATE (OBJ%Degeneracy(MaxL, MaxI, MaxZ))
            OBJ%Degeneracy = 0.0d0
            ALLOCATE (OBJ%Lev_energy(MaxL, MaxI, MaxZ))
            OBJ%Lev_energy = 0.0d0
            ALLOCATE (OBJ%Lev_is_mstbl(MaxL, MaxI, MaxZ))
            OBJ%Lev_is_mstbl = 1
            ALLOCATE (OBJ%Ion_energy(MaxI, MaxZ))
            OBJ%Ion_energy = 0.0d0
            ALLOCATE (OBJ%Zeta(MaxI, MaxZ))
            OBJ%Zeta = 0.0d0

            ! open and read atomnl2_meta_v2.0
            OPEN (UNIT=100, FILE=TRIM(DATA_DIR)//'/atomnl2_meta_v2.0', STATUS='OLD', IOSTAT=IO_status)
            IF (IO_status .NE. 0) THEN
                WRITE(stdout,*) '---- Error Opening the '//TRIM(DATA_DIR)//'/atomnl2_meta_v2.0 file STATUS ', &
                    IO_status, ' -----'
                STOP
            END IF

            DO ind = 1, Lines_genr_info
                ! get info from temproraly lists from Generatlinfo
                Z = Core_cahrge(ind) ! nuclar charge of the element we are reading
                I = Ionis_stage(ind) ! ionisation state of the element
                N_L = Numbr_level(ind) ! number of energy levels available for given spicies
                OBJ%List_L(I, Z) = N_L ! save the number of available energy levels for each spicies

                ! retriev data for the ground state ionisation
                READ (100, *) idum, idum, idum, idum, OBJ%Ion_energy(I, Z), rdum, rdum

                !!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!
                ! IF(Z .EQ. 26) THEN ! fortesting
                !     PRINT*,OBJ%Ion_energy(I,Z)
                ! END IF
                !!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!

                ! contenue reading indented lines for oll the other levels
                DO L = 1, N_L
                    READ (100, *) idum, idum, idum, &
                        OBJ%Lev_energy(L, I, Z), OBJ%Degeneracy(L, I, Z), &
                        Lev_flag, rdum, rdum, rdum, rdum, Ionisation_cs, rdum, rdum, cdum
                    
                    ! Computing the fraction of ground to total recombination 
                    ! processes leading to the ground level by temporarily storeing 
                    ! Sum_l( cross_section(l,i,z)) in Zeta
                    ! NOTE: here we do not destinguish between 'm' and 's' levels 
                    !       for future improvements need to consult Jo
                    OBJ%Zeta(I, Z) = OBJ%Zeta(I, Z) + Ionisation_cs

                    IF (TRIM(Lev_flag) .NE. 'm') THEN
                        OBJ%Lev_is_mstbl(L, I, Z) = 0
                    END IF
                    ! STOP 'stor in read data'
                    ! the data in Lev_energy(E,I,Z) and Degeneracy(L,I,Z) is organized for each level E of
                    ! given ionisation state I for every element Z as
                END DO

                ! Compution Zeta using te previously computed sum 
                ! NOTE: ground level recombination cross section is 1.0
                OBJ%Zeta(I, Z) = 1.0d0 / OBJ%Zeta(I, Z)

            END DO
            CLOSE (100)
            IF (ver) WRITE(stdout,*) '  - Atomic data read '
            ! atomnl2_meta_v2.0 read and closed

            ! open and read the partit file fontaining the degeneracies of 1+
            ! ionisation state of the maximum available ionisation state
            ! file containing the list of Max[I(Z)] and degeneracy for ionisation state Max[I(Z)] + 1

            OPEN (UNIT=100, FILE=TRIM(DATA_DIR)//'/partit', STATUS='OLD', IOSTAT=IO_status)
            IF (IO_status .NE. 0) THEN
                WRITE(stdout,*) '---- Error Opening the '//TRIM(DATA_DIR)//'/partit file STATUS  ', &
                    IO_status, ' -----'
                STOP
            END IF

            READ (100, *) (iduma(ind), ind=1, MaxZ)
            DO Z = 1, MaxZ
                OBJ%List_I(Z) = iduma(Z) + 1
            END DO

            READ (100, *) (iduma(ind), ind=1, MaxZ)
            DO Z = 1, MaxZ
                I = OBJ%List_I(Z)
                OBJ%Degeneracy(1, I, Z) = iduma(Z)
            END DO
            CLOSE (100)
            IF (ver) WRITE(stdout,*) '  - "partit" read'
            ! partit read and closed

            OBJ%Names = ['H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', &
                    & 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', &
                    & 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', &
                    & 'Ni', 'Cu', 'Zn']

            IF (ver) THEN
                WRITE(stdout,*) '  > Elements Elements:'
                WRITE(stdout,*) OBJ%Names(1:MaxZ)
                WRITE(stdout,*) '  - Atomic Data Type READY'
            END IF
        END SUBROUTINE Read_Atom_Data
        !========================================================================

        SUBROUTINE Print_ATOMIC_DATA_TYPE_info(OBJ)
            CLASS(ATOMIC_DATA_TYPE), INTENT(IN) ::  OBJ
            LOGICAL :: Initialised = .FALSE.

            WRITE(stdout,*) ''
            WRITE(stdout,*) 'Information on ATOMIC_DATA_TYPE requested'

            IF (.NOT. ALLOCATED(OBJ%List_I)) THEN
                WRITE(stdout,*) '  - ATOMIC_DATA_TYPE seems NOT Initialised '
                WRITE(stdout,*) '  > For Initialisation CALL ATOMIC_DATA_TYPE_NAME%Initialise()'
                WRITE(stdout,*) '  > Printing ONLY General information'
            ELSE
                WRITE(stdout,*) '  - ATOMIC_DATA_TYPE IS Initialised '
                WRITE(stdout,*) '  > Printing Full information '
                Initialised = .TRUE.
            END IF

            WRITE(stdout,*) ''
            WRITE(stdout,*) '  - Available variables:'
            WRITE(stdout,*) '  > MaxZ, MaxI, MaxL - Maximum number of Z,I,L '
            WRITE(stdout,*) '  > List_I(Z) - Maximum available ionisation stage for element Z '
            WRITE(stdout,*) '  > List_L(I,Z) - Maximum level for ionisation stage I of Z'
            WRITE(stdout,*) '  > Degeneracy(L,I,Z) - Degenerecy of level L of I,Z'
            WRITE(stdout,*) '  > Lev_energy(L,I,Z) - Excitation energy of level L of I,Z'
            WRITE(stdout,*) '  > Lev_is_mstbl(L,I,Z) - Status (1 if level is not m-stable) of level L of I,Z'
            WRITE(stdout,*) '  > Ion_energy(I,Z) - Ionisation energy from ground level of I,Z'
            WRITE(stdout,*) '  > Zeta(I,Z) - fraction of ground to total recombination processes', &
                                ' leading to the ground level of I,Z'
            WRITE(stdout,*) '  > Atomic_weight(Z) - Atomic weight of element Z (output)'
            WRITE(stdout,*) '  > Names(Z) - Name of element Z'
            WRITE(stdout,*) ''
            WRITE(stdout,*) '  - Available PROCEDURE:'
            WRITE(stdout,*) '  > Initialise(verbose) - Initialize the Atomic data, read files and assigne values'
            WRITE(stdout,*) '    > verbose - (OPTIONAL) enable/desable output '
            WRITE(stdout,*) '  > Info() - give Information about the Type/object '

            IF (Initialised) THEN
                WRITE(stdout,*) '  - Allocated data'
                WRITE(stdout,*) '  > MaxZ:', OBJ%MaxZ
                WRITE(stdout,*) '  > MaxI:', OBJ%MaxI
                WRITE(stdout,*) '  > MaxL:', OBJ%MaxL
                WRITE(stdout,*) '  > Available elements:'
                WRITE(stdout,*) OBJ%Names(1:OBJ%MaxZ)
            END IF

            WRITE(stdout,*) 'This PROCEDURE is under constraction'
        END SUBROUTINE Print_ATOMIC_DATA_TYPE_info
        !========================================================================
        !========================================================================
        !========================================================================
        !========================================================================

        SUBROUTINE Read_Composition(OBJ, ATOMIC_DATA, X_frac, Y_frac, verbose)
            ! In-Out and optional variables
            CLASS(NLTE_NUMBER_TYPE), INTENT(INOUT) :: OBJ
            CLASS(ATOMIC_DATA_TYPE), TARGET, INTENT(IN)   :: ATOMIC_DATA
            REAL(DP), OPTIONAL :: X_frac
            REAL(DP), OPTIONAL :: Y_frac
            LOGICAL, OPTIONAL :: verbose
            LOGICAL :: ver, USE_solar_composition, USE_H_normalisation

            REAL(DP), DIMENSION(:), ALLOCATABLE :: Mass_fractions
            REAL(DP), DIMENSION(:), ALLOCATABLE :: Aboundance
            REAL(DP) :: Z_frac = 0.0d0
            INTEGER(I4B) :: MaxZ = 0, MaxI = 0, MaxL = 0, Z = 0, IO_status = 0

            INTEGER(I4B) :: idum
            CHARACTER(LEN=10) :: cdum = '          '

            ! String containing the dectination to the project directiory
            CHARACTER(256) DATA_DIR

            ! Get the destination ot the data location (check if the enviroment was set)
            CALL get_environment_variable("MFORCE_DIR", DATA_DIR)
            IF (LEN_TRIM(DATA_DIR) .EQ. 0) STOP '---- Error: MFORCE_DIR not set -----'
            DATA_DIR = TRIM(DATA_DIR)//'/DATA'

            ! Pars the input optianal variables
            ! If verbose = .true. specified give datailed status
            ver = .FALSE.
            IF (PRESENT(verbose)) ver = verbose

            ! create the linc to atomic data for the later use
            OBJ%ATOMIC_DATA => ATOMIC_DATA

            IF (ver) WRITE(stdout,*) ' '
            IF (ver) WRITE(stdout,*) 'Initialise LTE Number Type Called'
            IF (ver) WRITE(stdout,*) '  - Initialising variables and setting composition'

            ! Chech if ATOMIC_DATA is allocated - if not Stop
            IF (.NOT. ALLOCATED(ATOMIC_DATA%List_I)) STOP '---- Terminating: Atomic data not Initialised ----'

            ! create short-hands
            MaxZ = ATOMIC_DATA%MaxZ
            MaxI = ATOMIC_DATA%MaxI
            MaxL = ATOMIC_DATA%MaxL

            ! Chech if LTE_type is already Initialised
            IF (ALLOCATED(OBJ%Aboundance)) THEN
                IF (ver) WRITE(stdout,*) '  - Already Initialised, Nothing to do here'
                RETURN
            END IF

            ! Check the input X_frac + Y_frac < 1 and Y_frac |= 0
            IF ((PRESENT(X_frac) .AND. PRESENT(Y_frac))) THEN
                IF ((X_frac + Y_frac) .LE. 1.0d0) THEN
                    IF (ver) WRITE(stdout,*) '  - Input mass fractions provided'
                    USE_solar_composition = .FALSE.
                    IF ( ABS(X_frac - 0.0d0) .GT. EPSILON(X_frac)) THEN
                        IF (ver) WRITE(stdout,*) '  - Scaling aboundances to H'
                        USE_H_normalisation = .TRUE.
                    ELSE
                        IF (ver) WRITE(stdout,*) '  - Scaling aboundances to He'
                        USE_H_normalisation = .FALSE.
                    END IF
                ELSE
                    STOP '---- Error in input X and Y mass fractions (X+Y > 1 or Y = 0) ----'
                END IF
            ELSE
                USE_solar_composition = .TRUE.
                USE_H_normalisation = .TRUE.
                IF (ver) WRITE(stdout,*) '  - No input mass fractions provided: Using solar values'
            END IF

            ! Allocate and zero out the occupation number arraies
            ALLOCATE (OBJ%Aboundance(MaxZ))
            OBJ%Aboundance = 0.0d0
            ALLOCATE (OBJ%Partition(MaxI, MaxZ))
            OBJ%Partition = 0.0d0
            ALLOCATE (OBJ%Ions(MaxI, MaxZ))
            OBJ%Ions = 0.0d0
            ALLOCATE (OBJ%Nuclei(MaxZ))
            OBJ%Nuclei = 0.0d0
            ALLOCATE (OBJ%Occupation(MaxL, MaxI, MaxZ))
            OBJ%Occupation = 0.0d0

            ! preper remaining variables
            OBJ%Electron = 0.0d0
            OBJ%X_frac = 0.0d0
            OBJ%Y_frac = 0.0d0
            OBJ%T = 0.0d0
            OBJ%rho = 0.0d0

            ! Allocate temp variables
            ALLOCATE (Mass_fractions(MaxZ))
            Mass_fractions = 0.0d0
            ALLOCATE (Aboundance(MaxZ))
            Aboundance = 0.0d0

            ! Clear the rho and T by settiong them to negative values to net
            OBJ%rho = -1.d0
            OBJ%T = -1.d0

            ! Open the line list info file
            OPEN (UNIT=100, FILE=TRIM(DATA_DIR)//'/abundan_solar_asplund_2009', STATUS='OLD', IOSTAT=IO_status)
            ! Check if opening was sucssesful
            IF (IO_status .NE. 0) THEN
                WRITE(stdout,*) '---- Error Opening the '//TRIM(DATA_DIR)//'/abundan_solar_asplund_2009 file STATUS ', &
                    IO_status, ' ----'
                STOP
            END IF

            DO Z = 1, MaxZ
                READ (100, *) idum, cdum, OBJ%Aboundance(Z)
            END DO
            CLOSE (100)
            IF (ver) WRITE(stdout,*) '  - Aboundances READ'

            ! Poniatowski+ A.1
            ! here we normalise abandences to H or He depending if there is Y aboundance and compute mass fractions for solar mixture.
            IF (USE_H_normalisation) THEN
                OBJ%Aboundance(:) = 10.0d0**(OBJ%Aboundance(:) - OBJ%Aboundance(1))
            ELSE
                OBJ%Aboundance(:) = 10.0d0**(OBJ%Aboundance(:) - OBJ%Aboundance(2))
            END IF

            Mass_fractions(:) = OBJ%Aboundance(:)*ATOMIC_DATA%Atomic_weight(:) ! compute direct masses
            Mass_fractions(:) = Mass_fractions(:)/SUM(Mass_fractions(:)) ! normalize to compute mass fractions

            ! If there is the  input X and Y fractions presentthan
            ! here we set the Mass fraction of H and He to input X and Y
            ! and renormilize metal mass fractions of solar mixture to agree with input X Y.
            IF (.NOT. USE_solar_composition) THEN
                Mass_fractions(1) = X_frac ! set the H mass fraction  to X
                Mass_fractions(2) = Y_frac ! set the He mass fraction to Y
                ! normalaze the summed mass fraction of the metals so that
                ! resulting Z = 1 - X - Y
                Z_frac = MAX(0.0d0, 1.0d0 - Mass_fractions(1) - Mass_fractions(2))
                IF (ver) WRITE(stdout,*) '  > Z fraction ', Z_frac
                Mass_fractions(3:MaxZ) = Mass_fractions(3:MaxZ)/SUM(Mass_fractions(3:MaxZ))*Z_frac

                ! Allocate the input X and Y mass fractions
                OBJ%X_frac = X_frac
                OBJ%Y_frac = Y_frac

                ! update the number dencities to with new mass fractions
                OBJ%Aboundance(:) = Mass_fractions(:)/ATOMIC_DATA%Atomic_weight(:)
                ! here we normalise abandences back to H or He depending if there is Y aboundance and compute mass fractions for solar mixture.
                IF (USE_H_normalisation) THEN
                    OBJ%Aboundance(:) = OBJ%Aboundance(:)/OBJ%Aboundance(1)
                ELSE
                    OBJ%Aboundance(:) = OBJ%Aboundance(:)/OBJ%Aboundance(2)
                END IF
            ELSE
                OBJ%X_frac = Mass_fractions(1)
                OBJ%Y_frac = Mass_fractions(2)
            END IF

            IF (ver) WRITE(stdout,*) '  - LTE Number Type READY'
        END SUBROUTINE Read_Composition
        !========================================================================

        SUBROUTINE Comput_NLTE_Occupation_Numbers(OBJ, rho, T, Te_to_T, Dilution, N_electron, verbose)
            CLASS(NLTE_NUMBER_TYPE), INTENT(INOUT) :: OBJ
            REAL(DP), INTENT(IN) :: rho
            REAL(DP), INTENT(IN) :: T
            REAL(DP), OPTIONAL :: Te_to_T
            REAL(DP), OPTIONAL :: Dilution
            REAL(DP), OPTIONAL :: N_electron
            LOGICAL, OPTIONAL :: verbose
            CLASS(ATOMIC_DATA_TYPE), POINTER :: ATOMIC_DATA  

            ! Temprorary variables
            LOGICAL  :: ver, iterative, use_nlte
            REAL(DP) :: N_electron_old = 0.0d0, N_electron_new = 0.0d0
            REAL(DP) :: Te, W
            INTEGER(I4B) :: I = 0, Z = 0, ind = 0

            ! Pars the input optianal variables
            ! If trerative = true recompute the electrone number density
            iterative = .TRUE.
            IF (PRESENT(N_electron)) iterative = .FALSE.
            ! If verbose = .true. specified give datailed status
            ver = .FALSE.
            IF (PRESENT(verbose)) ver = verbose
            ! If Tr to Te  ratio is present use NLTE option, else set Tr = T
            IF (PRESENT(Te_to_T)) THEN 
                use_nlte = .TRUE.
                Te = T * Te_to_T
            ELSE
                use_nlte = .FALSE.
                Te = T
            ENDIF

            IF (PRESENT(Dilution)) THEN 
                W  = Dilution
            ELSE
                W = 1.0d0
            ENDIF

            IF (ver) WRITE(stdout,*) ' '
            IF (ver) WRITE(stdout,*) 'Set LTE Temperature and Dencity CALLED'
            IF (ver) WRITE(stdout,*) '  - Computing Occupation numbers'
            IF (ver.AND.use_nlte) WRITE(stdout,*) '  > Using option NLTE'

            ! Chech if the NLTE_NUMBER_TYPE is allocated
            IF (.NOT. ALLOCATED(OBJ%Aboundance)) STOP '---- Terminating: LTE Numbre Type not Initialised ----'

            ! Check input dencity and temperature
            IF (.NOT. (rho .GT. 0.0d0)) STOP '---- Terminating: Inappropriate Dencity input ---- '
            IF (.NOT. (T .GT. 0.0d0)) STOP '---- Terminating: Inappropriate Temperature input ---- '
            If (.NOT.(Te .GT. 0.0d0)) STOP '---- Terminating: Inappropriate Temperature retio input ---- '
            If (.NOT.(W .GT. 0.0d0)) STOP '---- Terminating: Inappropriate Dilution factor input ---- '

            ! create the assosiation(Short-hand) to use the previous link to atomic data
            ATOMIC_DATA => OBJ%ATOMIC_DATA

            ! Store the inpude ,ass density and temperature
            OBJ%rho = rho
            OBJ%T = T
            OBJ%Te = Te
            OBJ%W = W

            ! To compute the Occupation numbers do the following:
            ! -- if the NLTE optimization is used skip steps marked with * --
            ! 1*) For  given Rho(D) and Aboundance Compute the numper density of each spicies
            ! 2*) If not given Compute the approximate numberdencity of electronse assuming fully ionised H and 1/2 rest
            ! 3*) For a given T_rad(T) and Atomic Data compute partition functions for each ion
            ! 4) Using the Atomic Data and partition functions compute the Ion number dencity (solve SAHA)
            ! 5*) If N_electon was not specified iterate 4) to compute the N_electon
            ! 6a*) Using the Atomic Data and Ion number dencity compute the occupation numbers (solve BOLTZMANN)
            ! 6b*) Sum up occupation numbers to check the mass conservation

            ! ----- Executing 1) !!!!
            ! compute the output total number densities for each spicies H or He depending if there is Y aboundance
            ! Poniatowski+ A.2
            IF (ABS( OBJ%X_frac  - 0.0d0) .GT. EPSILON(OBJ%X_frac) ) THEN
                IF (ver) WRITE(stdout,*) '  - Using the X mass fraction for computations'
                OBJ%Nuclei(:) = rho*OBJ%Aboundance(:)*OBJ%X_frac/OBJ%mH
            ELSE
                IF (ver) WRITE(stdout,*) '  - Using the Y mass fraction for computations'
                OBJ%Nuclei(:) = rho*OBJ%Aboundance(:)*OBJ%Y_frac/(4.0d0*OBJ%mH)
            END IF

            !!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!
            ! DO Z = 1,30
            !    PRINT*,OBJ%Nuclei(Z)
            ! END DO
            !!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!

            IF (ANY(ISNAN(OBJ%Nuclei(:)))) STOP '---- Nuclei Is NaN ----'
            IF (ver) WRITE(stdout,*) '  - Number densities of each Nuclei computed'

            ! ------ Executing 2)
            ! Here if N_electon was not given compute the Approximate
            IF (iterative) THEN
                ! here compute the initial guess of electron number dencity by assuming
                ! that hydrogen is fully ionised and elements form He and up contributes 1/2 of their electrons
                OBJ%Electron = OBJ%Nuclei(1)
                DO Z = 2, ATOMIC_DATA%MaxZ
                    OBJ%Electron = OBJ%Electron + Z/2.0d0*OBJ%Nuclei(Z)
                END DO

                IF (ver) WRITE(stdout,*) '   - No Input Electron Number Dencity found'
                IF (ver) WRITE(stdout,*) '   > Using iterative scheme, with initial guess:', OBJ%Electron
            ELSE
                OBJ%Electron = N_electron
                IF (ver) WRITE(stdout,*) '   - Input Electron Number Dencity found', OBJ%Electron
            END IF

            ! ----- Executing 3)
            CALL Comput_Partition_Functions(OBJ, ATOMIC_DATA, T, verbose=ver)
            IF (ANY(ISNAN(OBJ%Partition))) STOP '---- Partition functuon Is NaN ----'
            IF (ver) WRITE(stdout,*) '  - Partition Functions READY'

            ! prepare the old and new ststes for itteration if iterative option used 
            N_electron_old = OBJ%Electron
            N_electron_new = OBJ%Electron

            ! ----- Executing 4 % 5)
            ! If input electrone dencity was not specified do iteration
            IF (iterative) THEN

                DO ind = 1, 100 ! set the maximum number of iterations to 100
                    
                    ! update the electrone numbe dencities by Ne_new = sqrt(Ne_new * Ne_old)
                    ! see Lattimer & Cranmer 2020
                    OBJ%Electron = SQRT(N_electron_old*N_electron_new)
                    N_electron_old = OBJ%Electron

                    CALL Comput_Ionisation_structure(OBJ, ATOMIC_DATA, T, Te, W, iterative, verbose=ver)
                    IF (ISNAN(OBJ%Electron) .OR. &
                    ( ABS(OBJ%Electron - 0.0d0) .LE. EPSILON(OBJ%Electron))) THEN
                        WRITE(stdout,*) OBJ%Electron, ABS(OBJ%Electron - 0.0d0), EPSILON(OBJ%Electron)
                        STOP '---- Electrons Is NaN/Zero ----'
                    END IF
                    N_electron_new = OBJ%Electron

                    ! chech if the change in the numberdencity is less then 0.1%
                    ! is yes Ne is converged - exit the loop -- otherwithe continue
                    IF (ABS(N_electron_old/N_electron_new - 1) .LT. 1.0d-3) THEN
                        IF (ver) WRITE(stdout,*) '   > Electron Number Dencity converged', OBJ%Electron
                        EXIT
                    ELSE
                        IF (ver) WRITE(stdout,*) '   > Iteration:', ind, ' Electron Number Dencity:', OBJ%Electron
                    END IF
                END DO

                ! if the maximum number of iterations was reached terminate
                IF (ind .GT. 99) STOP '---- Electron Number Dencity NOT converged ----'
            ELSE

                ! compute the ionisation structure using the input electrone number density
                CALL Comput_Ionisation_structure(OBJ, ATOMIC_DATA, T, Te, W, iterative, verbose=ver)
            END IF

            !!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!
            ! DO Z = 1,30                                  
            !   PRINT*, Z
            !   DO I = 1,ATOMIC_DATA%List_I(Z)
            !    PRINT*,I, OBJ%Ions(I,Z)
            !  END DO
            ! END DO
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! Sanity Check
            IF (ANY(ISNAN(OBJ%Ions(:, :)))) STOP '---- Ions Is NaN ----'
            IF (ver) WRITE(stdout,*) '  - Ionisation Fructions READY'

            ! ----- Executing 6a)

            ! compute the occupation number solving Boltzmann + modified
            IF (ver) WRITE(stdout,*) '   - Computing the occupation numbers'
            CALL Comput_Excitation_structure(OBJ, ATOMIC_DATA, T, W, verbose=ver)

            ! ----- Executing 6b)
            ! cheking the total number dencity
            DO Z = 1, ATOMIC_DATA%MaxZ
                DO I = 1, ATOMIC_DATA%List_I(Z) - 1
                    IF (SUM(OBJ%Occupation(:, I, Z))/OBJ%Ions(I, Z) - 1 .GT. 1.0d-6) THEN
                        STOP '---- Number Dencities NOT conserved ----'
                    END IF
                END DO
            END DO

            IF (ver) WRITE(stdout,*) '  - NLTE Numbers READY'
        END SUBROUTINE Comput_NLTE_Occupation_Numbers
        !========================================================================

        SUBROUTINE Comput_Partition_Functions(OBJ, ATOMIC_DATA, T, verbose)
            CLASS(NLTE_NUMBER_TYPE), INTENT(INOUT) :: OBJ
            CLASS(ATOMIC_DATA_TYPE), INTENT(IN)    :: ATOMIC_DATA
            REAL(DP), INTENT(IN) :: T
            LOGICAL, OPTIONAL :: verbose

            ! Temprorary variables
            LOGICAL :: ver
            INTEGER(I4B) :: Z = 0, I = 0, L = 0
            REAL(DP) :: Bolz_fact = 0.0d0, Xg = 0.0d0, gg = 0.0d0 !, Xl = 0.0d0, gl = 0.0d0

            ! Pars the input optianal variables
            ! If verbose = .true. specified give datailed status
            ver = .FALSE.
            IF (PRESENT(verbose)) ver = verbose

            IF (ver) WRITE(stdout,*) '  - Computing Partition Functions'
            IF (ver) WRITE(stdout,*) '   - Partition fanctions of Fe:'
            ! compute the partition FUNCTION of a single ionisation state I of a single
            ! atom Z
            ! Partition_functuon = zeros(MaxI,MaxZ) is replaces by OBJ%Partition

            ! Poniatowski+ A.3
            ! For every available nuclear charge Z from 1 to Z_max
            DO Z = 1, ATOMIC_DATA%MaxZ ! Do for each Z

                ! % For every available ionization stage I of Z
                DO I = 1, ATOMIC_DATA%List_I(Z) ! Do for each Ionisation

                    Xg = ATOMIC_DATA%Lev_energy(1, I, Z) ! energy of the ground level
                    gg = ATOMIC_DATA%Degeneracy(1, I, Z) ! degenerecy of the ground level
                    OBJ%Partition(I, Z) = gg ! Partition FUNCTION of (I,Z) - ground level contributes degenerecy

                    ! For energy level E of Z WITH I
                    DO L = 2, ATOMIC_DATA%List_L(I, Z) ! Do for each level

                        ! Xl = ATOMIC_DATA%Lev_energy(L, I, Z) ! excitation energy of level E
                        ! gl = ATOMIC_DATA%Degeneracy(L, I, Z) ! degenerecy of the level E

                        ! Boltzmann factor EXP( -x/kT)
                        Bolz_fact = EXP(-OBJ%Bolz_norm*(ATOMIC_DATA%Lev_energy(L, I, Z) - Xg)/T) 
                        ! Partition FUNCTION of (I,Z)
                        OBJ%Partition(I, Z) = OBJ%Partition(I, Z) + &
                                            ATOMIC_DATA%Degeneracy(L, I, Z)*Bolz_fact 

                    END DO ! Do for each level

                    IF (ver .AND. (Z .EQ. 26)) WRITE(stdout,*) '   > ', I, OBJ%Partition(I, Z)
                END DO ! Do for each Ionisation
            END DO! Do for each Z
        END SUBROUTINE Comput_Partition_Functions
        !========================================================================

        SUBROUTINE Comput_Ionisation_structure(OBJ, ATOMIC_DATA, T, Te, W, iterative, verbose)
            CLASS(NLTE_NUMBER_TYPE), INTENT(INOUT) :: OBJ
            CLASS(ATOMIC_DATA_TYPE), INTENT(IN)    :: ATOMIC_DATA
            REAL(DP), INTENT(IN) :: T
            REAL(DP), INTENT(IN) :: Te
            REAL(DP), INTENT(IN) :: W
            LOGICAL, OPTIONAL :: verbose

            ! Temprorary variables
            LOGICAL :: iterative
            LOGICAL :: ver
            INTEGER(I4B) :: Z = 0, I = 0, k = 0
            REAL(DP) :: Saha_fact = 0.0d0, NLTE_fact = 0.0d0 !, Xi = 0.0d0
            ! REAL(DP) :: Upart1 = 0.0d0, Upart2 = 0.0d0
            REAL(DP) :: SUM_Fractions_to_nutral = 0.0d0
            REAL(DP), DIMENSION(:, :), ALLOCATABLE :: S
            REAL(DP), DIMENSION(:, :), ALLOCATABLE :: Fractions_to_nutral
            REAL(DP), DIMENSION(:, :), ALLOCATABLE :: Number_frac

            ! Pars the input optianal variables
            ver = .FALSE.
            IF (PRESENT(verbose)) ver = verbose

            IF (ver) WRITE(stdout,*) '   - Computing Ionisation Fractions'

            ALLOCATE(S(ATOMIC_DATA%MaxI, ATOMIC_DATA%MaxZ))
            S = 0.0d0

            ! For every available nuclear charge Z from 1 to Z_max
            ! compute the Saha fraction n_e * n(I+1)/n(I) and store in S(I,Z)
            DO Z = 1, ATOMIC_DATA%MaxZ
                !!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!
                ! Z = 26   ! for testing
                !!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!

                ! For every available ionization stage I of Z
                DO I = 1, ATOMIC_DATA%List_I(Z) - 1

                    ! ! partition function of the ionization state I 
                    ! Upart1 = OBJ%Partition(I, Z) 
                    ! ! partition function of the ionization state I+1 
                    ! Upart2 = OBJ%Partition(I + 1, Z) 
                    ! ! ion I ionization energy from the gorund level
                    ! Xi = ATOMIC_DATA%Ion_energy(I, Z) 

                    ! Saha exponentia factor
                    Saha_fact = EXP( -OBJ%Bolz_norm* ATOMIC_DATA%Ion_energy(I, Z)/T )
                    
                    ! NLTE cottection factor Lucy & Abbott 1993;  Puls et al. 2000 
                    NLTE_fact = W * SQRT(Te / T) * (ATOMIC_DATA%Zeta(I, Z) + W * (1.0d0 - ATOMIC_DATA%Zeta(I, Z)) )

                    ! Poniatowski+ A.4
                    ! S(I, Z) = 2.0d0 * Upart2 / Upart1 * OBJ%Saha_norm *T**1.5d0 * Saha_fact * NLTE_fact 
                    S(I, Z) = 2.0d0 * OBJ%Partition(I + 1, Z) / OBJ%Partition(I, Z) * &
                            OBJ%Saha_norm *T**1.5d0 * Saha_fact * NLTE_fact 
                            ! Saha fraction where S(I) = n_e * n(I+1)/n(I)

                END DO! For every available ionization stage
            END DO ! For every  nucleus

            IF (ver) THEN
                WRITE(stdout,*) '    - The S(I) = n_e * n(I+1)/n(I) of Fe:'
                DO I = 1, ATOMIC_DATA%List_I(26)
                    WRITE(stdout,*) '     > ', I, S(I, 26)
                END DO
            END IF

            !!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!
            ! % psi = log(S/n_elect); % is equvalent to psi in fastwind
            !!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!

            ! initiate temprorary memory
            ALLOCATE (Fractions_to_nutral(ATOMIC_DATA%MaxI, ATOMIC_DATA%MaxZ)) ! ionisation fractions with respect to neutral Ion n_{I}/n_{1} of Z
            Fractions_to_nutral = 0.0d0
            ALLOCATE (Number_frac(ATOMIC_DATA%MaxI, ATOMIC_DATA%MaxZ)) ! output number fraction of ions in with respect to the total number
            Number_frac = 0.0d0
            !  Here we conpute the number density of each ionisation state for everiy spicies Z by solving:
            ! \[ n_{z,\,1}\left(1 + \sum\limits_{k=1}^{I_{max}} \prod\limits_{i=1}^{k}
            !  \frac{S_{z,\,i}}{n_e}\right) = N_{z,\,tot} \]
            !  For that we first compute the fractions on ion number dencity
            !  with respec to nutral ions as:
            !  \[f_{z,\, i} = \frac{ n_{z,\, i} }{ n_{z,\, 1} } \Rightarrow
            !  \prod\limits_{i=1}^{I} \frac{S_{z,\,i}}{n_e} \]
            !  Then find number density of neutral ions
            ! \[ n_{z,\,1} = \frac{ N_{z,\,tot} }{ \sum\limits_{k=1}^{I_{max}}f_{z,\,i} }\]
            !  And finally we compute number density for each ionisation state as:
            ! \[ n_{z,\,i} = n_{z,\,1} f_{z,\,i} \]
            !  see the notes for definition

            ! cycle through all the spicies
            DO Z = 1, ATOMIC_DATA%MaxZ
                !  Z = 26;

                ! set n1/n1 = e^0 = 1
                Fractions_to_nutral(1, Z) = 1.0d0

                SUM_Fractions_to_nutral = 1.0d0

                ! do the sum over every ionisation state
                DO I = 2, ATOMIC_DATA%List_I(Z)
                    Fractions_to_nutral(I, Z) = 1.0d0

                    ! Poniatowski+ A.5
                    ! do the multiplication for all the
                    DO k = 1, I - 1
                        Fractions_to_nutral(I, Z) = Fractions_to_nutral(I, Z) * S(k, Z) /OBJ%Electron
                        ! here the Fractions_to_nutral(i,z) is equivalent of
                        ! exp(prod(i)) form Fastwind and SUM_Fractions_to_nutral
                        ! = summ for each atom
                    END DO ! do the multiplication

                    ! Poniatowski+ A.6 LHS brackets term
                    SUM_Fractions_to_nutral = SUM_Fractions_to_nutral + Fractions_to_nutral(I, Z) 
                END DO ! do the  sum 

                ! Poniatowski+ A.6
                Number_frac(:, Z) = Fractions_to_nutral(:, Z) / SUM_Fractions_to_nutral
                ! here Number_frac(:,Z) is equicalent to fkl(Z,I,?)  in fastwind
                OBJ%Ions(:, Z) = OBJ%Nuclei(Z) * Number_frac(:, Z)

            END DO ! cycle through all the spicies

            ! using the new atomic database Problem shows up in this routine when checking for the number dencity conservation 
            ! one can check thet using this bit of code:
            ! DO Z = 1,ATOMIC_DATA%MaxZ
            !    PRINT*, 'err = ' , ABS( SUM(OBJ%Ions(1:ATOMIC_DATA%List_I(Z),Z))/OBJ%Nuclei(Z) - 1 )
            !    IF( ABS( SUM(OBJ%Ions(1:ATOMIC_DATA%List_I(Z),Z))/OBJ%Nuclei(Z) - 1 ) .GT. 1.0d-6) THEN
            !       STOP '---- Number Densities in ionisation NOT conserved ----'
            !    END IF
            ! END DO
            ! As here we only compute relative number fractions
            ! Solution will be to multipy all of the fraction for each Z by the relative error in the number dencity
             
            IF (ver) THEN
                WRITE(stdout,*) '    - Ion number fractions relative to ground and ino densities of Fe'
                DO I = 1, ATOMIC_DATA%List_I(26)
                    WRITE(stdout,*) '     > ', I, Fractions_to_nutral(I, 26), OBJ%Ions(I,26)
                END DO
            END IF

            IF (iterative) THEN
                ! here we compute the total numberdencity of electrons by computing number of electrons given by each ion
                OBJ%Electron = 0
        
                ! Poniatowski+ A.7
                DO Z = 1, ATOMIC_DATA%MaxZ
                    DO I = 1, ATOMIC_DATA%List_I(Z)
                        OBJ%Electron = OBJ%Electron + OBJ%Ions(I, Z)*(I - 1)
                    END DO
                END DO
            END IF
            IF (ISNAN(OBJ%Electron)) STOP '---- Error With Electron number density (NAN) ----'
        END SUBROUTINE Comput_Ionisation_structure
        !========================================================================

        SUBROUTINE Comput_Excitation_structure(OBJ, ATOMIC_DATA, T, W, verbose)
            CLASS(NLTE_NUMBER_TYPE), INTENT(INOUT) :: OBJ
            CLASS(ATOMIC_DATA_TYPE), INTENT(IN)    :: ATOMIC_DATA
            REAL(DP), INTENT(IN) :: T
            REAL(DP), INTENT(IN) :: W
            LOGICAL, OPTIONAL :: verbose

            ! Temprorary variables
            LOGICAL :: ver
            REAL(DP), DIMENSION(:), ALLOCATABLE :: Excitation_numb_frac
            REAL(DP) :: gg = 0.0d0, Xg = 0.0d0 !, gl = 0.0d0, Xl = 0.0d0
            REAL(DP) :: Bolz_fact = 0.0d0,  SUM_Excitation_numb_frac = 0.0d0
            INTEGER(I4B) :: L = 0, I = 0, Z = 0

            ! Pars the input optianal variables
            ! If verbose = .true. specified give datailed status
            ver = .FALSE.
            IF (PRESENT(verbose)) ver = verbose

            ALLOCATE (Excitation_numb_frac(ATOMIC_DATA%MaxL))

            !!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!
            ! Z = 1; I = 1;
            !!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!
            DO Z = 1, ATOMIC_DATA%MaxZ ! for each element

                DO I = 1, ATOMIC_DATA%List_I(Z) - 1 ! for each ion of element

                    Excitation_numb_frac = 0.0d0
                    Excitation_numb_frac(1) = 1.0d0 ! set the fraction of groung state to 1
                    gg = ATOMIC_DATA%Degeneracy(1, I, Z) ! copy the statistical weight of ground level
                    Xg = ATOMIC_DATA%Lev_energy(1, I, Z) ! copy the excitation energy of ground level

                    SUM_Excitation_numb_frac = 1.0d0

                    DO L = 2, ATOMIC_DATA%List_L(I, Z) ! for every level higher then gorund

                        ! ! copy the statistical weight of the level
                        ! gl = ATOMIC_DATA%Degeneracy(L, I, Z) 
                        ! ! copy the excitation energy of the level
                        ! Xl = ATOMIC_DATA%Lev_energy(L, I, Z) 

                        ! compute the Boltzman factor
                        ! Bolz_fact = EXP(-OBJ%Bolz_norm*(Xl - Xg)/T) 
                        Bolz_fact = EXP(-OBJ%Bolz_norm*(ATOMIC_DATA%Lev_energy(L, I, Z) - Xg)/T) 

                        ! compute the number fraction of the level to ground
                        ! Accounting for the statis of the level
                        ! If level is M-stable then nl/n1 = boltz(T), else nl/n1 = W*boltz(T)
                        ! Excitation_numb_frac(L) = gl/gg*Bolz_fact
                        Excitation_numb_frac(L) = ATOMIC_DATA%Degeneracy(L,I,Z)/gg*Bolz_fact * &
                                                (W + (1.0d0-W)*ATOMIC_DATA%Lev_is_mstbl(L,I,Z))

                        ! ------------------------------------------------------
                        ! for some reason the SUM(Excitation_numb_frac) becoes nan
                        ! otherwithe here we use the dummy_SUM
                        ! ------------------------------------------------------
                        SUM_Excitation_numb_frac = SUM_Excitation_numb_frac + Excitation_numb_frac(L)

                        IF (ISNAN(Excitation_numb_frac(L))) STOP '---- Ef is nan ----'
                        IF (ISNAN(SUM_Excitation_numb_frac)) STOP '---- SUM Ef is nan ----'
                    END DO

                    !!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!
                    ! IF(ISNAN(SUM(Excitation_numb_frac)))  THEN
                    !   PRINT*,SUM(Excitation_numb_frac),dummy_SUM
                    !    DO L = 1,ATOMIC_DATA%List_L(I,Z)
                    !       PRINT*,L,Excitation_numb_frac(L)
                    !    END DO
                    !    STOP '---- Sum(Ef) is NaN ----'
                    ! END IF
                    ! IF(SUM(Excitation_numb_frac).EQ.0.0d0) &
                    !      STOP '---- Sum(Ef) is Zero ----'
                    !!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!


                    ! Compute the occupation number of each level using the number density of ion and fractions
                    OBJ%Occupation(1:ATOMIC_DATA%List_L(I, Z), I, Z) = &
                        OBJ%Ions(I, Z)*Excitation_numb_frac(1:ATOMIC_DATA%List_L(I, Z))/SUM_Excitation_numb_frac

                    IF (ver .AND. (Z .EQ. 26) .AND. (I .EQ. 1)) THEN
                        WRITE(stdout,*) '   - Number fractions of Fe I relative to ground '
                        DO L = 1, ATOMIC_DATA%List_L(1, 26)
                            WRITE(stdout,*) '    > ', L, OBJ%Occupation(L, 1, 26)
                        END DO
                    END IF

                END DO! for each ion of element
            END DO! for each element

            ! here the Occup_Number_denc(:,I,Z)/Atomic_data{:,I,Z}(2) should be
            ! equivalent of occng(i,:,radius), where i is record number from gen.info.
            
            !!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!
            ! for testing run:
            ! Z = 26; I = 1;
            !!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!
            IF (ver .AND. (Z .EQ. 26)) THEN
                WRITE(stdout,*) '   - Occupation number densities of Fe I '
                DO L = 1, ATOMIC_DATA%List_L(1, 26)
                    WRITE(stdout,*) '    > ', L, OBJ%Occupation(L, 1, 26)
                END DO
            END IF

        END SUBROUTINE Comput_Excitation_structure
        !========================================================================

        SUBROUTINE Clear_NLTE_NUMBER_TYPE(OBJ)
            CLASS(NLTE_NUMBER_TYPE), INTENT(INOUT) :: OBJ
            CLASS(ATOMIC_DATA_TYPE), POINTER    :: ATOMIC_DATA

            ! create the assosiation(Short-hand) to use the previous link to atomic data
            ATOMIC_DATA => OBJ%ATOMIC_DATA

            DEALLOCATE (OBJ%Partition)
            DEALLOCATE (OBJ%Ions)
            DEALLOCATE (OBJ%Nuclei)
            DEALLOCATE (OBJ%Occupation)

            ALLOCATE (OBJ%Partition(ATOMIC_DATA%MaxI, ATOMIC_DATA%MaxZ))
            OBJ%Partition = 0.0d0
            ALLOCATE (OBJ%Ions(ATOMIC_DATA%MaxI, ATOMIC_DATA%MaxZ))
            OBJ%Ions = 0.0d0
            ALLOCATE (OBJ%Nuclei(ATOMIC_DATA%MaxZ))
            OBJ%Nuclei = 0.0d0
            ALLOCATE (OBJ%Occupation(ATOMIC_DATA%MaxL, ATOMIC_DATA%MaxI, ATOMIC_DATA%MaxZ))
            OBJ%Occupation = 0.0d0

            ! Number dencity of electrones
            OBJ%Electron = 0.0d0

            ! Radiation+Equilibrium Temperature assigned to the type
            OBJ%T = -1.0d0
            ! Electorn Temperature assigned to the type
            OBJ%Te = -1.0d0
            ! Dilution factor assigned to the type
            OBJ%W = 1.0d0
            ! Mass dencity assigned to the type
            OBJ%rho = -1.0d0
        END SUBROUTINE Clear_NLTE_NUMBER_TYPE
        !========================================================================

        SUBROUTINE Print_NLTE_NUMBER_TYPE_info(OBJ)
            CLASS(NLTE_NUMBER_TYPE), INTENT(IN) ::  OBJ
            LOGICAL :: Initialised = .FALSE.
            LOGICAL :: Set = .FALSE.

            WRITE(stdout,*) ' '
            WRITE(stdout,*) 'Information on NLTE_NUMBER_TYPE requested'

            IF (.NOT. ALLOCATED(OBJ%Aboundance)) THEN
                WRITE(stdout,*) '  - NLTE_NUMBER_TYPE seems NOT Initialised '
                WRITE(stdout,*) '  > For Initialisation CALL NLTE_NUMBER_TYPE_NAME%Initialise()'
                WRITE(stdout,*) '  > Printing ONLY General information'
            ELSE
                WRITE(stdout,*) '  - NLTE_NUMBER_TYPE IS Initialised '
                Initialised = .TRUE.
                IF (OBJ%rho .GT. 0) THEN
                    WRITE(stdout,*) '  - NLTE_NUMBER_TYPE Temperature and Density are Set'
                    WRITE(stdout,*) '  > Printing Full information '
                    Set = .TRUE.
                ELSE
                    WRITE(stdout,*) '  - NLTE_NUMBER_TYPE Temperature and Density are NOT Set'
                    WRITE(stdout,*) '  > For Set T and rho CALL NLTE_NUMBER_TYPE_NAME%Set()'
                    WRITE(stdout,*) '  > Printing ONLY General information'
                END IF
            END IF

            WRITE(stdout,*) ''
            WRITE(stdout,*) '  - Available variables:'
            WRITE(stdout,*) '  > Aboundance(z) - Number aboundance of element Z'
            WRITE(stdout,*) '  > Nuclei(Z) - Number density of elemet Z'
            WRITE(stdout,*) '  > Ions(I,Z) - Number density of nucleis Z in ionisation stage I'
            WRITE(stdout,*) '  > Occupation(L,I,Z) - Number density of Z I in excitation state L'
            WRITE(stdout,*) '  > Electron - Number density of electrons (input see Set()/output)'
            WRITE(stdout,*) '  > X_frac - Hydrogen mass fraction (input see Initialise()/output)'
            WRITE(stdout,*) '  > Y_frac - Helium mass fraction (input see Initialise()/output)'
            WRITE(stdout,*) '  > T - Input radiatio(equilibrium) temperature in K (input)'
            WRITE(stdout,*) '  > Te_to_T - Input ratio of electron to radiation temperatures (OPTIONAL)'
            WRITE(stdout,*) '  > W - Input Dilution factor (OPTIONAL)'
            WRITE(stdout,*) '  > rho - Input mass dencity in CGS (input)'
            WRITE(stdout,*) '  > Bolz_norm - h_plank*cc/kb (constant, can be manually overridden)'
            WRITE(stdout,*) '  > Saha_norm - (2pi m_e k_B/h^2)^(3/2) (constant, can be manually overridden)'
            WRITE(stdout,*) '  > mH - mass of Hydrogen (constant, can be manually overridden)'
            WRITE(stdout,*) ''
            WRITE(stdout,*) '  - Available PROCEDURE'
            WRITE(stdout,*) '  > Initialise(ATOMIC_DATA,X_frac,Y_frac,verbose) - Initialise the occupation numbers '
            WRITE(stdout,*) '    > ATOMIC_DATA - input atomic data'
            WRITE(stdout,*) '    > X_frac - (OPTIONAL) input hydrogen mass fraction'
            WRITE(stdout,*) '    > Y_frac - (OPTIONAL) input helium mass fraction'
            WRITE(stdout,*) '    > verbose - (OPTIONAL) enable/desable output '
            WRITE(stdout,*) '  > Set(rho,T,N_electron,verbose) - '
            WRITE(stdout,*) '      Set temperature and density to initiate the computation of occupation numbers'
            WRITE(stdout,*) '    > rho - input density '
            WRITE(stdout,*) '    > T - input temperature '
            WRITE(stdout,*) '    > Te_to_T - (OPTIONAL) input gas-to-radiation temperaturs ratio'
            WRITE(stdout,*) '    > Dilution - (OPTIONAL) input dilution factro'
            WRITE(stdout,*) '    > N_electron - (OPTIONAL) input electrone number Density'
            WRITE(stdout,*) '       If no value provided N_electron is computed internaly'
            WRITE(stdout,*) '    > verbose - (OPTIONAL) enable/desable output '
            WRITE(stdout,*) '  > Info() - give Information about the Type/object'

            IF (Initialised) THEN
                WRITE(stdout,*) '  - Allocated data'
                WRITE(stdout,*) 'This PROCEDURE is under constraction'
            END IF
        END SUBROUTINE Print_NLTE_NUMBER_TYPE_info
        !=======================================================================
        
        SUBROUTINE LoGo()
    
            CALL EXECUTE_COMMAND_LINE('clear')
            WRITE(stdout,*) " "
            WRITE(stdout,*) "     ___           ___           ___           ___           ___           ___       "
            WRITE(stdout,*) "    /\__\         /\  \         /\  \         /\  \         /\  \         /\  \      "
            WRITE(stdout,*) "   /::|  |       /::\  \       /::\  \       /::\  \       /::\  \       /::\  \     "
            WRITE(stdout,*) "  /:|:|  |      /:/\:\  \     /:/\:\  \     /:/\:\  \     /:/\:\  \     /:/\:\  \    "
            WRITE(stdout,*) " /:/|:|__|__   /::\~\:\  \   /:/  \:\  \   /::\~\:\  \   /:/  \:\  \   /::\~\:\  \   "
            WRITE(stdout,*) "/:/ |::::\__\ /:/\:\ \:\__\ /:/__/ \:\__\ /:/\:\ \:\__\ /:/__/ \:\__\ /:/\:\ \:\__\  "
            WRITE(stdout,*) "\/__/~~/:/  / \/__\:\ \/__/ \:\  \ /:/  / \/_|::\/:/  / \:\  \  \/__/ \:\~\:\ \/__/  "
            WRITE(stdout,*) "      /:/  /       \:\__\    \:\  /:/  /     |:|::/  /   \:\  \        \:\ \:\__\    "
            WRITE(stdout,*) "     /:/  /         \/__/     \:\/:/  /      |:|\/__/     \:\  \        \:\ \/__/    "
            WRITE(stdout,*) "    /:/  /                     \::/  /       |:|  |        \:\__\        \:\__\      "
            WRITE(stdout,*) "    \/__/                       \/__/         \|__|         \/__/         \/__/      "
            WRITE(stdout,*) " Version: N.1.5 2023"
            WRITE(stdout,*) " "
            CALL  flush(stdout) 
        END SUBROUTINE LoGo

END MODULE mod_nlte
