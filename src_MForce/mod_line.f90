MODULE mod_line
    USE prog_type 
    IMPLICIT NONE
    ! Module for read put of theline data files from FASTWIND.
    !
    ! The Module was developed by Luka Poniatowski in 2021
    !   and contains routines ratially based on LTE_Line_module
    !
    ! This module contains the definitions of variables and variable types
    ! used to allocate the atomic, line data, partition functions, occupation numbers, etc.
    !
    ! version 1
    !
    ! Identifiers used in the description for definitions of dimensions and variable types:
    ! UI - Unsingen INTEGER
    ! DP - DOUBLE PRECISION
    ! 1D/2D/3D - 1/2/3 dimentional array
    ! LT - LINE_DATA_TYPE
    ! ----
    !
    ! Within this module variable names are reserd as follows:
    ! Z - Nuclear charge of the atom :: UI
    ! I - Ionisation state   :: UI
    ! L - Energy level  :: UI
    ! Ll / Lu - Lower / Upper Energy level :: UI
    ! LINE_DATA - structure containing line transition data :: LT
    ! ---

    ! Fixed parameters for accesing the line ID Z,I,Ll,Lu values
    INTEGER(I2B), PARAMETER :: Z_ = 1, I_ = 2, Ll_ = 3, Lu_ = 4

    TYPE LINE_DATA_TYPE
            ! This variable type is used to allocate the line transition data

            ! The variables are as follows:

            ! Total number of transitions avvailable
            INTEGER(I4B) :: Total_line_numb

            ! List of line ftransition Identifiers
            ! List must be 2D array of dimention (Total_line_numb,4)
            ! where:
            ! (Total_line_numb,1) - Nuclear charge Z of the element
            ! (Total_line_numb,2) - Ionisation stage I of the element
            ! (Total_line_numb,3) - Lower energy level
            ! (Total_line_numb,4) - Upper energy level
            INTEGER, DIMENSION(:, :), ALLOCATABLE :: ID

            ! Lists of transition wavelength (Angstrom)
            !  and oscilatior strength of specific transition with identifier ID
            ! Lists must be 1D array of dimention (Total_line_numb)
            REAL(SP), DIMENSION(:), ALLOCATABLE :: Lambda
            REAL(SP), DIMENSION(:), ALLOCATABLE :: gf_val

            ! THE Index variable stores the indexes founr by Find routine
            ! It is allocated every time it is called
            INTEGER, DIMENSION(:), ALLOCATABLE :: Index

        CONTAINS
            PROCEDURE, PASS :: Initialise => Read_Line_Data_ascii
            PROCEDURE, PASS :: Find => Find_line_index
            PROCEDURE, PASS :: Info => Print_LINE_DATA_TYPE_info
    END TYPE LINE_DATA_TYPE

    PRIVATE :: Read_Line_Data_ascii, Find_line_index, Print_LINE_DATA_TYPE_info

    !========================================================================
    !========================================================================
    !========================================================================
    !========================================================================

    CONTAINS

        SUBROUTINE Read_Line_Data_ascii(OBJ, verbose)
            ! In-Out and optional variables
            CLASS(LINE_DATA_TYPE), INTENT(INOUT) :: OBJ
            LOGICAL, OPTIONAL :: verbose
            LOGICAL :: ver

            ! Temproraly variables
            INTEGER(I4B) :: Total_line_numb = 0, IO_status = 0, ind = 0!, ind1 = 0, ind2 = 0, ind3 = 0
            INTEGER(I4B), DIMENSION(:), ALLOCATABLE ::  Line_ID
            REAL(SP), DIMENSION(:), ALLOCATABLE :: Lambda
            REAL(SP), DIMENSION(:), ALLOCATABLE :: gf_val

            ! Dummy variables
            INTEGER(I4B) :: idum = 0

            ! String containing the dectination to the project directiory
            CHARACTER(256) DATA_DIR

            ! Pars the input optianal variables
            ! If verbose = .true. specified give datailed status
            ver = .FALSE.
            IF (PRESENT(verbose)) ver = verbose

            IF (ver) WRITE(stdout,*) ' '
            IF (ver) WRITE(stdout,*) 'Initialising the Line Data Type:'

            ! Get the destination ot the data location (check if the enviroment was set)
            CALL get_environment_variable("MFORCE_DIR", DATA_DIR)
            IF (LEN_TRIM(DATA_DIR) .EQ. 0) STOP '---- Error: MFORCE_DIR not set -----'
            DATA_DIR = TRIM(DATA_DIR)//'/DATA'

            ! Chech if atomic data is already Initialised
            IF (ALLOCATED(OBJ%ID)) THEN
                IF (ver) WRITE(stdout,*) '  - Already Initialised, Nothing to do here'
                RETURN
            END IF

            ! Open the line list info file
            OPEN (UNIT=100, FILE=TRIM(DATA_DIR)//'/nl3info', STATUS='OLD', IOSTAT=IO_status)
            ! Check if opening was sucssesful
            IF (IO_status .NE. 0) THEN
                WRITE(stdout,*) '---- Error Opening the '//TRIM(DATA_DIR)//'/nl3info file STATUS ', &
                    IO_status, ' ----'
                STOP
            END IF

            ! Read the total number of lines
            READ (100, *) OBJ%Total_line_numb, idum, idum
            Total_line_numb = OBJ%Total_line_numb
            CLOSE (100)

            IF (ver) WRITE(stdout,*) '  - General info read'
            IF (ver) WRITE(stdout,*) '    > Total number of available transition', OBJ%Total_line_numb

            ! Allocating the variables
            ALLOCATE (Line_ID(Total_line_numb))
            Line_ID = 0
            ALLOCATE (Lambda(Total_line_numb))
            Lambda = 0.0d0
            ALLOCATE (gf_val(Total_line_numb))
            gf_val = 0.0d0
            ALLOCATE (OBJ%ID(Total_line_numb, 4))
            OBJ%ID = 0
            ALLOCATE (OBJ%Lambda(Total_line_numb))
            OBJ%Lambda = 0.0d0
            ALLOCATE (OBJ%gf_val(Total_line_numb))
            OBJ%gf_val = 0.0d0

            ! OPEN line identifiers
            OPEN (UNIT=100, FILE=TRIM(DATA_DIR)//'/nl3i_all_ascii', FORM='formatted', STATUS='OLD', IOSTAT=IO_status)
            ! Check if opening was sucssesful
            IF (IO_status .NE. 0) THEN
                WRITE(stdout,*) '---- Error Opening the '//TRIM(DATA_DIR)//'/nl3i_all_ascii file STATUS ', &
                    IO_status, ' ----'
                STOP
            END IF

            ! OPEN Transition wavelength
            OPEN (UNIT=101, FILE=TRIM(DATA_DIR)//'/nl3a_all_ascii', FORM='formatted', STATUS='OLD', IOSTAT=IO_status)
            ! Check if opening was sucssesful
            IF (IO_status .NE. 0) THEN
                WRITE(stdout,*) '---- Error Opening the '//TRIM(DATA_DIR)//'/nl3a_all_ascii file STATUS ', &
                    IO_status, ' ----'
                STOP
            END IF

            ! OPEN Transition wavelength
            OPEN (UNIT=102, FILE=TRIM(DATA_DIR)//'/nl3g_all_ascii', FORM='formatted', STATUS='OLD', IOSTAT=IO_status)
            ! Check if opening was sucssesful
            IF (IO_status .NE. 0) THEN
                WRITE(stdout,*) '---- Error Opening the '//TRIM(DATA_DIR)//'/nl3g_all_ascii file STATUS ', &
                    IO_status, ' ----'
                STOP
            END IF

            READ (100, *) (Line_ID(ind), ind=1, Total_line_numb)
            READ (101, *) (Lambda(ind), ind=1, Total_line_numb)
            READ (102, *) (gf_val(ind), ind=1, Total_line_numb)
            CLOSE (100)
            CLOSE (101)
            CLOSE (102)

            ! PARSE The Line ID of following format:
            ! zZI0LlLu  (z - is not fixed and can appear or disapperar
            ! depending on length,  0-has a fixd location and is not used)
            ! zZ - Nuclear charge Z
            ! I  - ionisation
            ! 0  - Not used
            ! Ll - Lower level
            ! Lu - Upper level
            ! Examples:
            !  1100102: Z = 1(Hydrogen), I = 1 (Neutral), 0, 01 = ground level, 02 = First excited
            ! 19500317: Z = 19(K), I = 5(4x Ionised), 0, 03 = 3rd excited, 17 = 17th excited
            ! indexis for storing  Z_ = 1, I_ = 2, Ll_ = 3, Lu = 4

            !!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!
            !Total_line_numb = 1
            !!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!
            DO ind = 1, Total_line_numb
                ! bunch of integer arithmetics to parts the IDs
                OBJ%ID(ind, 1) = Line_ID(ind)/1000000
                OBJ%ID(ind, 2) = (Line_ID(ind) - OBJ%ID(ind, 1)*1000000)/100000
                OBJ%ID(ind, 3) = MOD(Line_ID(ind), 10000)/100
                OBJ%ID(ind, 4) = MOD(Line_ID(ind), 100)

                ! check each translation
                IF ((OBJ%ID(ind, 1)*1000000 + OBJ%ID(ind, 2)*100000 + OBJ%ID(ind, 3)*100 &
                    + OBJ%ID(ind, 4)) .NE. Line_ID(ind)) &
                    STOP '---- InconsistencyIn reading the line ID ----'
            END DO
            IF (ver) WRITE(stdout,*) '  - Line IDs read'

            ! Thech that the last entry of the file is correctly read
            !  by comparing it to manually enterd value (Mast be changed id the list is changed)
            IF ((Lambda(Total_line_numb)/3.7174720d+07 - 1.0d0) .GT. 1.0d-8) THEN
                WRITE(stdout,*) '---- Inconsistency in the last wavelength of the list ----'
                WRITE(stdout,*) '  > Expected value:', 3.7174720d+07
                WRITE(stdout,*) '  > Read value:    ', Lambda(Total_line_numb)
                STOP
            END IF

            ! If read was sucssesful store the wavelength array
            OBJ%Lambda(:) = Lambda(:)
            IF (ver) WRITE(stdout,*) '  - Line wavelength read'

            ! Thech that the last entry of the file is correctly read
            !  by comparing it to manually enterd value (Mast be changed id the list is changed)
            IF ((gf_val(Total_line_numb)/3.4000002d-02 - 1.0d0) .GT. 1.0d-8) THEN
                WRITE(stdout,*) '---- Inconsistency in the last wavelength of the list ----'
                WRITE(stdout,*) '  > Expected value:', 3.4000002d-02
                WRITE(stdout,*) '  > Read value:    ', gf_val(Total_line_numb)
                STOP
            END IF

            ! If read was sucssesful store the gf values array
            OBJ%gf_val(:) = gf_val(:)
            IF (ver) WRITE(stdout,*) '  - Line oscilator strength read'

            IF (ver) WRITE(stdout,*) '  - Line Data Type READY'
        END SUBROUTINE Read_Line_Data_ascii
        !========================================================================

        SUBROUTINE Find_line_index(OBJ, Z, I, Ll, Lu, Lambda, Toler, verbose)
            ! In-Out and optional variables
            CLASS(LINE_DATA_TYPE), INTENT(INOUT) ::  OBJ
            ! INTEGER(I4B), INTENT(OUT) :: INDEX
            INTEGER(I4B), OPTIONAL :: Z
            INTEGER(I4B), OPTIONAL :: I
            INTEGER(I4B), OPTIONAL :: Ll
            INTEGER(I4B), OPTIONAL :: Lu
            REAL(SP), OPTIONAL :: Lambda
            REAL(SP), OPTIONAL :: Toler
            LOGICAL, OPTIONAL :: verbose
            LOGICAL :: ver

            INTEGER(I4B) :: ind, ind1, Total_line_numb, ID
            INTEGER, DIMENSION(:), ALLOCATABLE :: ID_Index

            ! Pars the input optianal variables
            ! If verbose = .true. specified give datailed status
            ver = .FALSE.
            IF (PRESENT(verbose)) ver = verbose

            IF (ver) WRITE(stdout,*) ' '
            IF (ver) WRITE(stdout,*) 'Find index CALLED'

            IF (.NOT. ALLOCATED(OBJ%Lambda)) STOP '---- Terminating: Line Data Not Initialized ----'

            ! Deallocate the OBJ%Index array for allocation of new matches
            IF (ALLOCATED(OBJ%Index)) DEALLOCATE(OBJ%Index)

            ! copy the number of lines for coding simplisity
            Total_line_numb = OBJ%Total_line_numb
            ALLOCATE (ID_Index(Total_line_numb))

            ! Check whot inputs are provided fo lookup
            IF (PRESENT(Lambda) .AND. PRESENT(Toler)) THEN
                ! Find by lambda
                IF (ver) WRITE(stdout,*) '  - Search index by input wavelength'

                ! start counting the number of matches
                ind1 = 0

                ! go trhough each line index
                DO ind = 1, Total_line_numb

                    ! If the
                    IF (ABS(OBJ%Lambda(ind)/Lambda - 1) .LE. Toler) THEN
                        ! increase the counter of found matches by 1 (a)
                        ind1 = ind1 + 1

                        ! if the id match save the index of the line
                        ID_Index(ind1) = ind
                    END IF
                END DO

            ELSE IF (PRESENT(Z) .AND. PRESENT(I) .AND. PRESENT(Ll) .AND. PRESENT(Lu)) THEN
                ! Find by ID
                IF (ver) WRITE(stdout,*) '  - Search index by input transition ID (Z,I,Ll,Lu)'

                ! convert input into single integer
                ID = Z*1000000 + I*100000 + Ll*100 + Lu
                IF (ver) WRITE(stdout,*) '   > Looking for:', ID

                ! start counting the number of matches
                ind1 = 0

                ! go trhough each line index
                DO ind = 1, Total_line_numb

                    ! Chech if the ID on the line matches the input ID
                    IF ((OBJ%ID(ind, 1)*1000000 + OBJ%ID(ind, 2)*100000 + &
                        OBJ%ID(ind, 3)*100 + OBJ%ID(ind, 4)) .EQ. ID) THEN

                        ! increase the counter of found matches by 1 (a)
                        ind1 = ind1 + 1

                        ! if the id match save the index of the line
                        ID_Index(ind1) = ind
                    END IF
                END DO
            ELSE
                ! No input provided - abbort
                IF (ver) WRITE(stdout,*) '  - No imput provided'
                ! allocate and store 0
                ALLOCATE (OBJ%Index(1))
                OBJ%Index = 0
            END IF

            IF (ver) WRITE(stdout,*) '  > Found matches:', ind1

            ! if any match is found allocate corresponding number of OBJ%Index
            ! and store the matchin indexes. If no match  found allocate and store 0
            IF (ind1 .GT. 0) THEN
                ALLOCATE (OBJ%Index(ind1))
                OBJ%Index(1:ind1) = ID_Index(1:ind1)
            ELSE
                ALLOCATE (OBJ%Index(1))
                OBJ%Index = 0
            END IF
        END SUBROUTINE Find_line_index
        !=============================================================================

        SUBROUTINE Print_LINE_DATA_TYPE_info(OBJ)
            CLASS(LINE_DATA_TYPE), INTENT(IN) ::  OBJ
            LOGICAL :: Initialised = .FALSE.

            WRITE(stdout,*) ' '
            WRITE(stdout,*) 'Information on LINE_DATA_TYPE requested'

            IF (.NOT. ALLOCATED(OBJ%ID)) THEN
                WRITE(stdout,*) '  - LINE_DATA_TYPE seems NOT Initialised '
                WRITE(stdout,*) '  > For Initialisation CALL LINE_DATA_TYPE_NAME%Initialise()'
                WRITE(stdout,*) '  > Printing ONLY General information'
            ELSE
                WRITE(stdout,*) '  - LINE_DATA_TYPE IS Initialised '
                WRITE(stdout,*) '  > Printing Full information '
                Initialised = .TRUE.
            END IF

            WRITE(stdout,*) '  - Available variables:'
            WRITE(stdout,*) '  > Total_line_numb - Number of available lines'
            WRITE(stdout,*) '  > ID - Transition identifier (output)'
            WRITE(stdout,*) '  > Lambda - Transition wavelength (output)'
            WRITE(stdout,*) '  > gf_val - Transition gf value (output)'
            WRITE(stdout,*) '  > Index - List of line index found (output) '
            WRITE(stdout,*) '  - Available PROCEDURE'
            WRITE(stdout,*) '  > Initialise(verbose)'
            WRITE(stdout,*) '    > verbose - (OPTIONAL) enable/desable output '
            WRITE(stdout,*) '  > Find(Z,I,Ll,Lu,Lambda,Toler,verbose) - '
            WRITE(stdout,*) '    > Z - (OPTIONAL) nuclear change Z of element of interest'
            WRITE(stdout,*) '    > I - (OPTIONAL) ionisation stage I of lelemnt of interest '
            WRITE(stdout,*) '    > Ll - (OPTIONAL) lower level of transition'
            WRITE(stdout,*) '    > Lu - (OPTIONAL) upper level of transition'
            WRITE(stdout,*) '    > Lambda - (OPTIONAL)'
            WRITE(stdout,*) '    > Toler - (OPTIONAL)'
            WRITE(stdout,*) '    > verbose - (OPTIONAL) enable/desable output '
            WRITE(stdout,*) '  > Info() - give Information about the Type/object'

            IF (Initialised) THEN
                WRITE(stdout,*) '  - Allocated data'
            END IF

            WRITE(stdout,*) 'This PROCEDURE is under constraction'
        END SUBROUTINE Print_LINE_DATA_TYPE_info

END MODULE mod_line
