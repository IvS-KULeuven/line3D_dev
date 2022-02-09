PROGRAM main
   USE f90getopt
   USE LTE_Line_module
   USE CGS_constants
   IMPLICIT NONE


   ! setup Fortran IO variables 
   INTEGER, PARAMETER :: stdin  = 5
   INTEGER, PARAMETER :: stdout = 6
   INTEGER, PARAMETER :: stderr = 0

   ! seting up the computationsl variables 
   REAL(DP), DIMENSION(:), ALLOCATABLE :: W_i
   REAL(DP), DIMENSION(:), ALLOCATABLE :: q_i
   REAL(DP), DIMENSION(:), ALLOCATABLE :: nu_i0

   REAL(DP), DIMENSION(:), ALLOCATABLE :: tt_list
   REAL(DP), DIMENSION(:), ALLOCATABLE :: lgT_list
   REAL(DP), DIMENSION(:), ALLOCATABLE :: lgD_list

   REAL(DP), DIMENSION(:), ALLOCATABLE :: kape_list
   REAL(DP), DIMENSION(:), ALLOCATABLE :: barQ_list

   CHARACTER(20) PAR_File
   REAL(DP) :: lgTmin,lgTmax,lgDmin, lgDmax, lgttmin, lgttmax
   INTEGER(I4B) :: N_tt, N_lgT, N_lgD
   CHARACTER(20) DIR
   LOGICAL :: ver

   INTEGER(I4B) :: T_gird_counter,D_gird_counter

   REAL(DP) :: D, T, Lam, gl, nl,  chi_i
   REAL(DP) :: tau_s, Mt
   INTEGER(I4B) :: Ll, I, Z
   
   ! INTEGER(I4B) :: exclude = 23
   INTEGER(I4B) :: Nlines(30,9) = 29

   INTEGER(I4B) :: ind1,ind2
   CHARACTER(10) str_2
   CHARACTER(10) str_1
   

   !Define types

   TYPE(ATOMIC_DATA_TYPE) :: ATOM_DATA
   TYPE(LINE_DATA_TYPE)   :: LINE_DATA
   TYPE(OCCUPATION_NUMBER_TYPE):: NUMB_DENS !(30)

   ! Input argument type
   TYPE(option_s):: opts(2)

   ! ... namelist
   NAMELIST / init_param / lgTmin, lgTmax, N_lgT, &
                           lgDmin, lgDmax, N_lgD, &
                           lgttmin, lgttmax, N_tt, &
                           ver, DIR
   !

   CALL EXECUTE_COMMAND_LINE('clear')
   WRITE(*,*) "   ___           ___  ,____          ,__      ,_________   ,_______  "
   WRITE(*,*) "   \  \    ^    /  /  |  _ \         | |      |___  ,___|  |  _____| "
   WRITE(*,*) "    \  \  / \  /  /   | | \ \        | |          | |      | |       "
   WRITE(*,*) "     \  \/ ^ \/  /    | |  | |  ==== | |          | |      |  ---,   "
   WRITE(*,*) "      \   / \   /     | |_/ /        | |          | |      |  ---'   "
   WRITE(*,*) "       \ /   \ /      |  _ _\        | |____      | |      | |_____  "
   WRITE(*,*) "        v     v       |_|  \_\       |,_____|     |_|      |_______| "
   WRITE(*,*) " "
   CALL  flush(stdout) 

   ! Parse the input argument
   opts(1) = option_s( "input", .TRUE., 'i' )

   ! If no options were committed
   IF (command_argument_count() .eq. 0 ) THEN
      STOP 'Please use option --input or -i Input_Parameter_file_name !'
   END IF

   ! Process options one by one
   DO 
      SELECT CASE( getopt( "i:", opts ) ) ! opts is optional (for longopts only)
      CASE( char(0) )
         EXIT
      CASE( 'i' )
         read (optarg, '(A)') PAR_File
         WRITE(*,*) 'Using input parameter file: ',TRIM(PAR_File)
         WRITE(*,*) " "
         CALL  flush(stdout) 
      END SELECT
   END DO


   OPEN(11, FILE=TRIM(PAR_File), STATUS='OLD', FORM='FORMATTED')
   READ(11, NML=init_param)
   CLOSE(11)

   WRITE(*,*) "----- Input Setup -----"
   WRITE(*,*)'         Min                       Max                       N'
   WRITE(*,*)'Lg T  ',lgTmin, lgTmax, N_lgT
   WRITE(*,*)'Lg D  ',lgDmin, lgDmax, N_lgD
   WRITE(*,*)'Lg t  ',lgttmin, lgttmax, N_tt
   WRITE(*,*)'Verbose = ',ver
   WRITE(*,*)'Output destionation -',DIR
   WRITE(*,*) " "
   CALL  flush(stdout) 

   WRITE(*,*)'Initialise Grid'
   CALL  flush(stdout) 
   ALLOCATE(lgT_list(N_lgT))
   CALL linspace(lgTmin, lgTmax,lgT_list)
   ALLOCATE(barQ_list(N_lgT))
   ALLOCATE(kape_list(N_lgT))
   ALLOCATE(lgD_list(N_lgT))
   CALL linspace(lgDmin, lgDmax,lgD_list)
   ALLOCATE(tt_list(N_tt))
   CALL linspace(lgttmin, lgttmax,tt_list)
   IF(ANY(ISNAN(tt_list))) STOP '---- tt is nan ----'
   IF(ANY(ISNAN(lgT_list))) STOP '---- tt is nan ----'
   IF(ANY(ISNAN(lgD_list))) STOP '---- tt is nan ----'
   WRITE(*,*)'Grid - dine'
   WRITE(*,*) " "
   CALL  flush(stdout) 
  
   
   WRITE(*,*)'Creating output destionation - ./',TRIM(DIR)
   CALL EXECUTE_COMMAND_LINE('mkdir -p '//TRIM(DIR), WAIT=.TRUE.)

   ! create Q outpute file
   OPEN (unit = 12, file = './'//TRIM(DIR)//'/Q_TD', FORM='formatted',STATUS='unknown', ACTION='write')
   WRITE(12,*) 0.0d0, lgT_list

   ! create Ke outpute file
   OPEN (unit = 13, file = './'//TRIM(DIR)//'/K_TD', FORM='formatted',STATUS='unknown', ACTION='write')
   WRITE(13,*) 0.0d0, lgT_list
   WRITE(*,*)'Creating output destionation - done'
   WRITE(*,*) " "
   CALL  flush(stdout) 

   ! Step 1) Get atomic data [done], aboundence data [done], Line data [done] (use X = 1 for H no Y|=0 !!)
   WRITE(*,*)'Initialise ATOM and NUMB'
   CALL ATOM_DATA%Initialise(verbose  = ver)
   CALL NUMB_DENS%Initialise(ATOM_DATA,verbose  = ver)
   WRITE(*,*)'ATOM and NUMB  - done'
   WRITE(*,*) " "
   CALL  flush(stdout) 

   ! Get Line data
   WRITE(*,*)'Initialise Line Data'
   CALL LINE_DATA%Initialise(verbose  = ver)

   WRITE(*,*)'  > Number of lines =',LINE_DATA%Total_line_numb
   ALLOCATE(q_i(LINE_DATA%Total_line_numb))
   ALLOCATE(nu_i0(LINE_DATA%Total_line_numb))
   ALLOCATE(W_i(LINE_DATA%Total_line_numb))
   nu_i0(:) = lght_speed /(LINE_DATA%Lambda(:) * Aa2cgs)
   WRITE(*,*)'Line Data - dine'
   WRITE(*,*) " "
   CALL  flush(stdout) 


   DO D_gird_counter = 1,N_lgD
      D = 10.0d0**lgD_list(D_gird_counter)

      DO T_gird_counter = 1,N_lgT
         T = 10.0d0**lgT_list(T_gird_counter)

         WRITE(*,*)'Set lgT =',lgT_list(T_gird_counter),' lgD =',lgD_list(D_gird_counter)
         CALL  flush(stdout)

         CALL Ilumination_finction(T,nu_i0,W_i)
         IF(ANY(ISNAN(W_i))) STOP '---- Wi is nan ----'
         CALL NUMB_DENS%Set(rho = D, T = T, verbose = ver)
         kape_list(T_gird_counter) = NUMB_DENS%Electron * sigma_Thom /D
         WRITE(*,*)'Kappa_e =', kape_list(T_gird_counter)

         WRITE(*,*)'Line strength'
         barQ_list(T_gird_counter) = 0.0d0
         DO ind1 = 1, LINE_DATA%Total_line_numb

            Ll = LINE_DATA%ID(ind1,Ll_)
            I = LINE_DATA%ID(ind1,I_)
            Z = LINE_DATA%ID(ind1,Z_)
            Lam = LINE_DATA%Lambda(ind1)* Aa2cgs

            ! IF(Z.EQ.exclude) CYCLE
            Nlines(Z,I) = Nlines(Z,I)+1

            IF (Ll.GT.ATOM_DATA%List_L(I,Z)) STOP '---- exceeding the available level ----'
            gl = ATOM_DATA%Degeneracy(Ll,I,Z)

            nl = NUMB_DENS%Occupation(Ll,I,Z)
            T =  NUMB_DENS%T

            chi_i = sigma_clas * nl * LINE_DATA%gf_val(ind1) / gl  &
               *(1.0d0 - EXP(-NUMB_DENS%Bolz_norm/(Lam * T)))

            q_i(ind1) = chi_i/(sigma_Thom * NUMB_DENS%Electron) * Lam/lght_speed

            IF(ISNAN(q_i(ind1)).OR.(q_i(ind1).LT.0.0d0)) THEN
               WRITE(*,*) chi_i,nl,T,gl,Lam
               STOP '---- qi is nan ----'
            ENDIF

            barQ_list(T_gird_counter) = barQ_list(T_gird_counter) + q_i(ind1)*W_i(ind1)
         END DO
         WRITE(*,*)'  > \barQ =',barQ_list(T_gird_counter)

         WRITE(*,*)'Compute Mt'

         ! Generate the  output name
         WRITE(str_1,'(F4.2)') LOG10(T)
         WRITE(str_2,'(F5.1)') LOG10(D)

         !create MT_ outpute file
         OPEN (unit = 11, file = './'//TRIM(DIR)//'/Mt_'//TRIM(str_1)//'_'//TRIM(str_2),& 
            FORM='formatted',STATUS='unknown', ACTION='write')

         ! fore each tt
         DO ind1 = 1,N_tt
            ! clean Mt an <tau>
            Mt = 0.0d0

            ! For each line
            DO ind2 = 1, LINE_DATA%Total_line_numb
               tau_s = q_i(ind2) * 10**tt_list(ind1)

               IF(tau_s .EQ. 0.0d0) THEN
                  Mt = Mt + q_i(ind2)*W_i(ind2)
               ELSE
                  Mt = Mt + q_i(ind2)*W_i(ind2) * (1.0d0 - EXP(-tau_s))/tau_s
               END IF
            END DO
            IF(ISNAN(Mt)) STOP '---- Mt is nan ----'
            ! write the CAK t, ansamble average <tau>, and M(t)
            WRITE(11,*) tt_list(ind1), Mt
            ! PRINT*,tt_list(ind1), Mt 
         END DO

         CLOSE(11)

         WRITE(*,*)'WRITE Mt - done'
         WRITE(*,*) " "
         CALL  flush(stdout) 

      END DO !T_gird_counter
      ! write dencity and all temperatures point of bar Q
      WRITE(12,*)lgD_list(D_gird_counter), barQ_list 

      ! write dencity and all temperatures point of kappa electrone
      WRITE(13,*)lgD_list(D_gird_counter), kape_list 
   END DO !D_gird_counter

   CLOSE(12)
   WRITE(*,*)'Write bar Q - done'

   CLOSE(13)
   WRITE(*,*)'Write kappa e - done'
   CALL  flush(stdout) 

CONTAINS

   SUBROUTINE Ilumination_finction(T_io,nu_i0_io,W_i_io)
      REAL(DP), INTENT(IN)  :: T_io
      REAL(DP), INTENT(IN)  :: nu_i0_io(:)
      REAL(DP), INTENT(OUT) :: W_i_io(:)

      REAL(DP), PARAMETER :: C1 = 2.0d0*pi*plnk_const/lght_speed**2.0d0
      REAL(DP), PARAMETER :: C2 = plnk_const/bolz_const
      REAL(DP) :: F

      IF(ANY( nu_i0_io.EQ.0)) STOP '---- Nu = 0 ----'

      F = sigma_stef*T**4.0d0

      W_i_io(:) = C1 * nu_i0_io(:)**4.0d0 / ( EXP( C2*nu_i0_io(:)/T_io ) - 1.0d0) / F
      ! print*,C1, C2, F
   END SUBROUTINE Ilumination_finction


   ! Generates evenly spaced numbers from `from` to `to` (inclusive).
   !
   ! Inputs:
   ! -------
   !
   ! from, to : the lower and upper boundaries of the numbers to generate
   !
  ! Outputs:
   ! -------
   !
   ! array : Array of evenly spaced numbers
   !
   SUBROUTINE linspace(from, to, array)
      REAL(dp), INTENT(in) :: from, to
      REAL(dp), INTENT(out) :: array(:)
      REAL(dp) :: range
      INTEGER :: n, ind
      n = SIZE(array)
      range = to - from

      IF (n == 0) RETURN

      IF (n == 1) THEN
         array(1) = from
         RETURN
      END IF


      DO ind=1, n
         array(ind) = from + range * (ind - 1) / (n - 1)
      END DO
   END SUBROUTINE linspace

END PROGRAM main
