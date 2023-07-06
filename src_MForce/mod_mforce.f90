MODULE mod_mforce
    use mod_nlte
    use mod_line
    use prog_type 
    use fund_const
    implicit none 


    TYPE(LINE_DATA_TYPE) :: LINES
    TYPE(ATOMIC_DATA_TYPE) :: ATOMS
    TYPE(NLTE_NUMBER_TYPE) :: NLTE

    integer(i4b) :: e_z_nlte, e_i_nlte, e_l_nlte, e_u_nlte, gl_nlte, gu_nlte, gf_nlte, lam_nlte

    private
    public :: Initialise_NLTE, Get_NLTE_line, Get_NLTE_kl
    contains

    SUBROUTINE Initialise_NLTE(X_frac, Y_frac)
        real(DP), optional, intent(in):: X_frac, Y_frac 

        call LINES%Initialise()
        call ATOMS%Initialise()
        ! if no input mass fractions provided use solar mixture
        if (PRESENT(X_frac).and.PRESENT(Y_frac)) then
            call NLTE%Initialise(ATOMS, X_frac=X_frac, Y_frac = Y_frac)
        else
            call NLTE%Initialise(ATOMS)
        end if 

    END SUBROUTINE Initialise_NLTE


    SUBROUTINE Get_NLTE_line(e_z, e_i, e_ll, e_lu, gl, gu, gf, lam)
        integer(I4B), intent(in) ::  e_z, e_i, e_ll, e_lu
        real(dp), intent(out) :: gl, gu, gf, lam

        integer(i4b) :: ind 

        call LINES%Find(e_z, e_i, e_ll, e_lu)

        gl = ATOMS%Degeneracy(e_z,e_i,e_ll)
        gu = ATOMS%Degeneracy(e_z,e_i,e_lu)

        gf  = LINES%gf_val(LINES%Index(1)) 
        lam = LINES%Lambda(LINES%Index(1)) 

        if (SIZE(LINES%Index).gt.1) then
            DO ind=1,SIZE(LINES%Index)
                if (gf.lt.LINES%gf_val(LINES%Index(ind))) then
                    gf  = LINES%gf_val(LINES%Index(ind))
                    lam = LINES%Lambda(LINES%Index(ind))
                end if 
            END DO 
        end if

        ! save variables for future use in the get kl function
        e_z_nlte = e_z 
        e_i_nlte = e_i
        e_l_nlte = e_ll
        e_u_nlte = e_lu

        gl_nlte = gl
        gu_nlte = gu
        gf_nlte = gf
        lam_nlte = lam

    END SUBROUTINE Get_NLTE_line

    REAL FUNCTION Get_NLTE_kl(rho, Tr, Tg, W) 
        real(dp), intent(in) :: rho, Tr, Tg, W 

        real(dp) :: TtoT, nl, nu
        real(dp), parameter :: c1 = pi*cgs_e**2 / cgs_me/cgs_clight

        TtoT = Tg/Tr
        call NLTE%Set(rho = rho , T = Tr, Te_to_T = TtoT, Dilution = W)

        nl = NLTE%Occupation(e_z_nlte,e_i_nlte,e_l_nlte)
        nu = NLTE%Occupation(e_z_nlte,e_i_nlte,e_u_nlte)

        get_nlte_kl = c1 * gf_nlte * (nl/gl_nlte - nu/gu_nlte)

    END FUNCTION Get_NLTE_kl

END MODULE mod_mforce