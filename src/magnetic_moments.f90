module magnetic_moments
    use para, only: dp, zi, mu_B, Echarge, hbar, Bohr_radius, eV2Hartree, band_degeneracy_threshold, &
    Num_wann, Nrpts, irvec, HmnR, Bx, By, Bz
    implicit none

    !> Lande g-factor
    real(dp), parameter :: Lande_g_S = 2d0
    real(dp), parameter :: Lande_g_L = 1d0

contains
    subroutine spin_operators(M_S_oper)
        !> extend the 2x2 pauli matrices to Num_wann x Num_wann, without -hbar/2

        use para, only: Package
        implicit none

        integer :: j, nwann
        complex(dp), intent(out) :: M_S_oper(Num_wann, Num_wann, 3)

        M_S_oper= 0d0
        nwann = Num_wann/2

        !> generate Pauli matrix
        !> set Package = 'VASP6' or 'QE' will automately change the spin order to uudd

        do j=1, nwann
            M_S_oper(j,       nwann+j, 1)=  1d0
            M_S_oper(j+nwann, j,       1)=  1d0
            M_S_oper(j,       nwann+j, 2)=  -zi
            M_S_oper(j+nwann, j,       2)=   zi
            M_S_oper(j,       j,       3)=  1d0
            M_S_oper(j+nwann, j+nwann, 3)= -1d0
        enddo
    end subroutine


    subroutine spin_magnetic_moments(UU, M_S)
        !> output <M_S>, in unit of mu_B

        implicit none

        complex(dp), intent(in)  :: UU(Num_wann, Num_wann)
        complex(dp), intent(out) :: M_S(Num_wann, Num_wann, 3)
        complex(dp), allocatable :: UU_dag(:, :), Amat(:, :)

        allocate( UU_dag(Num_wann, Num_wann), Amat(Num_wann, Num_wann) )
        UU_dag= conjg(transpose(UU))
        
        call spin_operators(M_S)

        call mat_mul(Num_wann, M_S(:,:,1), UU, Amat)
        call mat_mul(Num_wann, UU_dag, Amat, M_S(:,:,1)) 
        call mat_mul(Num_wann, M_S(:,:,2), UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, M_S(:,:,2))
        call mat_mul(Num_wann, M_S(:,:,3), UU, Amat) 
        call mat_mul(Num_wann, UU_dag, Amat, M_S(:,:,3))

        M_S = mu_B * (-0.5d0) * Lande_g_S * M_S

    end subroutine spin_magnetic_moments


    subroutine orbital_magnetic_moments(W, velocities, M_L) 
        !> ref: SciPost Phys. 14, 118 (2023), Eq 24b
        !> output <M_L>, in unit of mu_B

        implicit none

        integer :: l, n, p
        real(dp),    intent(in)  :: W(Num_wann)
        complex(dp), intent(in)  :: velocities(Num_wann, Num_wann,3)
        complex(dp), intent(out) :: M_L(Num_wann, Num_wann, 3)

        real(dp) :: dEpl, dEpn, inv_omega_plpn
        M_L = 0d0

        do l= 1, Num_wann
            do n= 1, Num_wann
                do p= 1, Num_wann
                    dEpl = W(p) - W(l)
                    dEpn = W(p) - W(n)
                    if ((ABS(dEpl) < band_degeneracy_threshold) .or. (ABS(dEpn) < band_degeneracy_threshold)) cycle
                    
                    inv_omega_plpn = (1/dEpl + 1/dEpn)
                    M_L(l,n,1) = M_L(l,n,1) + inv_omega_plpn * (velocities(l,p,2) * velocities(p,n,3) - velocities(l,p,3) * velocities(p,n,2))
                    M_L(l,n,2) = M_L(l,n,2) + inv_omega_plpn * (velocities(l,p,3) * velocities(p,n,1) - velocities(l,p,1) * velocities(p,n,3))
                    ! M_L(l,n,3) = M_L(l,n,3) + inv_omega_plpn * (velocities(l,p,1) * velocities(p,n,2) - velocities(l,p,2) * velocities(p,n,1))
                enddo !p
            enddo !n
        enddo !l
        
        M_L = Lande_g_L * M_L /zi/4 * Echarge / hbar * Bohr_radius**2 ! mu_B is already included here
        return
    end subroutine orbital_magnetic_moments


    subroutine orbital_operators(UU, W, velocities, M_L_oper) 
        !> output M_L, without unit

        implicit none

        complex(dp), intent(in)  :: UU(Num_wann, Num_wann)
        real(dp),    intent(in)  :: W(Num_wann)
        complex(dp), intent(in)  :: velocities(Num_wann, Num_wann,3)
        complex(dp), intent(out) :: M_L_oper(Num_wann, Num_wann, 3)

        integer :: l, n
        complex(dp), allocatable :: M_L(:, :, :), UU_dag(:, :), Amat(:, :)

        allocate( M_L(Num_wann, Num_wann, 3), UU_dag(Num_wann, Num_wann), Amat(Num_wann, Num_wann) )
        call orbital_magnetic_moments(W, velocities, M_L)

        M_L_oper = 0d0
        UU_dag= conjg(transpose(UU))

        do l= 1, Num_wann
            do n= 1, Num_wann
                call ZGEMM('N','N',  &
                Num_wann,Num_wann,1, & ! m,n,k
                1.0,                 & ! ALPHA
                UU(:,l:l),Num_wann,  & ! matA, m
                UU_dag(n:n,:),1,     & ! matB, k
                0.0,                 & ! BETA
                Amat,Num_wann)         ! matC, m

                M_L_oper(:,:,1) = M_L_oper(:,:,1) + M_L(l,n,1) * Amat
                M_L_oper(:,:,2) = M_L_oper(:,:,2) + M_L(l,n,2) * Amat
                M_L_oper(:,:,3) = M_L_oper(:,:,3) + M_L(l,n,3) * Amat 
            enddo !n
        enddo !l
                
        M_L_oper = M_L_oper / mu_B / Lande_g_L
        return
    end subroutine orbital_operators


    subroutine add_zeeman_spin() !> be called in readHmnR.f90, subroutine readNormalHmnR
        implicit none

        integer :: ir
        complex(dp), allocatable :: M_S_oper(:, :, :)
        allocate( M_S_oper(Num_wann, Num_wann, 3) )
        call spin_operators(M_S_oper)

        do ir=1, Nrpts
            if (irvec(1, ir)/=0.or.irvec(2, ir)/=0.or.irvec(3, ir)/=0) cycle
            HmnR(:, :, ir) = HmnR(:, :, ir) + mu_B * (-0.5d0) * Lande_g_S *( M_S_oper(:,:,1)*Bx + M_S_oper(:,:,2)*By + M_S_oper(:,:,3)*Bz )
        enddo ! ir
    end subroutine


    subroutine add_zeeman_orb(UU, W, velocities) !> be called in 
        implicit none

        complex(dp), intent(in)  :: UU(Num_wann, Num_wann)
        real(dp),    intent(in)  :: W(Num_wann)
        complex(dp), intent(in)  :: velocities(Num_wann, Num_wann,3)

        integer :: ir
        complex(dp), allocatable :: M_L_oper(:, :, :)
        allocate( M_L_oper(Num_wann, Num_wann, 3) )

        call orbital_operators(UU, W, velocities, M_L_oper) 

        do ir=1, Nrpts
            if (irvec(1, ir)/=0.or.irvec(2, ir)/=0.or.irvec(3, ir)/=0) cycle
            HmnR(:, :, ir) = HmnR(:, :, ir) + mu_B * Lande_g_L *( M_L_oper(:,:,1)*Bx + M_L_oper(:,:,2)*By + M_L_oper(:,:,3)*Bz )
        enddo ! ir

    end subroutine
end module
