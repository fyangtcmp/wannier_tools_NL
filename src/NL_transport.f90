!> Do not parallel these codes over multiple nodes, cause it does not support the OPENMP.
!> the calculations on nonlinear transport are heavy, so the parallel version of wt.x is needed.

module nonlinear_transport
    use para, only: dp, eV2Hartree, Num_wann, OmegaNum, zi, band_degeneracy_threshold, Eta_Arc
    implicit none

    !> magnetic moments in nonlinear planar Hall, see readinput.f90
    ! logical               :: include_m_spin = .false.
    ! logical               :: include_m_orb  = .true.

    !> temperature, eta = k_B*T = 8.617e-5*T
    !> temperature lists:                 20K       50K       70K      100K      200K      300K
    real(dp)     :: Eta_array_fixed(6) = [0.0017d0, 0.0043d0, 0.006d0, 0.0086d0, 0.0172d0, 0.0259d0]*eV2Hartree
    real(dp)     :: Eta_array      (7) !> in module, must be a constant, so 'Eta_array(Eta_number)' is wrong
    integer      :: Eta_number    = 7
    character*40 :: Eta_name, ahcfilename

    real(dp) :: time_start, time_end

    !> Fermi energy, dim= OmegaNum
    real(dp), allocatable :: energy(:)  
    
    !> loop 
    integer(8) :: ie, ik, ik2, knv3
    integer    :: n, m, l, ieta, ierr, icore, icut

    real(dp) :: k(3), kfine(3)

    !> Differential step size at the k-space, in unit of [Length]^-1
    !> Too small value may lead to large error, our tests show that 1e-5 is better than 1e-6 and 1e-7 
    real(dp) :: dkx = 1d-5
    real(dp) :: dky = 1d-5

    !> local fine k-grids
    integer  :: Nk_fine = 5
    integer  :: knv3_fine, ikfine, ikfinex, ikfiney, ikfinez, Nk_adapt
    real(dp) , allocatable :: k_fine_list(:,:) !> Warning! this array will not not explicit declared

    integer(8), allocatable :: ik_adapt_list(:)
    real(dp), allocatable   :: sk_list(:), sk_list_mpi(:) ! the maximum of abs(sigma) at every kpoint

    integer :: ready = 1d0

contains
    subroutine get_k_fine_list()
        use para, only: Nk1, Nk2, Nk3, K3D_vec1_cube, K3D_vec2_cube, K3D_vec3_cube
        implicit none

        if (Nk3<2) then
            knv3_fine = Nk_fine**2
        else 
            knv3_fine = Nk_fine**3
        endif
        allocate(k_fine_list(knv3_fine,3))
    
        k_fine_list = 0d0
        do ikfine=1, knv3_fine
            if (Nk3<2) then
                ikfinex= (ikfine-1)/(Nk_fine)+1
                ikfiney= (ikfine-1-(ikfinex-1)*Nk_fine)+1
                ikfinez= 1
            else 
                ikfinex= (ikfine-1)/(Nk_fine*Nk_fine)+1
                ikfiney= ((ikfine-1-(ikfinex-1)*Nk_fine*Nk_fine)/Nk_fine)+1
                ikfinez= (ikfine-(ikfiney-1)*Nk_fine- (ikfinex-1)*Nk_fine*Nk_fine)
            endif
            k_fine_list(ikfine,:) = K3D_vec1_cube*(ikfinex-1)/dble(Nk1*Nk_fine)  &
                + K3D_vec2_cube*(ikfiney-1)/dble(Nk2*Nk_fine)  &
                + K3D_vec3_cube*(ikfinez-1)/dble(Nk3*Nk_fine)
        enddo
    end subroutine


    !> fixed threshold(ongoing) or relative threshold(default)
    subroutine get_ik_adapt_list() 
        use wmpi, only: cpuid
        use para, only: stdout
        implicit none
        
        real(dp) :: sk_max
        logical, allocatable   :: sk_mask(:)
        allocate( sk_mask(knv3) )
    
        sk_list = log(1d0 + sk_list)
        sk_max  = maxval(sk_list)
    
        do icut = 1, 100
            sk_mask = ( sk_list>=(sk_max*(1d0-icut/100d0)) )
            if ( (sum(sk_list,mask=sk_mask)) / (sum(sk_list)) > 0.9 ) then
                Nk_adapt = count(sk_mask)
                if (cpuid .eq. 0) then
                    write(stdout, '(" ")')
                    write(stdout, '("max = ", E12.3e3, ",  threshold = ", E12.3e3)') exp(sk_max), exp((sk_max*(1d0-icut/100d0)))
                    write(stdout, '("There are ", i15, "/", i18, "  k-points hit the threshold")') Nk_adapt, knv3
                    write(stdout, '(" ")')
                    write(stdout, '("Start to scan the local fine k-grids")')
                endif
                exit
            endif
        enddo

        allocate( ik_adapt_list(Nk_adapt) )
        ik_adapt_list = 0
        l = 0
        do ik = 1,knv3
            if (sk_mask(ik)) then
                l = l + 1
                ik_adapt_list(l) = ik
            endif
        enddo    
    end subroutine

    
    subroutine sigma_SOAHC_int_single_k(k_in, sigma_xyy_k, sigma_yxx_k)
        use para, only:  OmegaMin, OmegaMax
        implicit none
    
        real(dp), intent(in)  :: k_in(3)
        real(dp), intent(out) :: sigma_xyy_k(OmegaNum, Eta_number)
        real(dp), intent(out) :: sigma_yxx_k(OmegaNum, Eta_number)
    
        ! eigen value of H
        real(dp),    allocatable :: W(:)
        complex(dp), allocatable :: Hamk_bulk(:, :)
        complex(dp), allocatable :: UU(:, :)
    
        !> velocities
        complex(dp), allocatable :: velocities(:,:,:)
        complex(dp), allocatable :: vx(:, :), vy(:, :)

        real(dp), allocatable :: diffFermi(:)
        real(dp) :: G_xy, G_yx, G_xx, G_yy     
    
        allocate( W(Num_wann))
        allocate( Hamk_bulk(Num_wann, Num_wann))
        allocate( UU(Num_wann, Num_wann))
        allocate( velocities(Num_wann, Num_wann, 3))

        allocate( diffFermi (OmegaNum))
        allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann) )
    
        Hamk_bulk= 0d0
        UU= 0d0
    
        ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
        call ham_bulk_latticegauge(k_in, Hamk_bulk)
    
        !> diagonalization by call zheev in lapack
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
    
        call velocity_latticegauge_simple(k_in, UU, velocities)
        vx = velocities(:,:,1)
        vy = velocities(:,:,2)

        sigma_xyy_k        = 0d0
        sigma_yxx_k        = 0d0

        do m= 1, Num_wann
            !> At room temperature, the derivatives of the Fermi-Dirac distribution
            !> at (E-Ef)=0.5 is seven orders smaller than (E-Ef)=0.1.
            !> So we choose the energy truncation as [-0.5eV 0.5eV]*eV2Hartree,
            !> it will not affect the precsion and will accelerate the calculations
            if (W(m)<OmegaMin- 2.d-2 .or. W(m)>OmegaMax+ 2.d-2) cycle

            !> calculate G for each band
            G_xx=0d0; G_xy=0d0; G_yx=0d0; G_yy=0d0

            do n= 1, Num_wann
                if (ABS(W(m)-W(n)) < band_degeneracy_threshold) cycle
                G_xx= G_xx+ 2.d0*real(vx(m, n)*vx(n, m)/((W(m)-W(n))**3))
                G_xy= G_xy+ 2.d0*real(vx(m, n)*vy(n, m)/((W(m)-W(n))**3))
                G_yx= G_yx+ 2.d0*real(vy(m, n)*vx(n, m)/((W(m)-W(n))**3))
                G_yy= G_yy+ 2.d0*real(vy(m, n)*vy(n, m)/((W(m)-W(n))**3))
            enddo ! n

            do ieta=1, Eta_number
                !> the outmost band index is m here
                !> this format is very important! prevent NaN error
                diffFermi = -1d0 / (Exp((W(m) - energy)/Eta_array(ieta))+1d0) / (Exp(-(W(m) - energy)/Eta_array(ieta))+1d0) / Eta_array(ieta)

                sigma_xyy_k(:,ieta)= sigma_xyy_k(:,ieta) + (G_yy*real(vx(m,m))-G_xy*real(vy(m,m)))*diffFermi
                sigma_yxx_k(:,ieta)= sigma_yxx_k(:,ieta) + (G_xx*real(vy(m,m))-G_yx*real(vx(m,m)))*diffFermi
            enddo ! ieta
        enddo ! m
    
    end subroutine sigma_SOAHC_int_single_k


    subroutine sigma_NPHC_int_single_k(k_in, Chi_xyyy_k, Chi_yxxx_k)
        use magnetic_moments
        use para
        implicit none
    
        real(dp), intent(in)  :: k_in(3)
        real(dp), intent(out) :: Chi_xyyy_k(OmegaNum, Eta_number, 2, 2) !> the third index: 1=spin, 2=orbital
        real(dp), intent(out) :: Chi_yxxx_k(OmegaNum, Eta_number, 2, 2)
    
        complex(dp), allocatable :: M_S(:, :, :) !> spin magnetic moments
        complex(dp), allocatable :: M_L(:, :, :) !> orbital magnetic moments
            
        real(dp) :: k_dkx(3)
        real(dp) :: k_dky(3)
        
        real(dp), allocatable :: diffFermi(:)

        ! eigen value of H
        real(dp),    allocatable :: W(:)
        complex(dp), allocatable :: Hamk_bulk(:, :)
        complex(dp), allocatable :: UU(:, :)
    
        real(dp), allocatable :: W_dkx(:)
        real(dp), allocatable :: W_dky(:)
    
        complex(dp), allocatable :: sx(:, :), sy(:, :)
        complex(dp), allocatable :: lx(:, :), ly(:, :)
        complex(dp), allocatable :: vx(:, :), vy(:, :)
        complex(dp), allocatable :: velocities(:,:,:)
    
        complex(dp), allocatable :: vx_dkx(:, :), vy_dkx(:, :)  
        complex(dp), allocatable :: vx_dky(:, :), vy_dky(:, :) 
    
        real(dp) :: Lambda_xyy_S, Lambda_yyy_S, Lambda_yxx_S, Lambda_xxx_S
        real(dp) :: Lambda_xyy_L, Lambda_yyy_L, Lambda_yxx_L, Lambda_xxx_L
        real(dp) :: G_xx, G_xy, G_yx, G_yy, G_yy_dkx, G_xy_dky, G_xx_dky, G_yx_dkx
        real(dp) :: dEnm, dEnm3, dEml, dEnl
        
        allocate( diffFermi (OmegaNum))

        !===========================================================================
        !> original kpoints
        allocate( W (Num_wann))
        allocate( Hamk_bulk (Num_wann, Num_wann))
        allocate( UU (Num_wann, Num_wann))
    
        allocate( velocities(Num_wann, Num_wann, 3))
        allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann))
        allocate( sx(Num_wann, Num_wann), sy(Num_wann, Num_wann))
        allocate( lx(Num_wann, Num_wann), ly(Num_wann, Num_wann))
    
        call ham_bulk_latticegauge(k_in, Hamk_bulk)
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
        call velocity_latticegauge_simple(k_in, UU, velocities)
        vx = velocities(:,:,1)
        vy = velocities(:,:,2)
    
        !===========================================================================
        !> magnetic operators 
        allocate( M_S(Num_wann, Num_wann,3) )
        allocate( M_L(Num_wann, Num_wann,3) )
    
        M_S = 0d0
        M_L = 0d0
        if (include_m_spin) then
            call spin_magnetic_moments(UU, M_S)
            sx = M_S(:,:,1)
            sy = M_S(:,:,2)
        endif
        if (include_m_orb) then
            call orbital_magnetic_moments(W, velocities, M_L)
            lx = M_L(:,:,1)
            ly = M_L(:,:,2)
        endif    
    
        !> k + dk_x <===============================================================
        allocate( W_dkx (Num_wann))   
        allocate( vx_dkx(Num_wann, Num_wann), vy_dkx(Num_wann, Num_wann))
    
        k_dkx = k_in+(/Origin_cell%Rua(1)*dkx , Origin_cell%Rub(1)*dkx , Origin_cell%Ruc(1)*dkx/)/twopi
    
        call ham_bulk_latticegauge(k_dkx, Hamk_bulk)
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W_dkx)
        call velocity_latticegauge_simple(k_dkx, UU, velocities)
        vx_dkx = velocities(:,:,1)
        vy_dkx = velocities(:,:,2)
        !===========================================================================
    
        !> k + dk_y <===============================================================
        allocate( W_dky (Num_wann))
        allocate( vx_dky(Num_wann, Num_wann), vy_dky(Num_wann, Num_wann))
    
        k_dky = k_in+(/Origin_cell%Rua(2)*dky , Origin_cell%Rub(2)*dky , Origin_cell%Ruc(2)*dky/)/twopi
    
        call ham_bulk_latticegauge(k_dky, Hamk_bulk)
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W_dky)
        call velocity_latticegauge_simple(k_dky, UU, velocities)
        vx_dky = velocities(:,:,1)
        vy_dky = velocities(:,:,2)
        !===========================================================================
    
        Chi_xyyy_k = 0d0
        Chi_yxxx_k = 0d0
    
        do n= 1, Num_wann
            if (W(n)<OmegaMin- 2.d-2 .or. W(n)>OmegaMax+ 2.d-2) cycle !> prevent NaN error
            G_xx= 0d0
            G_xy= 0d0
            G_yx= 0d0
            G_yy= 0d0
            G_yy_dkx= 0d0 
            G_xy_dky= 0d0 
            G_xx_dky= 0d0 
            G_yx_dkx= 0d0        
            Lambda_xyy_S = 0d0
            Lambda_yyy_S = 0d0
            Lambda_yxx_S = 0d0
            Lambda_xxx_S = 0d0
            Lambda_xyy_L = 0d0
            Lambda_yyy_L = 0d0
            Lambda_yxx_L = 0d0
            Lambda_xxx_L = 0d0
    
            do m= 1, Num_wann
                dEnm= W(n) - W(m)           
                if (ABS(dEnm) < band_degeneracy_threshold) cycle
    
                dEnm3= dEnm**3
                G_xx= G_xx+ 2.d0*real( vx(n, m)*vx(m, n) )/dEnm3
                G_xy= G_xy+ 2.d0*real( vx(n, m)*vy(m, n) )/dEnm3
                G_yx= G_yx+ 2.d0*real( vy(n, m)*vx(m, n) )/dEnm3
                G_yy= G_yy+ 2.d0*real( vy(n, m)*vy(m, n) )/dEnm3
    
                G_yy_dkx= G_yy_dkx + 2.d0*real( vy_dkx(n, m)*vy_dkx(m, n) )/(W_dkx(n) - W_dkx(m))**3
                G_yx_dkx= G_yx_dkx + 2.d0*real( vy_dkx(n, m)*vx_dkx(m, n) )/(W_dkx(n) - W_dkx(m))**3
                
                G_xy_dky= G_xy_dky + 2.d0*real( vx_dky(n, m)*vy_dky(m, n) )/(W_dky(n) - W_dky(m))**3
                G_xx_dky= G_xx_dky + 2.d0*real( vx_dky(n, m)*vx_dky(m, n) )/(W_dky(n) - W_dky(m))**3
    
                if (include_m_spin) then
                    Lambda_xyy_S = Lambda_xyy_S + 6.d0* real( vx(n, m)*vy(m, n)*(sy(n, n)-sy(m, m)) )/dEnm3/dEnm
                    Lambda_yyy_S = Lambda_yyy_S + 6.d0* real( vy(n, m)*vy(m, n)*(sy(n, n)-sy(m, m)) )/dEnm3/dEnm
                    Lambda_yxx_S = Lambda_yxx_S + 6.d0* real( vy(n, m)*vx(m, n)*(sx(n, n)-sx(m, m)) )/dEnm3/dEnm
                    Lambda_xxx_S = Lambda_xxx_S + 6.d0* real( vx(n, m)*vx(m, n)*(sx(n, n)-sx(m, m)) )/dEnm3/dEnm
                    
                    do l= 1, Num_wann
                        dEnl= W(n)-W(l)
                        dEml= W(m)-W(l)                    
                        if (ABS(dEnl) > band_degeneracy_threshold) then
                            Lambda_xyy_S = Lambda_xyy_S - 2.d0* real( (vx(l, m)*vy(m, n)+vy(l, m)*vx(m, n)) *sy(n, l)) /dEnm3/dEnl
                            Lambda_yyy_S = Lambda_yyy_S - 2.d0* real( (vy(l, m)*vy(m, n)+vy(l, m)*vy(m, n)) *sy(n, l)) /dEnm3/dEnl
                            Lambda_yxx_S = Lambda_yxx_S - 2.d0* real( (vy(l, m)*vx(m, n)+vx(l, m)*vy(m, n)) *sx(n, l)) /dEnm3/dEnl
                            Lambda_xxx_S = Lambda_xxx_S - 2.d0* real( (vx(l, m)*vx(m, n)+vx(l, m)*vx(m, n)) *sx(n, l)) /dEnm3/dEnl
                        endif
                        if (ABS(dEml) > band_degeneracy_threshold) then
                            Lambda_xyy_S = Lambda_xyy_S - 2.d0* real( (vx(l, n)*vy(n, m)+vy(l, n)*vx(n, m)) *sy(m, l)) /dEnm3/dEml
                            Lambda_yyy_S = Lambda_yyy_S - 2.d0* real( (vy(l, n)*vy(n, m)+vy(l, n)*vy(n, m)) *sy(m, l)) /dEnm3/dEml
                            Lambda_yxx_S = Lambda_yxx_S - 2.d0* real( (vy(l, n)*vx(n, m)+vx(l, n)*vy(n, m)) *sx(m, l)) /dEnm3/dEml
                            Lambda_xxx_S = Lambda_xxx_S - 2.d0* real( (vx(l, n)*vx(n, m)+vx(l, n)*vx(n, m)) *sx(m, l)) /dEnm3/dEml
                        endif
                    enddo ! l
                endif

                if (include_m_orb) then
                    Lambda_xyy_L = Lambda_xyy_L + 6.d0* real( vx(n, m)*vy(m, n)*(ly(n, n)-ly(m, m)) )/dEnm3/dEnm
                    Lambda_yyy_L = Lambda_yyy_L + 6.d0* real( vy(n, m)*vy(m, n)*(ly(n, n)-ly(m, m)) )/dEnm3/dEnm
                    Lambda_yxx_L = Lambda_yxx_L + 6.d0* real( vy(n, m)*vx(m, n)*(lx(n, n)-lx(m, m)) )/dEnm3/dEnm
                    Lambda_xxx_L = Lambda_xxx_L + 6.d0* real( vx(n, m)*vx(m, n)*(lx(n, n)-lx(m, m)) )/dEnm3/dEnm
                    
                    do l= 1, Num_wann
                        dEnl= W(n)-W(l)
                        dEml= W(m)-W(l)                    
                        if (ABS(dEnl) > band_degeneracy_threshold) then
                            Lambda_xyy_L = Lambda_xyy_L - 2.d0* real( (vx(l, m)*vy(m, n)+vy(l, m)*vx(m, n)) *ly(n, l)) /dEnm3/dEnl
                            Lambda_yyy_L = Lambda_yyy_L - 2.d0* real( (vy(l, m)*vy(m, n)+vy(l, m)*vy(m, n)) *ly(n, l)) /dEnm3/dEnl
                            Lambda_yxx_L = Lambda_yxx_L - 2.d0* real( (vy(l, m)*vx(m, n)+vx(l, m)*vy(m, n)) *lx(n, l)) /dEnm3/dEnl
                            Lambda_xxx_L = Lambda_xxx_L - 2.d0* real( (vx(l, m)*vx(m, n)+vx(l, m)*vx(m, n)) *lx(n, l)) /dEnm3/dEnl
                        endif
                        if (ABS(dEml) > band_degeneracy_threshold) then
                            Lambda_xyy_L = Lambda_xyy_L - 2.d0* real( (vx(l, n)*vy(n, m)+vy(l, n)*vx(n, m)) *ly(m, l)) /dEnm3/dEml
                            Lambda_yyy_L = Lambda_yyy_L - 2.d0* real( (vy(l, n)*vy(n, m)+vy(l, n)*vy(n, m)) *ly(m, l)) /dEnm3/dEml
                            Lambda_yxx_L = Lambda_yxx_L - 2.d0* real( (vy(l, n)*vx(n, m)+vx(l, n)*vy(n, m)) *lx(m, l)) /dEnm3/dEml
                            Lambda_xxx_L = Lambda_xxx_L - 2.d0* real( (vx(l, n)*vx(n, m)+vx(l, n)*vx(n, m)) *lx(m, l)) /dEnm3/dEml
                        endif
                    enddo ! l
                endif
    
            enddo ! m

            do ieta=1, Eta_number
                !> the outmost band index is n here
                !> this format is very important! prevent NaN error
                diffFermi = -1d0 / (Exp((W(n) - energy)/Eta_array(ieta))+1d0) / (Exp(-(W(n) - energy)/Eta_array(ieta))+1d0) / Eta_array(ieta)

                if (include_m_spin) then
                    Chi_xyyy_k(:,ieta, 1, 1) = Chi_xyyy_k(:,ieta, 1, 1) + real( (vx(n,n)*Lambda_yyy_S - vy(n,n)*Lambda_xyy_S) ) * diffFermi
                    Chi_xyyy_k(:,ieta, 1, 2) = Chi_xyyy_k(:,ieta, 1, 2) + real( ((G_yy_dkx - G_yy)/dkx - (G_xy_dky - G_xy)/dky)*sy(n,n) ) * diffFermi

                    Chi_yxxx_k(:,ieta, 1, 1) = Chi_yxxx_k(:,ieta, 1, 1) + real( (vy(n,n)*Lambda_xxx_S - vx(n,n)*Lambda_yxx_S) ) * diffFermi
                    Chi_yxxx_k(:,ieta, 1, 2) = Chi_yxxx_k(:,ieta, 1, 2) + real( ((G_xx_dky - G_xx)/dky - (G_yx_dkx - G_yx)/dkx)*sx(n,n) ) * diffFermi
                endif
                if (include_m_orb) then
                    Chi_xyyy_k(:,ieta, 2, 1) = Chi_xyyy_k(:,ieta, 2, 1) + real( (vx(n,n)*Lambda_yyy_L - vy(n,n)*Lambda_xyy_L) ) * diffFermi
                    Chi_xyyy_k(:,ieta, 2, 2) = Chi_xyyy_k(:,ieta, 2, 2) + real( ((G_yy_dkx - G_yy)/dkx - (G_xy_dky - G_xy)/dky)*ly(n,n) ) * diffFermi

                    Chi_yxxx_k(:,ieta, 2, 1) = Chi_yxxx_k(:,ieta, 2, 1) + real( (vy(n,n)*Lambda_xxx_L - vx(n,n)*Lambda_yxx_L) ) * diffFermi
                    Chi_yxxx_k(:,ieta, 2, 2) = Chi_yxxx_k(:,ieta, 2, 2) + real( ((G_xx_dky - G_xx)/dky - (G_yx_dkx - G_yx)/dkx)*lx(n,n) ) * diffFermi
                endif
            enddo ! ieta
        enddo ! n
    
    end subroutine


    subroutine sigma_NPHC_tau2_single_k(k_in, Chi_xyyy_k, Chi_yxxx_k)
        use magnetic_moments
        use para
        implicit none
    
        real(dp), intent(in)  :: k_in(3)
        real(dp), intent(out) :: Chi_xyyy_k(OmegaNum, Eta_number, 2, 2) !> the third index: 1=spin, 2=orbital
        real(dp), intent(out) :: Chi_yxxx_k(OmegaNum, Eta_number, 2, 2)
    
        complex(dp), allocatable :: M_S(:, :, :) !> spin magnetic moments
        complex(dp), allocatable :: M_L(:, :, :) !> orbital magnetic moments
            
        real(dp) :: k_dkx(3)
        real(dp) :: k_dky(3)
        
        real(dp), allocatable :: Fshort(:) !> short notation of the denominator of the Fermi distribution
        real(dp), allocatable :: diffFermi(:), diff2Fermi(:)

        ! eigen value of H
        real(dp),    allocatable :: W(:)
        complex(dp), allocatable :: Hamk_bulk(:, :)
        complex(dp), allocatable :: UU(:, :)
    
        real(dp), allocatable :: W_dkx(:)
        real(dp), allocatable :: W_dky(:)
    
        complex(dp), allocatable :: sx(:, :), sy(:, :)
        complex(dp), allocatable :: lx(:, :), ly(:, :)
        complex(dp), allocatable :: vx(:, :), vy(:, :)
        complex(dp), allocatable :: velocities(:,:,:)
    
        complex(dp), allocatable :: vx_dkx(:, :), vy_dkx(:, :)  
        complex(dp), allocatable :: vx_dky(:, :), vy_dky(:, :)
        complex(dp), allocatable :: sy_dkx(:, :), ly_dkx(:, :)  
        complex(dp), allocatable :: sx_dky(:, :), lx_dky(:, :) 
    
        real(dp) :: alpha_xyyy_S, alpha_xyyy_L, alpha_yxxx_S, alpha_yxxx_L
        
        allocate( Fshort(OmegaNum), diffFermi(OmegaNum), diff2Fermi(OmegaNum))

        !===========================================================================
        !> original kpoints
        allocate( W (Num_wann))
        allocate( Hamk_bulk (Num_wann, Num_wann))
        allocate( UU (Num_wann, Num_wann))
    
        allocate( velocities(Num_wann, Num_wann, 3))
        allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann))
        allocate( sx(Num_wann, Num_wann), sy(Num_wann, Num_wann))
        allocate( lx(Num_wann, Num_wann), ly(Num_wann, Num_wann))
    
        call ham_bulk_latticegauge(k_in, Hamk_bulk)
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
        call velocity_latticegauge_simple(k_in, UU, velocities)
        vx = velocities(:,:,1)
        vy = velocities(:,:,2)
    
        !> magnetic mats
        allocate( M_S(Num_wann, Num_wann,3) )
        allocate( M_L(Num_wann, Num_wann,3) )
    
        M_S = 0d0
        M_L = 0d0
        if (include_m_spin) then
            call spin_magnetic_moments(UU, M_S)
            sx = M_S(:,:,1)
            sy = M_S(:,:,2)
        endif
        if (include_m_orb) then
            call orbital_magnetic_moments(W, velocities, M_L)
            lx = M_L(:,:,1)
            ly = M_L(:,:,2)
        endif    
    
        !> k + dk_x <===============================================================
        allocate( W_dkx (Num_wann))   
        allocate( vx_dkx(Num_wann, Num_wann), vy_dkx(Num_wann, Num_wann))
        allocate( sy_dkx(Num_wann, Num_wann), ly_dkx(Num_wann, Num_wann))
    
        k_dkx = k_in+(/Origin_cell%Rua(1)*dkx , Origin_cell%Rub(1)*dkx , Origin_cell%Ruc(1)*dkx/)/twopi
    
        call ham_bulk_latticegauge(k_dkx, Hamk_bulk)
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W_dkx)
        call velocity_latticegauge_simple(k_dkx, UU, velocities)
        vx_dkx = velocities(:,:,1)
        vy_dkx = velocities(:,:,2)

        if (include_m_spin) then
            call spin_magnetic_moments(UU, M_S)
            sy_dkx = M_S(:,:,2)
        endif
        if (include_m_orb) then
            call orbital_magnetic_moments(W_dkx, velocities, M_L)
            ly_dkx = M_L(:,:,2)
        endif  
        !===========================================================================
    
        !> k + dk_y <===============================================================
        allocate( W_dky (Num_wann))
        allocate( vx_dky(Num_wann, Num_wann), vy_dky(Num_wann, Num_wann))
        allocate( sx_dky(Num_wann, Num_wann), lx_dky(Num_wann, Num_wann))
    
        k_dky = k_in+(/Origin_cell%Rua(2)*dky , Origin_cell%Rub(2)*dky , Origin_cell%Ruc(2)*dky/)/twopi
    
        call ham_bulk_latticegauge(k_dky, Hamk_bulk)
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W_dky)
        call velocity_latticegauge_simple(k_dky, UU, velocities)
        vx_dky = velocities(:,:,1)
        vy_dky = velocities(:,:,2)

        if (include_m_spin) then
            call spin_magnetic_moments(UU, M_S)
            sx_dky = M_S(:,:,1)
        endif
        if (include_m_orb) then
            call orbital_magnetic_moments(W_dky, velocities, M_L)
            lx_dky = M_L(:,:,1)
        endif  
        !===========================================================================
    
        Chi_xyyy_k = 0d0
        Chi_yxxx_k = 0d0
    
        do n= 1, Num_wann
            if (W(n)<OmegaMin- 2.d-2 .or. W(n)>OmegaMax+ 2.d-2) cycle !> prevent NaN error

            alpha_xyyy_S= 0d0
            alpha_yxxx_S= 0d0
            alpha_xyyy_L= 0d0
            alpha_yxxx_L= 0d0

            if (include_m_spin) then
                alpha_xyyy_S = real( (sy_dkx(n,n) - sy(n,n))/dkx*(vy(n,n)**2) - (vy_dky(n,n) - vy(n,n))/dky*sy(n,n)*vx(n,n) )
                alpha_yxxx_S = real( (sx_dky(n,n) - sx(n,n))/dky*(vx(n,n)**2) - (vx_dkx(n,n) - vx(n,n))/dkx*sx(n,n)*vy(n,n) )
            endif

            if (include_m_orb) then
                alpha_xyyy_L = real( (ly_dkx(n,n) - ly(n,n))/dkx*(vy(n,n)**2) - (vy_dky(n,n) - vy(n,n))/dky*ly(n,n)*vx(n,n) )
                alpha_yxxx_L = real( (lx_dky(n,n) - lx(n,n))/dky*(vx(n,n)**2) - (vx_dkx(n,n) - vx(n,n))/dkx*lx(n,n)*vy(n,n) )
            endif
    
            do ieta=1, Eta_number
                !> the outmost band index is n here
                Fshort = Exp((W(n)-energy)/Eta_array(ieta))
                !> this format is very important! prevent NaN error
                diffFermi  = -1d0 / (Fshort+1d0) / (1d0/Fshort+1d0) / Eta_array(ieta)
                diff2Fermi = diffFermi * ( 1d0 - 2d0 / (1d0/Fshort+1d0)) / Eta_array(ieta)

                if (include_m_spin) then
                    Chi_xyyy_k(:,ieta, 1, 1) = Chi_xyyy_k(:,ieta, 1, 1) + alpha_xyyy_S * diff2Fermi
                    Chi_yxxx_k(:,ieta, 1, 1) = Chi_yxxx_k(:,ieta, 1, 1) + alpha_yxxx_S * diff2Fermi
                endif
                if (include_m_orb) then
                    Chi_xyyy_k(:,ieta, 2, 1) = Chi_xyyy_k(:,ieta, 2, 1) + alpha_xyyy_L * diff2Fermi
                    Chi_yxxx_k(:,ieta, 2, 1) = Chi_yxxx_k(:,ieta, 2, 1) + alpha_yxxx_L * diff2Fermi
                endif
            enddo ! ieta
        enddo ! n
    
    end subroutine


    subroutine sigma_NPHC_int_single_k_2(k_in, Chi_xyyx_k, Chi_yxxy_k)
        use magnetic_moments
        use para
        implicit none
    
        real(dp), intent(in)  :: k_in(3)
        real(dp), intent(out) :: Chi_xyyx_k(OmegaNum, Eta_number, 2, 2) !> the third index: 1=spin, 2=orbital
        real(dp), intent(out) :: Chi_yxxy_k(OmegaNum, Eta_number, 2, 2)
    
        complex(dp), allocatable :: M_S(:, :, :) !> spin magnetic moments
        complex(dp), allocatable :: M_L(:, :, :) !> orbital magnetic moments
            
        real(dp) :: k_dkx(3)
        real(dp) :: k_dky(3)
        
        real(dp), allocatable :: diffFermi(:)

        ! eigen value of H
        real(dp),    allocatable :: W(:)
        complex(dp), allocatable :: Hamk_bulk(:, :)
        complex(dp), allocatable :: UU(:, :)
    
        real(dp), allocatable :: W_dkx(:)
        real(dp), allocatable :: W_dky(:)
    
        complex(dp), allocatable :: sx(:, :), sy(:, :)
        complex(dp), allocatable :: lx(:, :), ly(:, :)
        complex(dp), allocatable :: vx(:, :), vy(:, :)
        complex(dp), allocatable :: velocities(:,:,:)
    
        complex(dp), allocatable :: vx_dkx(:, :), vy_dkx(:, :)  
        complex(dp), allocatable :: vx_dky(:, :), vy_dky(:, :) 
    
        real(dp) :: Lambda_xxy_S, Lambda_yxy_S, Lambda_yyx_S, Lambda_xyx_S
        real(dp) :: Lambda_xxy_L, Lambda_yxy_L, Lambda_yyx_L, Lambda_xyx_L
        real(dp) :: G_xx, G_xy, G_yx, G_yy, G_yy_dkx, G_xy_dky, G_xx_dky, G_yx_dkx
        real(dp) :: dEnm, dEnm3, dEml, dEnl
        
        allocate( diffFermi (OmegaNum))

        !===========================================================================
        !> original kpoints
        allocate( W (Num_wann))
        allocate( Hamk_bulk (Num_wann, Num_wann))
        allocate( UU (Num_wann, Num_wann))
    
        allocate( velocities(Num_wann, Num_wann, 3))
        allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann))
        allocate( sx(Num_wann, Num_wann), sy(Num_wann, Num_wann))
        allocate( lx(Num_wann, Num_wann), ly(Num_wann, Num_wann))
    
        call ham_bulk_latticegauge(k_in, Hamk_bulk)
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
        call velocity_latticegauge_simple(k_in, UU, velocities)
        vx = velocities(:,:,1)
        vy = velocities(:,:,2)
    
        !===========================================================================
        !> magnetic operators 
        allocate( M_S(Num_wann, Num_wann,3) )
        allocate( M_L(Num_wann, Num_wann,3) )
    
        M_S = 0d0
        M_L = 0d0
        if (include_m_spin) then
            call spin_magnetic_moments(UU, M_S)
            sx = M_S(:,:,1)
            sy = M_S(:,:,2)
        endif
        if (include_m_orb) then
            call orbital_magnetic_moments(W, velocities, M_L)
            lx = M_L(:,:,1)
            ly = M_L(:,:,2)
        endif    
    
        !> k + dk_x <===============================================================
        allocate( W_dkx (Num_wann))   
        allocate( vx_dkx(Num_wann, Num_wann), vy_dkx(Num_wann, Num_wann))
    
        k_dkx = k_in+(/Origin_cell%Rua(1)*dkx , Origin_cell%Rub(1)*dkx , Origin_cell%Ruc(1)*dkx/)/twopi
    
        call ham_bulk_latticegauge(k_dkx, Hamk_bulk)
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W_dkx)
        call velocity_latticegauge_simple(k_dkx, UU, velocities)
        vx_dkx = velocities(:,:,1)
        vy_dkx = velocities(:,:,2)
        !===========================================================================
    
        !> k + dk_y <===============================================================
        allocate( W_dky (Num_wann))
        allocate( vx_dky(Num_wann, Num_wann), vy_dky(Num_wann, Num_wann))
    
        k_dky = k_in+(/Origin_cell%Rua(2)*dky , Origin_cell%Rub(2)*dky , Origin_cell%Ruc(2)*dky/)/twopi
    
        call ham_bulk_latticegauge(k_dky, Hamk_bulk)
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W_dky)
        call velocity_latticegauge_simple(k_dky, UU, velocities)
        vx_dky = velocities(:,:,1)
        vy_dky = velocities(:,:,2)
        !===========================================================================
    
        Chi_xyyx_k = 0d0
        Chi_yxxy_k = 0d0
    
        do n= 1, Num_wann
            if (W(n)<OmegaMin- 2.d-2 .or. W(n)>OmegaMax+ 2.d-2) cycle !> prevent NaN error
            G_xx= 0d0
            G_xy= 0d0
            G_yx= 0d0
            G_yy= 0d0
            G_yy_dkx= 0d0 
            G_xy_dky= 0d0 
            G_xx_dky= 0d0 
            G_yx_dkx= 0d0        
            Lambda_xxy_S = 0d0
            Lambda_yxy_S = 0d0
            Lambda_yyx_S = 0d0
            Lambda_xyx_S = 0d0
            Lambda_xxy_L = 0d0
            Lambda_yxy_L = 0d0
            Lambda_yyx_L = 0d0
            Lambda_xyx_L = 0d0
    
            do m= 1, Num_wann
                dEnm= W(n) - W(m)           
                if (ABS(dEnm) < band_degeneracy_threshold) cycle
    
                dEnm3= dEnm**3
                G_xx= G_xx+ 2.d0*real( vx(n, m)*vx(m, n) )/dEnm3
                G_xy= G_xy+ 2.d0*real( vx(n, m)*vy(m, n) )/dEnm3
                G_yx= G_yx+ 2.d0*real( vy(n, m)*vx(m, n) )/dEnm3
                G_yy= G_yy+ 2.d0*real( vy(n, m)*vy(m, n) )/dEnm3
    
                G_yy_dkx= G_yy_dkx + 2.d0*real( vy_dkx(n, m)*vy_dkx(m, n) )/(W_dkx(n) - W_dkx(m))**3
                G_yx_dkx= G_yx_dkx + 2.d0*real( vy_dkx(n, m)*vx_dkx(m, n) )/(W_dkx(n) - W_dkx(m))**3
                
                G_xy_dky= G_xy_dky + 2.d0*real( vx_dky(n, m)*vy_dky(m, n) )/(W_dky(n) - W_dky(m))**3
                G_xx_dky= G_xx_dky + 2.d0*real( vx_dky(n, m)*vx_dky(m, n) )/(W_dky(n) - W_dky(m))**3
    
                if (include_m_spin) then
                    Lambda_xxy_S = Lambda_xxy_S + 6.d0* real( vx(n, m)*vx(m, n)*(sy(n, n)-sy(m, m)) )/dEnm3/dEnm
                    Lambda_yxy_S = Lambda_yxy_S + 6.d0* real( vy(n, m)*vx(m, n)*(sy(n, n)-sy(m, m)) )/dEnm3/dEnm
                    Lambda_yyx_S = Lambda_yyx_S + 6.d0* real( vy(n, m)*vy(m, n)*(sx(n, n)-sx(m, m)) )/dEnm3/dEnm
                    Lambda_xyx_S = Lambda_xyx_S + 6.d0* real( vx(n, m)*vy(m, n)*(sx(n, n)-sx(m, m)) )/dEnm3/dEnm
                    
                    do l= 1, Num_wann
                        dEnl= W(n)-W(l)
                        dEml= W(m)-W(l)                    
                        if (ABS(dEnl) > band_degeneracy_threshold) then
                            Lambda_xxy_S = Lambda_xxy_S - 2.d0* real( (vx(l, m)*vx(m, n)*2                ) *sy(n, l)) /dEnm3/dEnl
                            Lambda_yxy_S = Lambda_yxy_S - 2.d0* real( (vy(l, m)*vx(m, n)+vx(l, m)*vy(m, n)) *sy(n, l)) /dEnm3/dEnl
                            Lambda_yyx_S = Lambda_yyx_S - 2.d0* real( (vy(l, m)*vy(m, n)*2                ) *sx(n, l)) /dEnm3/dEnl
                            Lambda_xyx_S = Lambda_xyx_S - 2.d0* real( (vx(l, m)*vy(m, n)+vy(l, m)*vx(m, n)) *sx(n, l)) /dEnm3/dEnl
                        endif
                        if (ABS(dEml) > band_degeneracy_threshold) then
                            Lambda_xxy_S = Lambda_xxy_S - 2.d0* real( (vx(l, n)*vx(n, m)*2                ) *sy(m, l)) /dEnm3/dEml
                            Lambda_yxy_S = Lambda_yxy_S - 2.d0* real( (vy(l, n)*vx(n, m)+vx(l, n)*vy(n, m)) *sy(m, l)) /dEnm3/dEml
                            Lambda_yyx_S = Lambda_yyx_S - 2.d0* real( (vy(l, n)*vy(n, m)*2                ) *sx(m, l)) /dEnm3/dEml
                            Lambda_xyx_S = Lambda_xyx_S - 2.d0* real( (vx(l, n)*vy(n, m)+vy(l, n)*vx(n, m)) *sx(m, l)) /dEnm3/dEml
                        endif
                    enddo ! l
                endif

                if (include_m_orb) then
                    Lambda_xxy_L = Lambda_xxy_L + 6.d0* real( vx(n, m)*vx(m, n)*(sy(n, n)-ly(m, m)) )/dEnm3/dEnm
                    Lambda_yxy_L = Lambda_yxy_L + 6.d0* real( vy(n, m)*vx(m, n)*(sy(n, n)-ly(m, m)) )/dEnm3/dEnm
                    Lambda_yyx_L = Lambda_yyx_L + 6.d0* real( vy(n, m)*vy(m, n)*(sx(n, n)-lx(m, m)) )/dEnm3/dEnm
                    Lambda_xyx_L = Lambda_xyx_L + 6.d0* real( vx(n, m)*vy(m, n)*(sx(n, n)-lx(m, m)) )/dEnm3/dEnm
                    
                    do l= 1, Num_wann
                        dEnl= W(n)-W(l)
                        dEml= W(m)-W(l)                    
                        if (ABS(dEnl) > band_degeneracy_threshold) then
                            Lambda_xxy_L = Lambda_xxy_L - 2.d0* real( (vx(l, m)*vx(m, n)*2                ) *ly(n, l)) /dEnm3/dEnl
                            Lambda_yxy_L = Lambda_yxy_L - 2.d0* real( (vy(l, m)*vx(m, n)+vx(l, m)*vy(m, n)) *ly(n, l)) /dEnm3/dEnl
                            Lambda_yyx_L = Lambda_yyx_L - 2.d0* real( (vy(l, m)*vy(m, n)*2                ) *lx(n, l)) /dEnm3/dEnl
                            Lambda_xyx_L = Lambda_xyx_L - 2.d0* real( (vx(l, m)*vy(m, n)+vy(l, m)*vx(m, n)) *lx(n, l)) /dEnm3/dEnl
                        endif
                        if (ABS(dEml) > band_degeneracy_threshold) then
                            Lambda_xxy_L = Lambda_xxy_L - 2.d0* real( (vx(l, n)*vx(n, m)*2                ) *ly(m, l)) /dEnm3/dEml
                            Lambda_yxy_L = Lambda_yxy_L - 2.d0* real( (vy(l, n)*vx(n, m)+vx(l, n)*vy(n, m)) *ly(m, l)) /dEnm3/dEml
                            Lambda_yyx_L = Lambda_yyx_L - 2.d0* real( (vy(l, n)*vy(n, m)*2                ) *lx(m, l)) /dEnm3/dEml
                            Lambda_xyx_L = Lambda_xyx_L - 2.d0* real( (vx(l, n)*vy(n, m)+vy(l, n)*vx(n, m)) *lx(m, l)) /dEnm3/dEml
                        endif
                    enddo ! l
                endif
    
            enddo ! m

            do ieta=1, Eta_number
                !> the outmost band index is n here
                !> this format is very important! prevent NaN error
                diffFermi = -1d0 / (Exp((W(n) - energy)/Eta_array(ieta))+1d0) / (Exp(-(W(n) - energy)/Eta_array(ieta))+1d0) / Eta_array(ieta)

                if (include_m_spin) then
                    Chi_xyyx_k(:,ieta, 1, 1) = Chi_xyyx_k(:,ieta, 1, 1) + real( (vx(n,n)*Lambda_yyx_S - vy(n,n)*Lambda_xyx_S) ) * diffFermi
                    Chi_xyyx_k(:,ieta, 1, 2) = Chi_xyyx_k(:,ieta, 1, 2) + real( ((G_yy_dkx - G_yy)/dkx - (G_xy_dky - G_xy)/dky)*sx(n,n) ) * diffFermi

                    Chi_yxxy_k(:,ieta, 1, 1) = Chi_yxxy_k(:,ieta, 1, 1) + real( (vy(n,n)*Lambda_xxy_S - vx(n,n)*Lambda_yxy_S) ) * diffFermi
                    Chi_yxxy_k(:,ieta, 1, 2) = Chi_yxxy_k(:,ieta, 1, 2) + real( ((G_xx_dky - G_xx)/dky - (G_yx_dkx - G_yx)/dkx)*sy(n,n) ) * diffFermi
                endif
                if (include_m_orb) then
                    Chi_xyyx_k(:,ieta, 2, 1) = Chi_xyyx_k(:,ieta, 2, 1) + real( (vx(n,n)*Lambda_yyx_L - vy(n,n)*Lambda_xyx_L) ) * diffFermi
                    Chi_xyyx_k(:,ieta, 2, 2) = Chi_xyyx_k(:,ieta, 2, 2) + real( ((G_yy_dkx - G_yy)/dkx - (G_xy_dky - G_xy)/dky)*lx(n,n) ) * diffFermi

                    Chi_yxxy_k(:,ieta, 2, 1) = Chi_yxxy_k(:,ieta, 2, 1) + real( (vy(n,n)*Lambda_xxy_L - vx(n,n)*Lambda_yxy_L) ) * diffFermi
                    Chi_yxxy_k(:,ieta, 2, 2) = Chi_yxxy_k(:,ieta, 2, 2) + real( ((G_xx_dky - G_xx)/dky - (G_yx_dkx - G_yx)/dkx)*ly(n,n) ) * diffFermi
                endif
            enddo ! ieta
        enddo ! n
    
    end subroutine


    subroutine drude_weight_single_k(k_in, drude_k)
        use para
        implicit none
    
        real(dp), intent(in)  :: k_in(3)
        real(dp), intent(out) :: drude_k(OmegaNum, Eta_number, 2)
                
        real(dp),    allocatable :: diffFermi(:)

        ! eigen value of H
        real(dp),    allocatable :: W(:)
        complex(dp), allocatable :: Hamk_bulk(:, :)
        complex(dp), allocatable :: UU(:, :)    
    
        complex(dp), allocatable :: velocities(:,:,:)
        complex(dp), allocatable :: vx(:, :), vy(:, :)

        allocate( diffFermi (OmegaNum))
        allocate( W         (Num_wann))
        allocate( Hamk_bulk (Num_wann, Num_wann))
        allocate( UU        (Num_wann, Num_wann))
        allocate( velocities(Num_wann, Num_wann, 3))
        allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann))
    
        call ham_bulk_latticegauge(k_in, Hamk_bulk)
        UU=Hamk_bulk
        call eigensystem_c( 'V', 'U', Num_wann, UU, W)
        call velocity_latticegauge_simple(k_in, UU, velocities)
        vx = velocities(:,:,1)
        vy = velocities(:,:,2)
    
        drude_k = 0d0
    
        do n= 1, Num_wann
            if (W(n)<OmegaMin- 2.d-2 .or. W(n)>OmegaMax+ 2.d-2) cycle
            do ieta=1, Eta_number
                !> this format is very important! prevent NaN error
                diffFermi = -1d0 / (Exp((W(n) - energy)/Eta_array(ieta))+1d0) / (Exp(-(W(n) - energy)/Eta_array(ieta))+1d0) / Eta_array(ieta)

                drude_k(:, ieta, 1) = drude_k(:, ieta, 1) + real(vx(n,n))**2 *diffFermi
                drude_k(:, ieta, 2) = drude_k(:, ieta, 2) + real(vy(n,n))**2 *diffFermi
            enddo ! ieta
        enddo ! n
    
    end subroutine drude_weight_single_k


    subroutine sigma_TRAHC_k(kin, sigma_tensor)

        use wmpi
        use para
        implicit none
    
        real(dp), intent(in)  :: kin(3)
        real(dp), intent(out) :: sigma_tensor(OmegaNum, 8, Eta_number)
    
        integer  :: iR, ikm
        real(dp) :: kdotr
        !real(dp) :: sqkvec1(3),sqkvec2(3),sqkvec3(3)
        !real(dp) :: kyp(3), kym(3)
        real(dp) :: kmat(3,3)
        real(dp) :: pyG_yy, pyG_yx, pyG_xx
        real(dp) :: pxG_xx, pxG_xy, pxG_yy
    
        !> Gij
        real(dp) :: G_xy(3),G_yx(3),G_xx(3),G_yy(3)
    
        !> eigen value of H
        real(dp), allocatable :: W(:),eig(:),eigmat(:,:)
        complex(dp), allocatable :: Hamk_bulk(:, :)
        complex(dp), allocatable :: Amat(:, :)
        complex(dp), allocatable :: UU(:, :)
        complex(dp), allocatable :: UU_dag(:, :)
    
        !> velocities
        complex(dp), allocatable :: vx(:, :), vy(:, :)
        complex(dp), allocatable :: vxmat(:, :, :), vymat(:, :, :)
        complex(dp), allocatable :: vxx(:, :), vxy(:, :), vyy(:, :)
    
        real(dp) :: vxx_2, vxy_2, vyy_2, exx, exy, eyy
        real(dp), allocatable :: Fshort(:) !> short notation of the denominator of the Fermi distribution
        real(dp), allocatable :: diffFermi(:), diff2Fermi(:)
    
        allocate( Fshort(OmegaNum), diffFermi(OmegaNum), diff2Fermi(OmegaNum))
    
        allocate( W (Num_wann))
        allocate( eig(Num_wann))
        allocate( eigmat(Num_wann,3))
        allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann))
        allocate( vxmat(Num_wann, Num_wann,3), vymat(Num_wann, Num_wann,3))
        allocate( vxx(Num_wann, Num_wann), vxy(Num_wann, Num_wann), vyy(Num_wann, Num_wann))
        allocate( Hamk_bulk(Num_wann, Num_wann))
        allocate( Amat(Num_wann, Num_wann))
        allocate( UU(Num_wann, Num_wann))
        allocate( UU_dag(Num_wann, Num_wann))
    
        Hamk_bulk=0d0
        Amat= 0d0
        UU_dag=0d0
        UU= 0d0
        eig=0d0
        eigmat=0d0
        W=0d0
    
        kmat(1,:)=kin - [Origin_cell%Rua(1) , Origin_cell%Rub(1) , Origin_cell%Ruc(1)]*dkx/twopi
        kmat(2,:)=kin - [Origin_cell%Rua(2) , Origin_cell%Rub(2) , Origin_cell%Ruc(2)]*dky/twopi
        kmat(3,:)=kin
    
        do ikm= 1,3
            k=kmat(ikm,:)
            ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
            call ham_bulk_latticegauge(k, Hamk_bulk)
            
            !> diagonalization by call zheev in lapack
            UU=Hamk_bulk
            call eigensystem_c( 'V', 'U', Num_wann, UU, W)
            UU_dag= conjg(transpose(UU))
            eigmat(:,ikm)=W
            vx= 0d0; vy= 0d0!; vz= 0d0
            vxx= 0d0; vxy =0d0; vyy=0d0
            do iR= 1, Nrpts
                kdotr= k(1)*irvec(1,iR) + k(2)*irvec(2,iR) + k(3)*irvec(3,iR)
                vx= vx+ zi*crvec(1, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
                vy= vy+ zi*crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
                if (ikm == 3) then
                    vxx= vxx - crvec(1, iR)*crvec(1, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
                    vxy= vxy - crvec(1, iR)*crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
                    vyy= vyy - crvec(2, iR)*crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
                endif
            enddo ! iR
        
            !> unitility rotate velocity
            UU_dag= conjg(transpose(UU))
            call mat_mul(Num_wann, vx, UU, Amat)
            call mat_mul(Num_wann, UU_dag, Amat, vx)
            call mat_mul(Num_wann, vy, UU, Amat)
            call mat_mul(Num_wann, UU_dag, Amat, vy)
            vxmat(:,:,ikm)=vx
            vymat(:,:,ikm)=vy
            if (ikm == 3) then
                call mat_mul(Num_wann, vxx, UU, Amat)
                call mat_mul(Num_wann, UU_dag, Amat, vxx)
                call mat_mul(Num_wann, vxy, UU, Amat)
                call mat_mul(Num_wann, UU_dag, Amat, vxy)
                call mat_mul(Num_wann, vyy, UU, Amat)
                call mat_mul(Num_wann, UU_dag, Amat, vyy)
            endif
        enddo !ikm
    
        eig= eigmat(:,3)
        vx = vxmat(:,:,3)
        vy = vymat(:,:,3)
        sigma_tensor=0d0
    
        do m= 1, Num_wann
            if (eig(m)<OmegaMin- 2.d-2 .or. eig(m)>OmegaMax+ 2.d-2) cycle !> prevent NaN error
            !> calculate G for each band
            !> G_xy == G_yx, we calculate both of them to check it
            G_xx=0d0; G_xy=0d0; G_yx=0d0; G_yy=0d0
            vxx_2=0d0; vxy_2=0d0; vyy_2=0d0
            
            !> sum of all other bands n
            do n= 1, Num_wann
                if (ABS(eig(m)-eig(n)) < band_degeneracy_threshold) cycle
    
                do ikm = 1,3
                    G_xx(ikm)= G_xx(ikm)+ 2.d0*real(vxmat(m, n, ikm)*vxmat(n, m, ikm)*((eigmat(m,ikm)-eigmat(n,ikm))/((eigmat(m,ikm)-eigmat(n,ikm))**2))**3)
                    G_xy(ikm)= G_xy(ikm)+ 2.d0*real(vxmat(m, n, ikm)*vymat(n, m, ikm)*((eigmat(m,ikm)-eigmat(n,ikm))/((eigmat(m,ikm)-eigmat(n,ikm))**2))**3)
                    G_yx(ikm)= G_yx(ikm)+ 2.d0*real(vymat(m, n, ikm)*vxmat(n, m, ikm)*((eigmat(m,ikm)-eigmat(n,ikm))/((eigmat(m,ikm)-eigmat(n,ikm))**2))**3)
                    G_yy(ikm)= G_yy(ikm)+ 2.d0*real(vymat(m, n, ikm)*vymat(n, m, ikm)*((eigmat(m,ikm)-eigmat(n,ikm))/((eigmat(m,ikm)-eigmat(n,ikm))**2))**3)
                enddo
    
                vxx_2= vxx_2+2d0*real(vx(m, n)*vx(n, m)*(eig(m)-eig(n))/((eig(m)-eig(n))**2))
                vxy_2= vxy_2+2d0*real(vx(m, n)*vy(n, m)*(eig(m)-eig(n))/((eig(m)-eig(n))**2))
                vyy_2= vyy_2+2d0*real(vy(m, n)*vy(n, m)*(eig(m)-eig(n))/((eig(m)-eig(n))**2))
            enddo ! n
        
            pyG_yy=(G_yy(3)-G_yy(2))/dky
            pyG_xx=(G_xx(3)-G_xx(2))/dky
            pyG_yx=(G_yx(3)-G_yx(2))/dky  
    
            pxG_xx=(G_xx(3)-G_xx(1))/dkx
            pxG_xy=(G_xy(3)-G_xy(1))/dkx
            pxG_yy=(G_yy(3)-G_yy(1))/dkx
    
            exx=real(vxx(m,m))+vxx_2
            exy=real(vxy(m,m))+vxy_2
            eyy=real(vyy(m,m))+vyy_2
            
            do ieta=1, Eta_number
                !> the outmost band index is m here
                Fshort = Exp((eig(m)-energy)/Eta_array(ieta))
                !> this format is very important! prevent NaN error
                diffFermi  = -1d0 / (Fshort+1d0) / (1d0/Fshort+1d0) / Eta_array(ieta)
                diff2Fermi = diffFermi * ( 1d0 - 2d0 / (1d0/Fshort+1d0)) / Eta_array(ieta)
    
                ! xxxx/tau
                sigma_tensor(:,1,ieta) = sigma_tensor(:,1,ieta) + pxG_xx*real(vx(m,m))*diffFermi + real(vx(m,m)*vx(m,m))*G_xx(3)*diff2Fermi/2d0
    
                ! xxyy/tau
                sigma_tensor(:,2,ieta) = sigma_tensor(:,2,ieta) +  real(vy(m,m))* (2d0*pxG_xy + pyG_xx) *diffFermi/3d0 &
                    + ( real(vx(m,m)*vx(m,m))*G_yy(3) + 2d0*real(vx(m,m)*vy(m,m))*G_xy(3) )*diff2Fermi/6d0
    
                ! yyxx/tau
                sigma_tensor(:,3,ieta) = sigma_tensor(:,3,ieta) +  real(vx(m,m))* (2d0*pyG_yx + pxG_yy) *diffFermi/3d0 &
                    + ( real(vy(m,m)*vy(m,m))*G_xx(3) + 2d0*real(vy(m,m)*vx(m,m))*G_yx(3) )*diff2Fermi/6d0
    
                ! yyyy/tau
                sigma_tensor(:,4,ieta) = sigma_tensor(:,4,ieta) + pyG_yy*real(vy(m,m))*diffFermi + real(vy(m,m)*vy(m,m))*G_yy(3)*diff2Fermi/2d0
    
                ! xxxx/tau**3
                sigma_tensor(:,5,ieta) = sigma_tensor(:,5,ieta) + exx*exx*diffFermi + exx*real(vx(m,m)*vx(m,m))*diff2Fermi
    
                ! xxyy/tau**3
                sigma_tensor(:,6,ieta) = sigma_tensor(:,6,ieta) + (exx*eyy*diffFermi + exx*real(vy(m,m)*vy(m,m))*diff2Fermi &
                    +2d0*exy*exy*diffFermi + 2d0*exy*real(vx(m,m)*vy(m,m))*diff2Fermi)/3d0
    
                ! yyxx/tau**3
                sigma_tensor(:,7,ieta) = sigma_tensor(:,7,ieta) + (eyy*exx*diffFermi + eyy*real(vx(m,m)*vx(m,m))*diff2Fermi &
                    +2d0*exy*exy*diffFermi + 2d0*exy*real(vx(m,m)*vy(m,m))*diff2Fermi)/3d0
                
                ! yyyy/tau**3
                sigma_tensor(:,8,ieta) = sigma_tensor(:,8,ieta) + eyy*eyy*diffFermi + eyy*real(vy(m,m)*vy(m,m))*diff2Fermi
                
            enddo ! ieta
        enddo !m
    
    end subroutine sigma_TRAHC_k    
end module


subroutine velocity_latticegauge_simple(k_in, UU, velocities) !> dH_dk, without 1/hbar
    use para, only: dp, irvec, crvec, HmnR, pi2zi, ndegen, Nrpts, Num_wann, zi
    implicit none

    real(dp),    intent(in)  :: k_in(3)
    complex(dp), intent(in)  :: UU(Num_wann, Num_wann)
    complex(dp), intent(out) :: velocities(Num_wann, Num_wann, 3)

    real(dp):: kdotr
    integer :: iR
    complex(dp), allocatable :: Amat(:, :), UU_dag(:,:)
    allocate( Amat(Num_wann, Num_wann), UU_dag(Num_wann, Num_wann))

    velocities= 0d0
    do iR= 1, Nrpts
        kdotr= k_in(1)*irvec(1,iR) + k_in(2)*irvec(2,iR) + k_in(3)*irvec(3,iR)
        velocities(:,:,1)= velocities(:,:,1) + zi*crvec(1, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
        velocities(:,:,2)= velocities(:,:,2) + zi*crvec(2, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
        velocities(:,:,3)= velocities(:,:,3) + zi*crvec(3, iR)*HmnR(:,:,iR)*Exp(pi2zi*kdotr)/ndegen(iR)
    enddo ! iR

    UU_dag= conjg(transpose(UU))
    !> unitility rotate velocity
    call mat_mul(Num_wann, velocities(:,:,1), UU, Amat)
    call mat_mul(Num_wann, UU_dag, Amat, velocities(:,:,1))
    call mat_mul(Num_wann, velocities(:,:,2), UU, Amat)
    call mat_mul(Num_wann, UU_dag, Amat, velocities(:,:,2))
    call mat_mul(Num_wann, velocities(:,:,3), UU, Amat)
    call mat_mul(Num_wann, UU_dag, Amat, velocities(:,:,3))

end subroutine velocity_latticegauge_simple


subroutine ik_to_kpoint(ik,k)
    use para, only: dp, Nk1, Nk2, Nk3, K3D_start_cube, K3D_vec1_cube, K3D_vec2_cube, K3D_vec3_cube
    implicit none
    integer(8),intent(in)  :: ik
    real(dp),  intent(out) :: k(3)
    integer(8)             :: ikx, iky, ikz

    ikx= (ik-1)/(Nk2*Nk3)+1
    iky= ((ik-1-(ikx-1)*Nk2*Nk3)/Nk3)+1
    ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
    k= K3D_start_cube + K3D_vec1_cube*(ikx-1)/dble(Nk1) + K3D_vec2_cube*(iky-1)/dble(Nk2) + K3D_vec3_cube*(ikz-1)/dble(Nk3)
end subroutine ik_to_kpoint


subroutine get_Fermi_energy_list(energy) 
    !> return Fermi energy in Hatree energy, not eV
    use para, only: dp, OmegaNum, OmegaMin, OmegaMax
    implicit none
    real(dp),  intent(inout) :: energy(OmegaNum)
    integer :: ie

    energy = 0d0
    do ie=1, OmegaNum
        if (OmegaNum>1) then
            energy(ie)= OmegaMin+ (OmegaMax-OmegaMin)* (ie- 1d0)/(OmegaNum- 1d0)
        else
            energy= OmegaMin
        endif
    enddo ! ie
end subroutine get_Fermi_energy_list

