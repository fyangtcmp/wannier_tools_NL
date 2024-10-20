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

    integer :: ready = 1

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


    subroutine sigma_NPHC_int_single_k(k_in, Chi_xyyy_k, Chi_yxxx_k, Chi_xyyx_k, Chi_yxxy_k)
        use magnetic_moments
        use para
        implicit none
    
        real(dp), intent(in)  :: k_in(3)
        real(dp), intent(out) :: Chi_xyyy_k(OmegaNum, Eta_number, 2, 2) !> the third index: 1=spin, 2=orbital
        real(dp), intent(out) :: Chi_yxxx_k(OmegaNum, Eta_number, 2, 2)
        real(dp), intent(out) :: Chi_xyyx_k(OmegaNum, Eta_number, 2, 2) !> the third index: 1=spin, 2=orbital
        real(dp), intent(out) :: Chi_yxxy_k(OmegaNum, Eta_number, 2, 2)

        integer :: L1, L2, L3
        real(dp) :: Lambda_S(2,2,2)
        real(dp) :: Lambda_L(2,2,2)
    
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
    
        complex(dp), allocatable :: vx(:, :), vy(:, :)
        complex(dp), allocatable :: velocities(:,:,:)
    
        complex(dp), allocatable :: vx_dkx(:, :), vy_dkx(:, :)  
        complex(dp), allocatable :: vx_dky(:, :), vy_dky(:, :) 
    
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
        endif
        if (include_m_orb) then
            call orbital_magnetic_moments(W, velocities, M_L)
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

            Lambda_S = 0d0
            Lambda_L = 0d0
    
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
    
                do L1 = 1,2
                    do L2 = 1,2
                        do L3 = 1,2
                            Lambda_S(L1,L2,L3) = Lambda_S(L1,L2,L3) + 6.d0* real( velocities(n,m,L1)*velocities(m,n,L2)*(M_S(n,n,L3)-M_S(m,m,L3)) )/dEnm3/dEnm
                            Lambda_L(L1,L2,L3) = Lambda_L(L1,L2,L3) + 6.d0* real( velocities(n,m,L1)*velocities(m,n,L2)*(M_L(n,n,L3)-M_L(m,m,L3)) )/dEnm3/dEnm
                        enddo
                    enddo
                enddo
                    
                    do l= 1, Num_wann
                        dEnl= W(n)-W(l)
                        dEml= W(m)-W(l)                
                        if (ABS(dEnl) > band_degeneracy_threshold) then
                            do L1 = 1,2
                                do L2 = 1,2
                                    do L3 = 1,2
                                        Lambda_S(L1,L2,L3) = Lambda_S(L1,L2,L3) - 2.d0* real( (velocities(l,m,L1)*velocities(m,n,L2) + velocities(l,m,L2)*velocities(m,n,L1))* M_S(n,l,L3) )/dEnm3/dEnl
                                        Lambda_L(L1,L2,L3) = Lambda_L(L1,L2,L3) - 2.d0* real( (velocities(l,m,L1)*velocities(m,n,L2) + velocities(l,m,L2)*velocities(m,n,L1))* M_L(n,l,L3) )/dEnm3/dEnl
                                    enddo
                                enddo
                            enddo
                        endif
                        if (ABS(dEml) > band_degeneracy_threshold) then
                            do L1 = 1,2
                                do L2 = 1,2
                                    do L3 = 1,2
                                        Lambda_S(L1,L2,L3) = Lambda_S(L1,L2,L3) - 2.d0* real( (velocities(l,n,L1)*velocities(n,m,L2) + velocities(l,n,L2)*velocities(n,m,L1)) * M_S(m,l,L3) )/dEnm3/dEml
                                        Lambda_L(L1,L2,L3) = Lambda_L(L1,L2,L3) - 2.d0* real( (velocities(l,n,L1)*velocities(n,m,L2) + velocities(l,n,L2)*velocities(n,m,L1)) * M_L(m,l,L3) )/dEnm3/dEml
                                    enddo
                                enddo
                            enddo
                        endif

                    enddo ! l
            enddo ! m

            do ieta=1, Eta_number
                !> the outmost band index is n here
                !> this format is very important! prevent NaN error
                diffFermi = -1d0 / (Exp((W(n) - energy)/Eta_array(ieta))+1d0) / (Exp(-(W(n) - energy)/Eta_array(ieta))+1d0) / Eta_array(ieta)

                if (include_m_spin) then
                    Chi_xyyy_k(:,ieta, 1, 1) = Chi_xyyy_k(:,ieta, 1, 1) + ( real(velocities(n,n,1))*Lambda_S(2,2,2) - real(velocities(n,n,2))*Lambda_S(1,2,2) ) * diffFermi
                    Chi_xyyy_k(:,ieta, 1, 2) = Chi_xyyy_k(:,ieta, 1, 2) + ((G_yy_dkx - G_yy)/dkx - (G_xy_dky - G_xy)/dky)*real(M_S(n,n,2)) * diffFermi

                    Chi_yxxx_k(:,ieta, 1, 1) = Chi_yxxx_k(:,ieta, 1, 1) + ( real(velocities(n,n,2))*Lambda_S(1,1,1) - real(velocities(n,n,1))*Lambda_S(2,1,1) ) * diffFermi
                    Chi_yxxx_k(:,ieta, 1, 2) = Chi_yxxx_k(:,ieta, 1, 2) + ((G_xx_dky - G_xx)/dky - (G_yx_dkx - G_yx)/dkx)*real(M_S(n,n,1)) * diffFermi

                    Chi_xyyx_k(:,ieta, 1, 1) = Chi_xyyx_k(:,ieta, 1, 1) + ( real(velocities(n,n,1))*Lambda_S(2,2,1) - real(velocities(n,n,2))*Lambda_S(1,2,1) ) * diffFermi
                    Chi_xyyx_k(:,ieta, 1, 2) = Chi_xyyx_k(:,ieta, 1, 2) + ((G_yy_dkx - G_yy)/dkx - (G_xy_dky - G_xy)/dky)*real(M_S(n,n,1)) * diffFermi

                    Chi_yxxy_k(:,ieta, 1, 1) = Chi_yxxy_k(:,ieta, 1, 1) + ( real(velocities(n,n,2))*Lambda_S(1,1,2) - real(velocities(n,n,1))*Lambda_S(2,1,2) ) * diffFermi
                    Chi_yxxy_k(:,ieta, 1, 2) = Chi_yxxy_k(:,ieta, 1, 2) + ((G_xx_dky - G_xx)/dky - (G_yx_dkx - G_yx)/dkx)*real(M_S(n,n,2)) * diffFermi
                endif
                if (include_m_orb) then
                    Chi_xyyy_k(:,ieta, 2, 1) = Chi_xyyy_k(:,ieta, 2, 1) + ( real(velocities(n,n,1))*Lambda_L(2,2,2) - real(velocities(n,n,2))*Lambda_L(1,2,2) ) * diffFermi
                    Chi_xyyy_k(:,ieta, 2, 2) = Chi_xyyy_k(:,ieta, 2, 2) + ((G_yy_dkx - G_yy)/dkx - (G_xy_dky - G_xy)/dky)*real(M_L(n,n,2)) * diffFermi

                    Chi_yxxx_k(:,ieta, 2, 1) = Chi_yxxx_k(:,ieta, 2, 1) + ( real(velocities(n,n,2))*Lambda_L(1,1,1) - real(velocities(n,n,1))*Lambda_L(2,1,1) ) * diffFermi
                    Chi_yxxx_k(:,ieta, 2, 2) = Chi_yxxx_k(:,ieta, 2, 2) + ((G_xx_dky - G_xx)/dky - (G_yx_dkx - G_yx)/dkx)*real(M_L(n,n,1)) * diffFermi

                    Chi_xyyx_k(:,ieta, 2, 1) = Chi_xyyx_k(:,ieta, 2, 1) + ( real(velocities(n,n,1))*Lambda_L(2,2,1) - real(velocities(n,n,2))*Lambda_L(1,2,1) ) * diffFermi
                    Chi_xyyx_k(:,ieta, 2, 2) = Chi_xyyx_k(:,ieta, 2, 2) + ((G_yy_dkx - G_yy)/dkx - (G_xy_dky - G_xy)/dky)*real(M_L(n,n,1)) * diffFermi

                    Chi_yxxy_k(:,ieta, 2, 1) = Chi_yxxy_k(:,ieta, 2, 1) + ( real(velocities(n,n,2))*Lambda_L(1,1,2) - real(velocities(n,n,1))*Lambda_L(2,1,2) ) * diffFermi
                    Chi_yxxy_k(:,ieta, 2, 2) = Chi_yxxy_k(:,ieta, 2, 2) + ((G_xx_dky - G_xx)/dky - (G_yx_dkx - G_yx)/dkx)*real(M_L(n,n,2)) * diffFermi
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


subroutine drude_weight

    use wmpi
    use para
    use nonlinear_transport
    implicit none
    
    real(dp), allocatable :: drude      (:,:,:)
    real(dp), allocatable :: drude_mpi  (:,:,:)
    real(dp), allocatable :: drude_k    (:,:,:)

    allocate( drude    (OmegaNum, Eta_number,2))
    allocate( drude_mpi(OmegaNum, Eta_number,2))
    allocate( drude_k  (OmegaNum, Eta_number,2))

    Eta_array(1) = Eta_Arc !> from wt.in
    Eta_array(2:Eta_number) = Eta_array_fixed(1:Eta_number-1)

    allocate( energy(OmegaNum))
    call get_Fermi_energy_list(energy)

    knv3= int8(Nk1)*Nk2*Nk3
    drude = 0d0
    drude_mpi = 0d0

    call now(time_start)
    do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0 .and. mod(ik/num_cpu, 2000).eq.0) then
            call now(time_end)
            write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/knv3', &
                ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/2000d0/60d0
            time_start= time_end
        endif

        call ik_to_kpoint(ik,k)

        call drude_weight_single_k(k, drude_k)

        drude_mpi = drude_mpi + drude_k
    enddo ! ik

#if defined (MPI)
    call mpi_barrier(mpi_cmw,ierr)
    call mpi_reduce(drude_mpi, drude, size(drude_mpi), mpi_dp, mpi_sum, 0, mpi_cmw, ierr)
    call mpi_reduce(drude_mpi, drude, size(drude_mpi), mpi_dp, mpi_sum, 0, mpi_cmw, ierr)
#endif

    drude = - drude * Echarge**2/hbar**2 *Hartree2J/Bohr_radius /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume

    outfileindex= outfileindex+ 1
    if (cpuid.eq.0) then
        do ieta=1, Eta_number
            write(Eta_name, '(f12.2)') Eta_array(ieta)*1000d0/eV2Hartree

            write(ahcfilename, '(7a)')'drude_weight_eta', trim(adjustl(Eta_name)), 'meV.dat'
            open(unit=outfileindex, file=ahcfilename)
            write(outfileindex, '("#",a)')' Drude weight, in unit of S/m/(relaxation time)'
            write(outfileindex, '("#",a)')'  '

            write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', 'xx', 'yy' !, 'zz'
            do ie=1, OmegaNum
                write(outfileindex, '(20E16.5e3)')energy(ie)/eV2Hartree, drude(ie,ieta,1), drude(ie,ieta,2)
            enddo
            close(outfileindex)
        enddo
    endif

end subroutine drude_weight


subroutine sigma_NPHC_int ! dynamical mpi, auto adapted k-mesh
    !> Calculate the intrinsic nonlinear planar Hall conductivity, the xyyy and yxxx elements
    !
    !> usage: sigma_NPHC_int_calc = T
    !
    !> ref : 10.1103/PhysRevLett.130.126303
    !
    !> 2023/10/31 Fan Yang
    !

    use wmpi
    use para
    use nonlinear_transport
    implicit none

    real(dp), parameter :: NPHC_int_unit_factor = Echarge**3/hbar/Hartree2J

    real(dp), allocatable :: Chi_xyyy_k         (:,:,:,:)
    real(dp), allocatable :: Chi_yxxx_k         (:,:,:,:)
    real(dp), allocatable :: Chi_xyyy_tensor    (:,:,:,:)
    real(dp), allocatable :: Chi_yxxx_tensor    (:,:,:,:)
    real(dp), allocatable :: Chi_xyyy_tensor_mpi(:,:,:,:)
    real(dp), allocatable :: Chi_yxxx_tensor_mpi(:,:,:,:)

    real(dp), allocatable :: Chi_xyyx_k         (:,:,:,:)
    real(dp), allocatable :: Chi_yxxy_k         (:,:,:,:)
    real(dp), allocatable :: Chi_xyyx_tensor    (:,:,:,:)
    real(dp), allocatable :: Chi_yxxy_tensor    (:,:,:,:)
    real(dp), allocatable :: Chi_xyyx_tensor_mpi(:,:,:,:)
    real(dp), allocatable :: Chi_yxxy_tensor_mpi(:,:,:,:)

    real(dp) :: max_tmp(4)

    Eta_array(1) = Eta_Arc !> from wt.in
    Eta_array(2:Eta_number) = Eta_array_fixed(1:Eta_number-1)

    allocate( Chi_xyyy_k         (OmegaNum, Eta_number,2,2))
    allocate( Chi_yxxx_k         (OmegaNum, Eta_number,2,2))
    allocate( Chi_xyyy_tensor    (OmegaNum, Eta_number,2,2))
    allocate( Chi_yxxx_tensor    (OmegaNum, Eta_number,2,2))
    allocate( Chi_xyyy_tensor_mpi(OmegaNum, Eta_number,2,2))
    allocate( Chi_yxxx_tensor_mpi(OmegaNum, Eta_number,2,2))

    allocate( Chi_xyyx_k         (OmegaNum, Eta_number,2,2))
    allocate( Chi_yxxy_k         (OmegaNum, Eta_number,2,2))
    allocate( Chi_xyyx_tensor    (OmegaNum, Eta_number,2,2))
    allocate( Chi_yxxy_tensor    (OmegaNum, Eta_number,2,2))
    allocate( Chi_xyyx_tensor_mpi(OmegaNum, Eta_number,2,2))
    allocate( Chi_yxxy_tensor_mpi(OmegaNum, Eta_number,2,2))

    Chi_xyyy_tensor     = 0d0
    Chi_yxxx_tensor     = 0d0
    Chi_xyyy_tensor_mpi = 0d0
    Chi_yxxx_tensor_mpi = 0d0

    Chi_xyyx_tensor     = 0d0
    Chi_yxxy_tensor     = 0d0
    Chi_xyyx_tensor_mpi = 0d0
    Chi_yxxy_tensor_mpi = 0d0

    allocate( energy(OmegaNum))
    call get_Fermi_energy_list(energy)

    !> each k point
    knv3= int8(Nk1)*Nk2*Nk3

    call get_k_fine_list()

    allocate( sk_list_mpi(knv3))
    allocate( sk_list    (knv3))
    sk_list_mpi = 0
    sk_list = 0

#if defined (MPI)
    if (cpuid==0) then ! dispatcher
        call now(time_start)
        do ik= 1, (knv3+num_cpu-1)
            if (mod(ik, 2000*(num_cpu-1))==0) then
                call now(time_end)
                write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/knv3', &
                    ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/(num_cpu-1)/2000d0/60d0
                time_start= time_end
            endif
    
            call MPI_Recv(ready, 1, mpi_in, mpi_any_source,       0, mpi_cmw, mpistatus, ierr)
            call MPI_Send(ik,    1, mpi_in, mpistatus.MPI_SOURCE, 0, mpi_cmw, ierr)       
        enddo ! ik
    else ! workers
        ! loop until all the kpoints have been scanned
        do while (.True.)
            call MPI_Send(ready, 1, mpi_in, 0, 0, mpi_cmw, ierr)
            call MPI_Recv(ik,    1, mpi_in, 0, 0, mpi_cmw, mpistatus, ierr)

            if (ik>knv3) exit
            call ik_to_kpoint(ik,k)
            call sigma_NPHC_int_single_k(k, Chi_xyyy_k, Chi_yxxx_k, Chi_xyyx_k, Chi_yxxy_k)

            max_tmp(1) = maxval(abs(Chi_xyyy_k))
            max_tmp(2) = maxval(abs(Chi_yxxx_k))
            max_tmp(3) = maxval(abs(Chi_xyyx_k))
            max_tmp(4) = maxval(abs(Chi_yxxy_k))
            sk_list_mpi(ik) = maxval( max_tmp )
        enddo
    endif

    call mpi_barrier(mpi_cmw,ierr)
    call mpi_reduce(sk_list_mpi,sk_list,size(sk_list),mpi_dp,mpi_sum,0,mpi_cmw,ierr)

    if (cpuid==0) then ! dispatcher
        call get_ik_adapt_list()
        call now(time_start)
        do ik2= 1, (Nk_adapt+num_cpu-1)
            if (mod(ik2, 200*(num_cpu-1))==0) then
                call now(time_end)
                write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/Nk_adapt', &
                    ik2, Nk_adapt, '  time left', (Nk_adapt-ik2)*(time_end-time_start)/(num_cpu-1)/200d0/60d0
                time_start= time_end
            endif

            if (ik2>Nk_adapt) then !> exit workers
                ik = knv3 + 1
            else
                ik = ik_adapt_list(ik2)
            endif
                
            call MPI_Recv(ready, 1, mpi_in, mpi_any_source,    0, mpi_cmw, mpistatus, ierr)    
            call MPI_Send(ik,    1, mpi_in, mpistatus.MPI_SOURCE, 0, mpi_cmw, ierr)  
        enddo ! ik

    else ! workers
        do while (.True.)
            call MPI_Send(ready, 1, mpi_in, 0, 0, mpi_cmw, ierr)  
            call MPI_Recv(ik,    1, mpi_in, 0, 0, mpi_cmw, mpistatus, ierr)
            !> MPI_sendrecv slower than individual MPI_send and MPI_recv
            ! call MPI_sendrecv(ready, 1, mpi_in, 0, 0, ik, 1, mpi_in, 0, 0, mpi_cmw, mpistatus, ierr) 

            if (ik>knv3) exit
            call ik_to_kpoint(ik,k)
            do ikfine=1, knv3_fine
                call sigma_NPHC_int_single_k(k + k_fine_list(ikfine,:), Chi_xyyy_k, Chi_yxxx_k, Chi_xyyx_k, Chi_yxxy_k)
    
                Chi_xyyy_tensor_mpi = Chi_xyyy_tensor_mpi + Chi_xyyy_k/dble(knv3_fine)
                Chi_yxxx_tensor_mpi = Chi_yxxx_tensor_mpi + Chi_yxxx_k/dble(knv3_fine)
                Chi_xyyx_tensor_mpi = Chi_xyyx_tensor_mpi + Chi_xyyx_k/dble(knv3_fine)
                Chi_yxxy_tensor_mpi = Chi_yxxy_tensor_mpi + Chi_yxxy_k/dble(knv3_fine)
            enddo
        enddo
    endif

    call mpi_reduce(Chi_xyyy_tensor_mpi, Chi_xyyy_tensor, size(Chi_xyyy_tensor), mpi_dp,mpi_sum, 0, mpi_cmw,ierr)
    call mpi_reduce(Chi_yxxx_tensor_mpi, Chi_yxxx_tensor, size(Chi_yxxx_tensor), mpi_dp,mpi_sum, 0, mpi_cmw,ierr)
    call mpi_reduce(Chi_xyyx_tensor_mpi, Chi_xyyx_tensor, size(Chi_xyyx_tensor), mpi_dp,mpi_sum, 0, mpi_cmw,ierr)
    call mpi_reduce(Chi_yxxy_tensor_mpi, Chi_yxxy_tensor, size(Chi_yxxy_tensor), mpi_dp,mpi_sum, 0, mpi_cmw,ierr)
#endif

    if (cpuid==0) then
        Chi_xyyy_tensor = Chi_xyyy_tensor * NPHC_int_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume
        Chi_yxxx_tensor = Chi_yxxx_tensor * NPHC_int_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume
        Chi_xyyx_tensor = Chi_xyyx_tensor * NPHC_int_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume
        Chi_yxxy_tensor = Chi_yxxy_tensor * NPHC_int_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume

        if (Nk3==1) then
            Chi_xyyy_tensor = Chi_xyyy_tensor* Origin_cell%Ruc(3)/Ang2Bohr
            Chi_yxxx_tensor = Chi_yxxx_tensor* Origin_cell%Ruc(3)/Ang2Bohr
            Chi_xyyx_tensor = Chi_xyyx_tensor* Origin_cell%Ruc(3)/Ang2Bohr
            Chi_yxxy_tensor = Chi_yxxy_tensor* Origin_cell%Ruc(3)/Ang2Bohr
        endif

        outfileindex= outfileindex+ 1
        do ieta=1, Eta_number
            write(Eta_name, '(f12.2)') Eta_array(ieta)*1000d0/eV2Hartree

            if (include_m_spin) then
                write(ahcfilename, '(7a)')'sigma_NPHC_int_S_eta', trim(adjustl(Eta_name)), 'meV.dat'
                open(unit=outfileindex, file=ahcfilename)
                write(outfileindex, '("#",a)')' Intrinsic nonlinear planar hall effect, in unit of A*V^-2*T^-1 for 3D case, Ang*A*V^-2*T^-1 for 2D cases'
                write(outfileindex, '("#",a)')' Please refer to the Sec. III of the supplementary materials of 10.1103/PhysRevLett.130.126303, for the definition of term I and term II of the INPHE conductivities'

                write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', '\chi_{xyyy,I}', '\chi_{xyyy,II}', '\chi_{xyyy}', '\chi_{yxxx,I}', '\chi_{yxxx,II}','\chi_{yxxx}', '\chi_{xyyx,I}', '\chi_{xyyx,II}', '\chi_{xyyx}', '\chi_{yxxy,I}', '\chi_{yxxy,II}','\chi_{yxxy}'
                do ie=1, OmegaNum
                    write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, &
                        Chi_xyyy_tensor(ie,ieta,1,1), Chi_xyyy_tensor(ie,ieta,1,2), Chi_xyyy_tensor(ie,ieta,1,1) + Chi_xyyy_tensor(ie,ieta,1,2),&
                        Chi_yxxx_tensor(ie,ieta,1,1), Chi_yxxx_tensor(ie,ieta,1,2), Chi_yxxx_tensor(ie,ieta,1,1) + Chi_yxxx_tensor(ie,ieta,1,2),&
                        Chi_xyyx_tensor(ie,ieta,1,1), Chi_xyyx_tensor(ie,ieta,1,2), Chi_xyyx_tensor(ie,ieta,1,1) + Chi_xyyx_tensor(ie,ieta,1,2),&
                        Chi_yxxy_tensor(ie,ieta,1,1), Chi_yxxy_tensor(ie,ieta,1,2), Chi_yxxy_tensor(ie,ieta,1,1) + Chi_yxxy_tensor(ie,ieta,1,2)
                enddo
                close(outfileindex)
            endif

            if (include_m_orb ) then
                write(ahcfilename, '(7a)')'sigma_NPHC_int_L_eta', trim(adjustl(Eta_name)), 'meV.dat'
                open(unit=outfileindex, file=ahcfilename)
                write(outfileindex, '("#",a)')' Intrinsic nonlinear planar hall effect, in unit of A*V^-2*T^-1 for 3D case, Ang*A*V^-2*T^-1 for 2D cases'
                write(outfileindex, '("#",a)')' Please refer to the Sec. III of the supplementary materials of 10.1103/PhysRevLett.130.126303, for the definition of term I and term II of the INPHE conductivities'

                write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', '\chi_{xyyy,I}', '\chi_{xyyy,II}', '\chi_{xyyy}', '\chi_{yxxx,I}', '\chi_{yxxx,II}','\chi_{yxxx}', '\chi_{xyyx,I}', '\chi_{xyyx,II}', '\chi_{xyyx}', '\chi_{yxxy,I}', '\chi_{yxxy,II}','\chi_{yxxy}'
                do ie=1, OmegaNum
                    write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, &
                        Chi_xyyy_tensor(ie,ieta,2,1), Chi_xyyy_tensor(ie,ieta,2,2), Chi_xyyy_tensor(ie,ieta,2,1) + Chi_xyyy_tensor(ie,ieta,2,2),&
                        Chi_yxxx_tensor(ie,ieta,2,1), Chi_yxxx_tensor(ie,ieta,2,2), Chi_yxxx_tensor(ie,ieta,2,1) + Chi_yxxx_tensor(ie,ieta,2,2),&
                        Chi_xyyx_tensor(ie,ieta,2,1), Chi_xyyx_tensor(ie,ieta,2,2), Chi_xyyx_tensor(ie,ieta,2,1) + Chi_xyyx_tensor(ie,ieta,2,2),&
                        Chi_yxxy_tensor(ie,ieta,2,1), Chi_yxxy_tensor(ie,ieta,2,2), Chi_yxxy_tensor(ie,ieta,2,1) + Chi_yxxy_tensor(ie,ieta,2,2)
                enddo
                close(outfileindex)
            endif
        enddo
    endif

end subroutine


subroutine sigma_NPHC_tau2 ! dynamical mpi, auto adapted k-mesh
    !> Calculate the tau^2 nonlinear planar Hall conductivity, the xyyy and yxxx elements
    !
    !> ref: DOI: 10.1103/PhysRevB.108.075155, eq(18)
    !
    !> usage: sigma_NPHC_tau2_calc = T
    !
    !> 2024/05/30 Fan Yang
    !

    use wmpi
    use para
    use nonlinear_transport
    implicit none

    real(dp), parameter :: NPHC_tau2_unit_factor = Echarge**3/hbar**3 *Hartree2J

    real(dp), allocatable :: Chi_xyyy_k         (:,:,:,:)
    real(dp), allocatable :: Chi_yxxx_k         (:,:,:,:)
    real(dp), allocatable :: Chi_xyyy_tensor    (:,:,:,:)
    real(dp), allocatable :: Chi_yxxx_tensor    (:,:,:,:)
    real(dp), allocatable :: Chi_xyyy_tensor_mpi(:,:,:,:)
    real(dp), allocatable :: Chi_yxxx_tensor_mpi(:,:,:,:)

    real(dp) :: max_tmp(2)

    Eta_array(1) = Eta_Arc !> from wt.in
    Eta_array(2:Eta_number) = Eta_array_fixed(1:Eta_number-1)

    allocate( Chi_xyyy_k         (OmegaNum, Eta_number,2,2))
    allocate( Chi_yxxx_k         (OmegaNum, Eta_number,2,2))
    allocate( Chi_xyyy_tensor    (OmegaNum, Eta_number,2,2))
    allocate( Chi_yxxx_tensor    (OmegaNum, Eta_number,2,2))
    allocate( Chi_xyyy_tensor_mpi(OmegaNum, Eta_number,2,2))
    allocate( Chi_yxxx_tensor_mpi(OmegaNum, Eta_number,2,2))

    Chi_xyyy_tensor     = 0d0
    Chi_yxxx_tensor     = 0d0
    Chi_xyyy_tensor_mpi = 0d0
    Chi_yxxx_tensor_mpi = 0d0

    allocate( energy(OmegaNum))
    call get_Fermi_energy_list(energy)

    !> each k point
    knv3= int8(Nk1)*Nk2*Nk3

    call get_k_fine_list()

    allocate( sk_list_mpi(knv3))
    allocate( sk_list    (knv3))
    sk_list_mpi = 0
    sk_list = 0

#if defined (MPI)
    if (cpuid==0) then ! dispatcher
        call now(time_start)
        do ik= 1, (knv3+num_cpu-1)
            if (mod(ik, 2000*(num_cpu-1))==0) then
                call now(time_end)
                write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/knv3', &
                    ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/(num_cpu-1)/2000d0/60d0
                time_start= time_end
            endif
    
            call MPI_Recv(ready, 1, mpi_in, mpi_any_source,       0, mpi_cmw, mpistatus, ierr)
            call MPI_Send(ik,    1, mpi_in, mpistatus.MPI_SOURCE, 0, mpi_cmw, ierr)       
        enddo ! ik
    else ! workers
        ! loop until all the kpoints have been scanned
        do while (.True.)
            call MPI_Send(ready, 1, mpi_in, 0, 0, mpi_cmw, ierr)
            call MPI_Recv(ik,    1, mpi_in, 0, 0, mpi_cmw, mpistatus, ierr)

            if (ik>knv3) exit
            call ik_to_kpoint(ik,k)
            call sigma_NPHC_tau2_single_k(k, Chi_xyyy_k, Chi_yxxx_k)

            max_tmp(1) = maxval(abs(Chi_xyyy_k))
            max_tmp(2) = maxval(abs(Chi_yxxx_k))
            sk_list_mpi(ik) = maxval( max_tmp )
        enddo
    endif

    call mpi_barrier(mpi_cmw,ierr)
    call mpi_reduce(sk_list_mpi,sk_list,size(sk_list),mpi_dp,mpi_sum,0,mpi_cmw,ierr)

    if (cpuid==0) then ! dispatcher
        call get_ik_adapt_list()
        call now(time_start)
        do ik2= 1, (Nk_adapt+num_cpu-1)
            if (mod(ik2, 200*(num_cpu-1))==0) then
                call now(time_end)
                write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/Nk_adapt', &
                    ik2, Nk_adapt, '  time left', (Nk_adapt-ik2)*(time_end-time_start)/(num_cpu-1)/200d0/60d0
                time_start= time_end
            endif

            if (ik2>Nk_adapt) then !> exit workers
                ik = knv3 + 1
            else
                ik = ik_adapt_list(ik2)
            endif
                
            call MPI_Recv(ready, 1, mpi_in, mpi_any_source,    0, mpi_cmw, mpistatus, ierr)    
            call MPI_Send(ik,    1, mpi_in, mpistatus.MPI_SOURCE, 0, mpi_cmw, ierr)  
        enddo ! ik

    else ! workers
        do while (.True.)
            call MPI_Send(ready, 1, mpi_in, 0, 0, mpi_cmw, ierr)  
            call MPI_Recv(ik,    1, mpi_in, 0, 0, mpi_cmw, mpistatus, ierr)
            !> MPI_sendrecv slower than individual MPI_send and MPI_recv
            ! call MPI_sendrecv(ready, 1, mpi_in, 0, 0, ik, 1, mpi_in, 0, 0, mpi_cmw, mpistatus, ierr) 

            if (ik>knv3) exit
            call ik_to_kpoint(ik,k)
            do ikfine=1, knv3_fine
                call sigma_NPHC_tau2_single_k(k + k_fine_list(ikfine,:), Chi_xyyy_k, Chi_yxxx_k)
    
                Chi_xyyy_tensor_mpi = Chi_xyyy_tensor_mpi + Chi_xyyy_k/dble(knv3_fine)
                Chi_yxxx_tensor_mpi = Chi_yxxx_tensor_mpi + Chi_yxxx_k/dble(knv3_fine)
            enddo
        enddo
    endif

    call mpi_reduce(Chi_xyyy_tensor_mpi, Chi_xyyy_tensor, size(Chi_xyyy_tensor), mpi_dp,mpi_sum, 0, mpi_cmw,ierr)
    call mpi_reduce(Chi_yxxx_tensor_mpi, Chi_yxxx_tensor, size(Chi_yxxx_tensor), mpi_dp,mpi_sum, 0, mpi_cmw,ierr)
#endif

    if (cpuid==0) then
        Chi_xyyy_tensor = Chi_xyyy_tensor * NPHC_tau2_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume
        Chi_yxxx_tensor = Chi_yxxx_tensor * NPHC_tau2_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume

        if (Nk3==1) then
            Chi_xyyy_tensor = Chi_xyyy_tensor* Origin_cell%Ruc(3)/Ang2Bohr
            Chi_yxxx_tensor = Chi_yxxx_tensor* Origin_cell%Ruc(3)/Ang2Bohr
        endif

        outfileindex= outfileindex+ 1
        do ieta=1, Eta_number
            write(Eta_name, '(f12.2)') Eta_array(ieta)*1000d0/eV2Hartree

            if (include_m_spin) then
                write(ahcfilename, '(7a)')'sigma_NPHC_tau2_S_eta', trim(adjustl(Eta_name)), 'meV.dat'
                open(unit=outfileindex, file=ahcfilename)
                write(outfileindex, '("#",a)')' tau2 nonlinear planar hall effect, in unit of A*V^-2*T^-1 for 3D case, Ang*A*V^-2*T^-1 for 2D cases'
                write(outfileindex, '("#",a)')' Without Tau'

                write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', '\sigma_xyyy_I', '\sigma_xyyy_II', '\sigma_xyyy_tot', '\sigma_yxxx_I', '\sigma_yxxx_II','\sigma_yxxx_tot'
                do ie=1, OmegaNum
                    write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, &
                        Chi_xyyy_tensor(ie,ieta,1,1), Chi_xyyy_tensor(ie,ieta,1,2), Chi_xyyy_tensor(ie,ieta,1,1) + Chi_xyyy_tensor(ie,ieta,1,2),&
                        Chi_yxxx_tensor(ie,ieta,1,1), Chi_yxxx_tensor(ie,ieta,1,2), Chi_yxxx_tensor(ie,ieta,1,1) + Chi_yxxx_tensor(ie,ieta,1,2)
                enddo
                close(outfileindex)
            endif

            if (include_m_orb ) then
                write(ahcfilename, '(7a)')'sigma_NPHC_tau2_L_eta', trim(adjustl(Eta_name)), 'meV.dat'
                open(unit=outfileindex, file=ahcfilename)
                write(outfileindex, '("#",a)')' tau2 nonlinear planar hall effect, in unit of A*V^-2*T^-1 for 3D case, Ang*A*V^-2*T^-1 for 2D cases'
                write(outfileindex, '("#",a)')' Without Tau'

                write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', '\sigma_xyyy_I', '\sigma_xyyy_II', '\sigma_xyyy_tot', '\sigma_yxxx_I', '\sigma_yxxx_II','\sigma_yxxx_tot'
                do ie=1, OmegaNum
                    write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, &
                        Chi_xyyy_tensor(ie,ieta,2,1), Chi_xyyy_tensor(ie,ieta,2,2), Chi_xyyy_tensor(ie,ieta,2,1) + Chi_xyyy_tensor(ie,ieta,2,2),&
                        Chi_yxxx_tensor(ie,ieta,2,1), Chi_yxxx_tensor(ie,ieta,2,2), Chi_yxxx_tensor(ie,ieta,2,1) + Chi_yxxx_tensor(ie,ieta,2,2)
                enddo
                close(outfileindex)
            endif
        enddo
    endif

end subroutine


subroutine sigma_TRAHC ! dynamical mpi, auto adapted k-mesh
    ! 对于对称性约化的带括号的项，因子3已经乘进去和1/3抵消了，计算电流时不需要再乘3了
    !> DOI: 10.1103/PhysRevB.105.045118

    use wmpi
    use para
    use nonlinear_transport
    implicit none

    real(dp), parameter :: TRAHC_factor_tau1 = Echarge**4/hbar**2 *Bohr_radius /Hartree2J
    real(dp), parameter :: TRAHC_factor_tau3 = Echarge**4/hbar**4 *Bohr_radius *Hartree2J
    ! real(dp), parameter :: sigma_threshold = 1d7

    real(dp), allocatable :: sigma_k         (:,:,:)   !> the second index = xxxx xxyy yyxx yyyy
    real(dp), allocatable :: sigma_tensor    (:,:,:)
    real(dp), allocatable :: sigma_tensor_mpi(:,:,:)

    allocate( energy(OmegaNum))

    allocate( sigma_k            (OmegaNum, 8, Eta_number))
    allocate( sigma_tensor       (OmegaNum, 8, Eta_number))
    allocate( sigma_tensor_mpi   (OmegaNum, 8, Eta_number))
    
    sigma_tensor     = 0d0
    sigma_tensor_mpi = 0d0

    call get_Fermi_energy_list(energy)

    Eta_array(1) = Eta_Arc !> from wt.in
    Eta_array(2:Eta_number) = Eta_array_fixed(1:Eta_number-1)

    !> each k point
    knv3= int8(Nk1)*Nk2*Nk3

    call get_k_fine_list()

    allocate( sk_list_mpi(knv3))
    allocate( sk_list(knv3) )
    sk_list_mpi = 0
    sk_list     = 0

#if defined (MPI)
    if (cpuid==0) then ! dispatcher
        call now(time_start)
        do ik= 1, (knv3+num_cpu-1)
            if (mod(ik, 200*(num_cpu-1))==0) then
                call now(time_end)
                write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/knv3', &
                    ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/(num_cpu-1)/200d0/60d0
                time_start= time_end
            endif
    
            call MPI_Recv(ready, 1, mpi_in, mpi_any_source,       0, mpi_cmw, mpistatus, ierr)
            call MPI_Send(ik,    1, mpi_in, mpistatus.MPI_SOURCE, 0, mpi_cmw, ierr)       
        enddo ! ik
    else ! workers
        ! loop until all the kpoints have been scanned
        do while (.True.)
            call MPI_Send(ready, 1, mpi_in, 0, 0, mpi_cmw, ierr)
            call MPI_Recv(ik,    1, mpi_in, 0, 0, mpi_cmw, mpistatus, ierr)

            if (ik>knv3) exit
            call ik_to_kpoint(ik,k)
            call sigma_TRAHC_k(k, sigma_k)

            sk_list_mpi(ik) = maxval(abs(sigma_k(:,1:4,:)))
        enddo
    endif

    call mpi_barrier(mpi_cmw,ierr)
    call mpi_reduce(sk_list_mpi,sk_list,size(sk_list),mpi_dp,mpi_sum,0,mpi_cmw,ierr)

    if (cpuid==0) then ! dispatcher
        call get_ik_adapt_list()
        call now(time_start)
        do ik2= 1, (Nk_adapt+num_cpu-1)
            if (mod(ik2, 200*(num_cpu-1))==0) then
                call now(time_end)
                write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/Nk_adapt', &
                    ik2, Nk_adapt, '  time left', (Nk_adapt-ik2)*(time_end-time_start)/(num_cpu-1)/200d0/60d0
                time_start= time_end
            endif

            if (ik2>Nk_adapt) then !> exit workers
                ik = knv3 + 1
            else
                ik = ik_adapt_list(ik2)
            endif
                
            call MPI_Recv(ready, 1, mpi_in, mpi_any_source,    0, mpi_cmw, mpistatus, ierr)    
            call MPI_Send(ik,    1, mpi_in, mpistatus.MPI_SOURCE, 0, mpi_cmw, ierr)  
        enddo ! ik

    else ! workers
        do while (.True.)
            call MPI_Send(ready, 1, mpi_in, 0, 0, mpi_cmw, ierr)  
            call MPI_Recv(ik,    1, mpi_in, 0, 0, mpi_cmw, mpistatus, ierr)
            !> MPI_sendrecv slower than individual MPI_send and MPI_recv
            ! call MPI_sendrecv(ready, 1, mpi_in, 0, 0, ik, 1, mpi_in, 0, 0, mpi_cmw, mpistatus, ierr) 

            if (ik>knv3) exit
            call ik_to_kpoint(ik,k)
            do ikfine=1, knv3_fine
                call sigma_TRAHC_k( k + k_fine_list(ikfine,:), sigma_k )
                sigma_tensor_mpi = sigma_tensor_mpi + sigma_k/dble(knv3_fine)
            enddo
        enddo
    endif

    call mpi_barrier(mpi_cmw,ierr)
    call mpi_reduce(sigma_tensor_mpi,sigma_tensor,size(sigma_tensor),mpi_dp,mpi_sum, 0, mpi_cmw,ierr)
#endif

    if (cpuid==0) then
        sigma_tensor(:,1:4,:)= sigma_tensor(:,1:4,:)/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume*TRAHC_factor_tau1 
        sigma_tensor(:,5:8,:)= sigma_tensor(:,5:8,:)/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume*TRAHC_factor_tau3
        
        if (Nk3==1) then
            sigma_tensor(:,1:4,:)= sigma_tensor(:,1:4,:)* Origin_cell%Ruc(3)/Ang2Bohr
            sigma_tensor(:,5:8,:)= sigma_tensor(:,5:8,:)* Origin_cell%Ruc(3)/Ang2Bohr
        endif

        outfileindex= outfileindex+ 1
        do ieta=1, Eta_number
            write(Eta_name, '(f12.2)') Eta_array(ieta)*1000d0/eV2Hartree
            write(ahcfilename, '(7a)')'sigma_TRAHC_eta', trim(adjustl(Eta_name)), 'meV.dat'
            open(unit=outfileindex, file=ahcfilename)
            write(outfileindex, '("#",a)')' Third-order anomalous hall conductivity (without tau)'
            write(outfileindex, '("#",a)')' in unit of A*V^-3 for 3D cases, Ang*A*V^-3 for 2D cases (with tau)'
            write(outfileindex, '("#",a13, 20a16)') 'Energy (eV) ', 'xxxx/tau', ' xxyy/tau', ' yyxx/tau',  ' yyyy/tau', &
                'xxxx/tau**3', ' xxyy/tau**3', ' yyxx/tau**3',  ' yyyy/tau**3'
            do ie=1, OmegaNum
                write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, &
                    sigma_tensor(ie,1,ieta), sigma_tensor(ie,2,ieta), sigma_tensor(ie,3,ieta), sigma_tensor(ie,4,ieta),&
                    sigma_tensor(ie,5,ieta), sigma_tensor(ie,6,ieta), sigma_tensor(ie,7,ieta), sigma_tensor(ie,8,ieta)
            enddo
            close(outfileindex)
        enddo
    endif

end subroutine


subroutine sigma_SOAHC_int ! static mpi, fixed k-mesh

    !> Calculate the intrinsic second order anomalous hall conductivity, the xyy and yxx elements
    !
    !> usage: sigma_SOAHC_int_calc = T
    !
    !> ref1 : 10.1103/PhysRevLett.127.277201
    !> ref2 : 10.1103/PhysRevLett.127.277202
    !
    !> Original developed by Huiying Liu
    !> 2022/07/15 Fan Yang, correct the units
    !> 2023/10/30 Fan Yang, update to wannier tools 2.7.0

    use wmpi
    use para
    use nonlinear_transport
    implicit none

    real(dp), parameter :: SOAHC_unit_factor =  Echarge**3/hbar/Hartree2J

    real(dp), allocatable :: sigma_xyy    (:,:)
    real(dp), allocatable :: sigma_yxx    (:,:)
    real(dp), allocatable :: sigma_xyy_k  (:,:)
    real(dp), allocatable :: sigma_yxx_k  (:,:)
    real(dp), allocatable :: sigma_xyy_mpi(:,:)
    real(dp), allocatable :: sigma_yxx_mpi(:,:)

    Eta_array(1) = Eta_Arc !> from wt.in
    Eta_array(2:Eta_number) = Eta_array_fixed(1:Eta_number-1)

    allocate( sigma_xyy        (OmegaNum, Eta_number))
    allocate( sigma_yxx        (OmegaNum, Eta_number))
    allocate( sigma_xyy_k      (OmegaNum, Eta_number))
    allocate( sigma_yxx_k      (OmegaNum, Eta_number))
    allocate( sigma_xyy_mpi    (OmegaNum, Eta_number))
    allocate( sigma_yxx_mpi    (OmegaNum, Eta_number))

    allocate( energy(OmegaNum))
    call get_Fermi_energy_list(energy)

    knv3= int8(Nk1)*Nk2*Nk3

    sigma_xyy_mpi    = 0.d0
    sigma_yxx_mpi    = 0.d0

    call now(time_start)
    do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0.and. mod(ik/num_cpu, 2000).eq.0) then
            call now(time_end)
            write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/knv3', &
                ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/2000d0/60d0
            time_start= time_end
        endif

        call ik_to_kpoint(ik,k)

        call sigma_SOAHC_int_single_k(k, sigma_xyy_k, sigma_yxx_k)

        sigma_xyy_mpi = sigma_xyy_mpi + sigma_xyy_k
        sigma_yxx_mpi = sigma_yxx_mpi + sigma_yxx_k
    enddo ! ik

#if defined (MPI)
    call mpi_barrier(mpi_cmw,ierr)
    call mpi_allreduce(sigma_xyy_mpi,sigma_xyy,size(sigma_xyy),&
        mpi_dp,mpi_sum,mpi_cmw,ierr)
    call mpi_allreduce(sigma_yxx_mpi,sigma_yxx,size(sigma_yxx),&
        mpi_dp,mpi_sum,mpi_cmw,ierr)
#endif

    !> the sigma_xyy contains an additional [energy]^-1 dimension, so besides e^3/hbar, we need to convert hartree to joule
    sigma_xyy= sigma_xyy * SOAHC_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume
    sigma_yxx= sigma_yxx * SOAHC_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume
    
    if (Nk3==1) then
        sigma_xyy= sigma_xyy * Origin_cell%Ruc(3)/Ang2Bohr
        sigma_yxx= sigma_yxx * Origin_cell%Ruc(3)/Ang2Bohr
    endif

    outfileindex= outfileindex+ 1
    if (cpuid.eq.0) then
        do ieta=1, Eta_number
            write(Eta_name, '(f12.2)') Eta_array(ieta)*1000d0/eV2Hartree
            write(ahcfilename, '(7a)')'sigma_SOAHC_int_eta', trim(adjustl(Eta_name)), 'meV.dat'
            open(unit=outfileindex, file=ahcfilename)
            write(outfileindex, '("#",a)')' Intrinsic 2nd anomalous hall conductivity'
            write(outfileindex, '("#",a)')' in unit of A*V^-2 for 3D cases, Ang*A*V^-2 for 2D cases'
            write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', '\sigma_xyy', '\sigma_yxx'
            do ie=1, OmegaNum
                write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, sigma_xyy(ie,ieta), &
                    sigma_yxx(ie,ieta)
            enddo
            close(outfileindex)
        enddo
    endif

end subroutine sigma_SOAHC_int