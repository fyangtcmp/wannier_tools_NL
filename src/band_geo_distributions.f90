!> Calculate the distributions of some band geometry properties which are used in the nonlinear transport
!> Do not parallel these codes over multiple nodes, cause it does not support the OPENMP.
!> the calculations on nonlinear transport are heavy, so the parallel version of wt.x is needed.


subroutine band_geo_props_kplane
    
    use wmpi
    use para
    use nonlinear_transport
    implicit none
   
    integer  :: i, j, nkmesh

    real(dp) :: kxy_plane(3)

    !> k points slice
    real(dp), allocatable :: kslice(:, :), kslice_xyz(:, :)
  
    real(dp), allocatable :: props(:,:), props_mpi(:,:)

    nkmesh= Nk1*Nk2
    allocate( props    (nkmesh, 6))
    allocate( props_mpi(nkmesh, 6))

    allocate( kslice(nkmesh, 3))
    allocate( kslice_xyz(nkmesh, 3))

    props= 0d0
    props_mpi= 0d0

    kslice=0d0
    kslice_xyz=0d0
   
    ik =0
    do i= 1, nk1
        do j= 1, nk2
            ik =ik +1
            kslice(ik, :)= K3D_start+ K3D_vec1*(i-1)/dble(nk1-1)  &
                        + K3D_vec2*(j-1)/dble(nk2-1) - (K3D_vec1+ K3D_vec2)/2d0
            kslice_xyz(ik, :)= kslice(ik, 1)* Origin_cell%Kua+ kslice(ik, 2)* Origin_cell%Kub+ kslice(ik, 3)* Origin_cell%Kuc 
        enddo
    enddo

    call now(time_start)
    do ik= 1+ cpuid, nkmesh, num_cpu
        if (cpuid.eq.0 .and. mod(ik/num_cpu, 200).eq.0) then
            call now(time_end)
            write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/knv3', &
                ik, nkmesh, '  time left', (nkmesh-ik)*(time_end-time_start)/num_cpu/200d0/60d0
            time_start= time_end
        endif 

        k= kslice(ik, :)

        select case (option_prop)
        case ("ISOAHC")
            call ISOAHC_dist_single_k_Ef(k,  props_mpi(ik, 1:6))
        case ("INPHC")
            call INPHC_dist_single_k_Ef(k,  props_mpi(ik, 1:6))
        case ("SHC")
            call sigma_SHC_single_k_EF(k,  props_mpi(ik, 1:6))
        case Default
            write(stdout,*) "ERROR: option_prop was wrongly set or did not been provided, must be in ""ISOAHC"" ""INPHC"" ""SHC"" "
            stop
        end select

    enddo ! ik

#if defined (MPI)
    call mpi_allreduce(props_mpi, props, size(props), mpi_dp,mpi_sum,mpi_cmw,ierr)
#endif

    !> write the Berry curvature to file
    outfileindex= outfileindex+ 1
    if (cpuid==0) then
        select case (option_prop)
        case ("ISOAHC")
            open(unit=outfileindex, file='ISOAHC_kplane.dat')
            write(outfileindex, '("# sumover: Fermi ", a15)') option_sumover
            write(outfileindex, '("#")')
            write(outfileindex, '("#",1a11, 2a12, 15a18)') "kx(1/A)", "ky(1/A)", "kz(1/A)", 'G_{xx}', 'G_{xy}', 'G_{yx}', 'G_{yy}', 'Λ_{xyy}', 'Λ_{yxx}'
        case ("INPHC")
            open(unit=outfileindex, file='INPHC_kplane.dat')
            write(outfileindex, '("# sumover: Fermi ", a15)') option_sumover
            write(outfileindex, '("# unit: (e^2)/hbar * Angstorm^3 * (Ohm * V * T)^-1 ")') 
            write(outfileindex, '("#",1a11, 2a12, 17a20)') "kx(1/A)", "ky(1/A)", "kz(1/A)", '\Upsilon_{xyyy}', '\Upsilon_{yxxx}', '\Upsilon_{xyyx}', '\Upsilon_{yxxy}' 
        case ("SHC")
            open(unit=outfileindex, file='SHC_kplane.dat')
            write(outfileindex, '("# sumover: Fermi ", a15)') option_sumover
            write(outfileindex, '("#")')
            write(outfileindex, '("#",1a11, 2a12, 15a18)') "kx (1/A)", "ky (1/A)", "kz (1/A)", 'zx^x', 'zx^y', 'zx^z', 'prop4', 'prop5', 'prop6'           
        case Default
            ! close(outfileindex)
            return
        end select

        ik= 0
        do i= 1, nk1
            do j= 1, nk2
                ik= ik+ 1
                ! call rotate_k3_to_kplane(kslice_xyz(ik, :), kxy_plane)
                write(outfileindex, '(3f12.6, 12E15.4e3)') kslice_xyz(ik, :)*Ang2Bohr, props(ik, :) ! kxy_plane*Angstrom2atomic, &
            enddo
            ! write(outfileindex, *) ' '
        enddo

        close(outfileindex)

    endif

end subroutine


subroutine ISOAHC_dist_single_k_Ef(k_in, props)
    
    use nonlinear_transport
    use para
    implicit none
   
    real(dp), intent(in)  :: k_in(3)
    real(dp), intent(out) :: props(6)
    real(dp) :: diffFermi
    integer  :: sum_end

    ! eigen value of H
    real(dp),    allocatable :: W(:)
    complex(dp), allocatable :: Hamk_bulk(:, :)
    complex(dp), allocatable :: Amat(:, :)
    complex(dp), allocatable :: UU(:, :)
    complex(dp), allocatable :: UU_dag(:, :)

    complex(dp), allocatable :: vx(:, :), vy(:, :)
    complex(dp), allocatable :: velocities(:,:,:)

    real(dp) :: G_xy, G_yx, G_xx, G_yy

    allocate( W (Num_wann))
    allocate( Hamk_bulk (Num_wann, Num_wann))
    allocate( Amat (Num_wann, Num_wann))
    allocate( UU (Num_wann, Num_wann))
    allocate( UU_dag (Num_wann, Num_wann))

    allocate( velocities(Num_wann, Num_wann, 3))
    allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann))

    call ham_bulk_latticegauge(k_in, Hamk_bulk)
    UU=Hamk_bulk
    call eigensystem_c( 'V', 'U', Num_wann, UU, W)
    UU_dag= conjg(transpose(UU))
    call velocity_latticegauge_simple(k_in, UU, velocities)
    vx = velocities(:,:,1)
    vy = velocities(:,:,2)

    props = 0d0

    if (option_sumover=='sea') then ! 'sea'=Fermi sea
        sum_end = Numoccupied
    elseif (option_sumover=='surface') then ! 'surface'=Fermi surface
        sum_end = Num_wann
    else
        write(stdout,*) "ERROR: option_sumover was wrongly set, must be in ""sea"" ""surface"" "
        stop
    endif

    do n= 1, sum_end
        !> calculate G for each band
        G_xx=0d0; G_xy=0d0; G_yx=0d0; G_yy=0d0

        do m= 1, Num_wann
            if (ABS(W(m)-W(n)) < band_degeneracy_threshold) cycle
            G_xx= G_xx+ 2.d0*real(vx(n, m)*vx(m, n)/((W(n)-W(m))**3))
            G_xy= G_xy+ 2.d0*real(vx(n, m)*vy(m, n)/((W(n)-W(m))**3))
            G_yx= G_yx+ 2.d0*real(vy(n, m)*vx(m, n)/((W(n)-W(m))**3))
            G_yy= G_yy+ 2.d0*real(vy(n, m)*vy(m, n)/((W(n)-W(m))**3))
        enddo ! m
        
        if (option_sumover=='sea') then ! 'sea'=Fermi sea
            diffFermi = 1d0
        elseif (option_sumover=='surface') then ! 'surface'=Fermi surface
            diffFermi = -1d0 / (Exp((W(n) - E_arc)/Eta_Arc)+1d0) / (Exp(-(W(n) - E_arc)/Eta_Arc)+1d0) / Eta_Arc
        endif

        props(1) = props(1) + G_xx * diffFermi
        props(2) = props(2) + G_xy * diffFermi
        props(3) = props(3) + G_yx * diffFermi
        props(4) = props(4) + G_yy * diffFermi
        props(5) = props(5) + (G_yy*real(vx(n,n))-G_xy*real(vy(n,n))) * diffFermi
        props(6) = props(6) + (G_xx*real(vy(n,n))-G_yx*real(vx(n,n))) * diffFermi
        
    enddo ! n

    deallocate(W, vx, vy, Hamk_bulk, Amat, UU, UU_dag, velocities)
    return
end subroutine


subroutine INPHC_dist_single_k_Ef(k_in, props) 

    use nonlinear_transport
    use magnetic_moments
    use para
    implicit none

    real(dp), intent(in)  :: k_in(3)
    real(dp), intent(out) :: props(6)
    integer  :: sum_end

    real(dp), allocatable :: Chi_xyyy_k         (:,:,:,:)
    real(dp), allocatable :: Chi_yxxx_k         (:,:,:,:)
    real(dp), allocatable :: Chi_xyyx_k         (:,:,:,:)
    real(dp), allocatable :: Chi_yxxy_k         (:,:,:,:)

    complex(dp), allocatable :: M_S(:, :, :) !> spin magnetic moments
    complex(dp), allocatable :: M_L(:, :, :) !> orbital magnetic moments
        
    real(dp) :: k_dkx(3)
    real(dp) :: k_dky(3)
    
    real(dp), allocatable :: diffFermi(:)

    ! eigen value of H
    real(dp),    allocatable :: W(:)
    complex(dp), allocatable :: Hamk_bulk(:, :)
    complex(dp), allocatable :: Amat(:, :)
    complex(dp), allocatable :: UU(:, :)
    complex(dp), allocatable :: UU_dag(:, :)

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
    real(dp) :: Lambda_xxy_S, Lambda_yxy_S, Lambda_yyx_S, Lambda_xyx_S
    real(dp) :: Lambda_xxy_L, Lambda_yxy_L, Lambda_yyx_L, Lambda_xyx_L
    real(dp) :: G_xx, G_xy, G_yx, G_yy, G_yy_dkx, G_xy_dky, G_xx_dky, G_yx_dkx
    real(dp) :: dEnm, dEnm3, dEml, dEnl

    allocate( Chi_xyyy_k         (OmegaNum, Eta_number,2,2))
    allocate( Chi_yxxx_k         (OmegaNum, Eta_number,2,2))
    allocate( Chi_xyyx_k         (OmegaNum, Eta_number,2,2))
    allocate( Chi_yxxy_k         (OmegaNum, Eta_number,2,2))
    
    allocate( diffFermi (OmegaNum))

    !===========================================================================
    !> original kpoints
    allocate( W (Num_wann))
    allocate( Hamk_bulk (Num_wann, Num_wann))
    allocate( Amat (Num_wann, Num_wann))
    allocate( UU (Num_wann, Num_wann))
    allocate( UU_dag (Num_wann, Num_wann))

    allocate( velocities(Num_wann, Num_wann, 3))
    allocate( vx(Num_wann, Num_wann), vy(Num_wann, Num_wann))
    allocate( sx(Num_wann, Num_wann), sy(Num_wann, Num_wann))
    allocate( lx(Num_wann, Num_wann), ly(Num_wann, Num_wann))

    call ham_bulk_latticegauge(k_in, Hamk_bulk)
    UU=Hamk_bulk
    call eigensystem_c( 'V', 'U', Num_wann, UU, W)
    UU_dag= conjg(transpose(UU))
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
    Chi_xyyx_k = 0d0
    Chi_yxxy_k = 0d0

    if (option_sumover=='sea') then ! 'sea'=Fermi sea
        sum_end = Numoccupied
    elseif (option_sumover=='surface') then ! 'surface'=Fermi surface
        sum_end = Num_wann
    else
        write(stdout,*)  "ERROR: option_sumover was wrongly set, must be in ""sea"" ""surface"" "
        stop
    endif

    do n= 1, sum_end
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
                Lambda_xyy_S = Lambda_xyy_S + 6.d0* real( vx(n, m)*vy(m, n)*(sy(n, n)-sy(m, m)) )/dEnm3/dEnm
                Lambda_yyy_S = Lambda_yyy_S + 6.d0* real( vy(n, m)*vy(m, n)*(sy(n, n)-sy(m, m)) )/dEnm3/dEnm
                Lambda_yxx_S = Lambda_yxx_S + 6.d0* real( vy(n, m)*vx(m, n)*(sx(n, n)-sx(m, m)) )/dEnm3/dEnm
                Lambda_xxx_S = Lambda_xxx_S + 6.d0* real( vx(n, m)*vx(m, n)*(sx(n, n)-sx(m, m)) )/dEnm3/dEnm
                Lambda_xxy_S = Lambda_xxy_S + 6.d0* real( vx(n, m)*vx(m, n)*(sy(n, n)-sy(m, m)) )/dEnm3/dEnm
                Lambda_yxy_S = Lambda_yxy_S + 6.d0* real( vy(n, m)*vx(m, n)*(sy(n, n)-sy(m, m)) )/dEnm3/dEnm
                Lambda_yyx_S = Lambda_yyx_S + 6.d0* real( vy(n, m)*vy(m, n)*(sx(n, n)-sx(m, m)) )/dEnm3/dEnm
                Lambda_xyx_S = Lambda_xyx_S + 6.d0* real( vx(n, m)*vy(m, n)*(sx(n, n)-sx(m, m)) )/dEnm3/dEnm
                
                do l= 1, Num_wann
                    dEnl= W(n)-W(l)
                    dEml= W(m)-W(l)                    
                    if (ABS(dEnl) > band_degeneracy_threshold) then
                        Lambda_xyy_S = Lambda_xyy_S - 2.d0* real( (vx(l, m)*vy(m, n)+vy(l, m)*vx(m, n)) *sy(n, l)) /dEnm3/dEnl
                        Lambda_yyy_S = Lambda_yyy_S - 2.d0* real( (vy(l, m)*vy(m, n)+vy(l, m)*vy(m, n)) *sy(n, l)) /dEnm3/dEnl
                        Lambda_yxx_S = Lambda_yxx_S - 2.d0* real( (vy(l, m)*vx(m, n)+vx(l, m)*vy(m, n)) *sx(n, l)) /dEnm3/dEnl
                        Lambda_xxx_S = Lambda_xxx_S - 2.d0* real( (vx(l, m)*vx(m, n)+vx(l, m)*vx(m, n)) *sx(n, l)) /dEnm3/dEnl
                        Lambda_xxy_S = Lambda_xxy_S - 2.d0* real( (vx(l, m)*vx(m, n)*2                ) *sy(n, l)) /dEnm3/dEnl
                        Lambda_yxy_S = Lambda_yxy_S - 2.d0* real( (vy(l, m)*vx(m, n)+vx(l, m)*vy(m, n)) *sy(n, l)) /dEnm3/dEnl
                        Lambda_yyx_S = Lambda_yyx_S - 2.d0* real( (vy(l, m)*vy(m, n)*2                ) *sx(n, l)) /dEnm3/dEnl
                        Lambda_xyx_S = Lambda_xyx_S - 2.d0* real( (vx(l, m)*vy(m, n)+vy(l, m)*vx(m, n)) *sx(n, l)) /dEnm3/dEnl
                    endif
                    if (ABS(dEml) > band_degeneracy_threshold) then
                        Lambda_xyy_S = Lambda_xyy_S - 2.d0* real( (vx(l, n)*vy(n, m)+vy(l, n)*vx(n, m)) *sy(m, l)) /dEnm3/dEml
                        Lambda_yyy_S = Lambda_yyy_S - 2.d0* real( (vy(l, n)*vy(n, m)+vy(l, n)*vy(n, m)) *sy(m, l)) /dEnm3/dEml
                        Lambda_yxx_S = Lambda_yxx_S - 2.d0* real( (vy(l, n)*vx(n, m)+vx(l, n)*vy(n, m)) *sx(m, l)) /dEnm3/dEml
                        Lambda_xxx_S = Lambda_xxx_S - 2.d0* real( (vx(l, n)*vx(n, m)+vx(l, n)*vx(n, m)) *sx(m, l)) /dEnm3/dEml
                        Lambda_xxy_S = Lambda_xxy_S - 2.d0* real( (vx(l, n)*vx(n, m)*2                ) *sy(m, l)) /dEnm3/dEml
                        Lambda_yxy_S = Lambda_yxy_S - 2.d0* real( (vy(l, n)*vx(n, m)+vx(l, n)*vy(n, m)) *sy(m, l)) /dEnm3/dEml
                        Lambda_yyx_S = Lambda_yyx_S - 2.d0* real( (vy(l, n)*vy(n, m)*2                ) *sx(m, l)) /dEnm3/dEml
                        Lambda_xyx_S = Lambda_xyx_S - 2.d0* real( (vx(l, n)*vy(n, m)+vy(l, n)*vx(n, m)) *sx(m, l)) /dEnm3/dEml
                    endif
                enddo ! l
            endif

            if (include_m_orb) then
                Lambda_xyy_L = Lambda_xyy_L + 6.d0* real( vx(n, m)*vy(m, n)*(ly(n, n)-ly(m, m)) )/dEnm3/dEnm
                Lambda_yyy_L = Lambda_yyy_L + 6.d0* real( vy(n, m)*vy(m, n)*(ly(n, n)-ly(m, m)) )/dEnm3/dEnm
                Lambda_yxx_L = Lambda_yxx_L + 6.d0* real( vy(n, m)*vx(m, n)*(lx(n, n)-lx(m, m)) )/dEnm3/dEnm
                Lambda_xxx_L = Lambda_xxx_L + 6.d0* real( vx(n, m)*vx(m, n)*(lx(n, n)-lx(m, m)) )/dEnm3/dEnm
                Lambda_xxy_L = Lambda_xxy_L + 6.d0* real( vx(n, m)*vx(m, n)*(sy(n, n)-ly(m, m)) )/dEnm3/dEnm
                Lambda_yxy_L = Lambda_yxy_L + 6.d0* real( vy(n, m)*vx(m, n)*(sy(n, n)-ly(m, m)) )/dEnm3/dEnm
                Lambda_yyx_L = Lambda_yyx_L + 6.d0* real( vy(n, m)*vy(m, n)*(sx(n, n)-lx(m, m)) )/dEnm3/dEnm
                Lambda_xyx_L = Lambda_xyx_L + 6.d0* real( vx(n, m)*vy(m, n)*(sx(n, n)-lx(m, m)) )/dEnm3/dEnm
                
                do l= 1, Num_wann
                    dEnl= W(n)-W(l)
                    dEml= W(m)-W(l)                    
                    if (ABS(dEnl) > band_degeneracy_threshold) then
                        Lambda_xyy_L = Lambda_xyy_L - 2.d0* real( (vx(l, m)*vy(m, n)+vy(l, m)*vx(m, n)) *ly(n, l)) /dEnm3/dEnl
                        Lambda_yyy_L = Lambda_yyy_L - 2.d0* real( (vy(l, m)*vy(m, n)+vy(l, m)*vy(m, n)) *ly(n, l)) /dEnm3/dEnl
                        Lambda_yxx_L = Lambda_yxx_L - 2.d0* real( (vy(l, m)*vx(m, n)+vx(l, m)*vy(m, n)) *lx(n, l)) /dEnm3/dEnl
                        Lambda_xxx_L = Lambda_xxx_L - 2.d0* real( (vx(l, m)*vx(m, n)+vx(l, m)*vx(m, n)) *lx(n, l)) /dEnm3/dEnl
                        Lambda_xxy_L = Lambda_xxy_L - 2.d0* real( (vx(l, m)*vx(m, n)*2                ) *ly(n, l)) /dEnm3/dEnl
                        Lambda_yxy_L = Lambda_yxy_L - 2.d0* real( (vy(l, m)*vx(m, n)+vx(l, m)*vy(m, n)) *ly(n, l)) /dEnm3/dEnl
                        Lambda_yyx_L = Lambda_yyx_L - 2.d0* real( (vy(l, m)*vy(m, n)*2                ) *lx(n, l)) /dEnm3/dEnl
                        Lambda_xyx_L = Lambda_xyx_L - 2.d0* real( (vx(l, m)*vy(m, n)+vy(l, m)*vx(m, n)) *lx(n, l)) /dEnm3/dEnl
                    endif
                    if (ABS(dEml) > band_degeneracy_threshold) then
                        Lambda_xyy_L = Lambda_xyy_L - 2.d0* real( (vx(l, n)*vy(n, m)+vy(l, n)*vx(n, m)) *ly(m, l)) /dEnm3/dEml
                        Lambda_yyy_L = Lambda_yyy_L - 2.d0* real( (vy(l, n)*vy(n, m)+vy(l, n)*vy(n, m)) *ly(m, l)) /dEnm3/dEml
                        Lambda_yxx_L = Lambda_yxx_L - 2.d0* real( (vy(l, n)*vx(n, m)+vx(l, n)*vy(n, m)) *lx(m, l)) /dEnm3/dEml
                        Lambda_xxx_L = Lambda_xxx_L - 2.d0* real( (vx(l, n)*vx(n, m)+vx(l, n)*vx(n, m)) *lx(m, l)) /dEnm3/dEml
                        Lambda_xxy_L = Lambda_xxy_L - 2.d0* real( (vx(l, n)*vx(n, m)*2                ) *ly(m, l)) /dEnm3/dEml
                        Lambda_yxy_L = Lambda_yxy_L - 2.d0* real( (vy(l, n)*vx(n, m)+vx(l, n)*vy(n, m)) *ly(m, l)) /dEnm3/dEml
                        Lambda_yyx_L = Lambda_yyx_L - 2.d0* real( (vy(l, n)*vy(n, m)*2                ) *lx(m, l)) /dEnm3/dEml
                        Lambda_xyx_L = Lambda_xyx_L - 2.d0* real( (vx(l, n)*vy(n, m)+vy(l, n)*vx(n, m)) *lx(m, l)) /dEnm3/dEml
                    endif
                enddo ! l
            endif
        enddo ! m

        if (option_sumover=='sea') then ! 'sea'=Fermi sea
            diffFermi = 1d0
        elseif (option_sumover=='surface') then ! 'surface'=Fermi surface
            diffFermi = -1d0 / (Exp((W(n) - E_arc)/Eta_Arc)+1d0) / (Exp(-(W(n) - E_arc)/Eta_Arc)+1d0) / Eta_Arc
        endif

        ieta = 1
        if (include_m_spin) then
            Chi_xyyy_k(:,ieta, 1, 1) = Chi_xyyy_k(:,ieta, 1, 1) + real( (vx(n,n)*Lambda_yyy_S - vy(n,n)*Lambda_xyy_S) ) * diffFermi
            Chi_xyyy_k(:,ieta, 1, 2) = Chi_xyyy_k(:,ieta, 1, 2) + real( ((G_yy_dkx - G_yy)/dkx - (G_xy_dky - G_xy)/dky)*sy(n,n) ) * diffFermi

            Chi_yxxx_k(:,ieta, 1, 1) = Chi_yxxx_k(:,ieta, 1, 1) + real( (vy(n,n)*Lambda_xxx_S - vx(n,n)*Lambda_yxx_S) ) * diffFermi
            Chi_yxxx_k(:,ieta, 1, 2) = Chi_yxxx_k(:,ieta, 1, 2) + real( ((G_xx_dky - G_xx)/dky - (G_yx_dkx - G_yx)/dkx)*sx(n,n) ) * diffFermi

            Chi_xyyx_k(:,ieta, 1, 1) = Chi_xyyx_k(:,ieta, 1, 1) + real( (vx(n,n)*Lambda_yyx_S - vy(n,n)*Lambda_xyx_S) ) * diffFermi
            Chi_xyyx_k(:,ieta, 1, 2) = Chi_xyyx_k(:,ieta, 1, 2) + real( ((G_yy_dkx - G_yy)/dkx - (G_xy_dky - G_xy)/dky)*sx(n,n) ) * diffFermi

            Chi_yxxy_k(:,ieta, 1, 1) = Chi_yxxy_k(:,ieta, 1, 1) + real( (vy(n,n)*Lambda_xxy_S - vx(n,n)*Lambda_yxy_S) ) * diffFermi
            Chi_yxxy_k(:,ieta, 1, 2) = Chi_yxxy_k(:,ieta, 1, 2) + real( ((G_xx_dky - G_xx)/dky - (G_yx_dkx - G_yx)/dkx)*sy(n,n) ) * diffFermi
        endif
        if (include_m_orb) then
            Chi_xyyy_k(:,ieta, 2, 1) = Chi_xyyy_k(:,ieta, 2, 1) + real( (vx(n,n)*Lambda_yyy_L - vy(n,n)*Lambda_xyy_L) ) * diffFermi
            Chi_xyyy_k(:,ieta, 2, 2) = Chi_xyyy_k(:,ieta, 2, 2) + real( ((G_yy_dkx - G_yy)/dkx - (G_xy_dky - G_xy)/dky)*ly(n,n) ) * diffFermi

            Chi_yxxx_k(:,ieta, 2, 1) = Chi_yxxx_k(:,ieta, 2, 1) + real( (vy(n,n)*Lambda_xxx_L - vx(n,n)*Lambda_yxx_L) ) * diffFermi
            Chi_yxxx_k(:,ieta, 2, 2) = Chi_yxxx_k(:,ieta, 2, 2) + real( ((G_xx_dky - G_xx)/dky - (G_yx_dkx - G_yx)/dkx)*lx(n,n) ) * diffFermi

            Chi_xyyx_k(:,ieta, 2, 1) = Chi_xyyx_k(:,ieta, 2, 1) + real( (vx(n,n)*Lambda_yyx_L - vy(n,n)*Lambda_xyx_L) ) * diffFermi
            Chi_xyyx_k(:,ieta, 2, 2) = Chi_xyyx_k(:,ieta, 2, 2) + real( ((G_yy_dkx - G_yy)/dkx - (G_xy_dky - G_xy)/dky)*lx(n,n) ) * diffFermi

            Chi_yxxy_k(:,ieta, 2, 1) = Chi_yxxy_k(:,ieta, 2, 1) + real( (vy(n,n)*Lambda_xxy_L - vx(n,n)*Lambda_yxy_L) ) * diffFermi
            Chi_yxxy_k(:,ieta, 2, 2) = Chi_yxxy_k(:,ieta, 2, 2) + real( ((G_xx_dky - G_xx)/dky - (G_yx_dkx - G_yx)/dkx)*ly(n,n) ) * diffFermi
        endif
    enddo ! n

    props(1) = sum(Chi_xyyy_k(1,1,:,:))
    props(2) = sum(Chi_yxxx_k(1,1,:,:))
    props(3) = sum(Chi_xyyx_k(1,1,:,:))
    props(4) = sum(Chi_yxxy_k(1,1,:,:))
    props(5) = 0d0
    props(6) = 0d0

    props = props * (Echarge/Hartree2J * mu_B) / (Ang2Bohr)**3

end subroutine


subroutine sigma_SHC_single_k_EF(k_in, props)

    use wmpi
    use para
    use magnetic_moments
    implicit none
    
    real(dp), intent(in)  :: k_in(3)
    real(dp), intent(out) :: props(6)

    integer :: m, n, ialpha, ibeta, igamma

    real(dp) :: deno_fac

    ! eigen value of H
    real(dp), allocatable :: W(:)
    complex(dp), allocatable :: Hamk_bulk(:, :)
    complex(dp), allocatable :: UU(:, :)

    !> sigma^gamma_{alpha, beta}, alpha, beta, gamma=1,2,3 for x, y, z
    !>  sigma_tensor_shc(igamma, ialpha, ibeta)
    real(dp), allocatable :: sigma_tensor_shc(:, :, :)
    
    !> Fermi-Dirac distribution
    real(dp), external :: fermi

    complex(dp), allocatable :: Vmn_Ham(:, :, :)
    complex(dp), allocatable :: Vmn_wann(:, :, :)
    complex(dp), allocatable :: spin_sigma(:, :)
    complex(dp), allocatable :: j_spin_gamma_alpha(:, :)
    complex(dp), allocatable :: mat_t(:, :)

    !> Berry curvature vectors for all bands
    real(dp),allocatable :: Omega_spin(:)

    ! spin operator matrix spin_sigma_x,spin_sigma_y in spin_sigma_z representation
    complex(Dp),allocatable :: pauli_matrices(:, :, :) 
    
    allocate(Vmn_Ham(Num_wann, Num_wann, 3))
    allocate(Vmn_wann(Num_wann, Num_wann, 3))
    allocate(spin_sigma(Num_wann, Num_wann))
    allocate(j_spin_gamma_alpha(Num_wann, Num_wann))
    allocate(Omega_spin(Num_wann))

    allocate( W (Num_wann))
    allocate( Hamk_bulk(Num_wann, Num_wann))
    allocate( UU(Num_wann, Num_wann))
    allocate( mat_t(Num_wann, Num_wann))
    allocate( sigma_tensor_shc   (3, 3, 3))

    allocate(pauli_matrices(Num_wann, Num_wann, 3))
    spin_sigma= 0d0
    sigma_tensor_shc    = 0d0
    Hamk_bulk=0d0
    UU= 0d0
    Vmn_wann= 0d0
    pauli_matrices= 0d0

    call spin_operators(pauli_matrices)

    ! calculation bulk hamiltonian by a direct Fourier transformation of HmnR
    call ham_bulk_atomicgauge(k_in, Hamk_bulk)
    ! call ham_bulk_latticegauge(k, Hamk_bulk)

    !> diagonalization by call zheev in lapack
    UU= Hamk_bulk
    call eigensystem_c( 'V', 'U', Num_wann, UU, W)

    !> get velocity operator in Wannier basis
    !> \partial_k H_nm
    call dHdk_atomicgauge(k_in, Vmn_wann)

    do ialpha= 1, 3
        call rotation_to_Ham_basis(UU, Vmn_wann(:, :, ialpha), Vmn_Ham(:, :, ialpha))
    enddo

    !> spin axis igamma= x, y, z
    do igamma= 1, 3
        !> set Pauli matrix
        spin_sigma= pauli_matrices(:, :, igamma)

        !> calculate spin current operator j_spin_gamma_alpha^l_alpha= 1/2*{Sigma_gamma, v_alpha} 
        ! in order to calculate SHC^z_xy
        do ialpha= 1, 3
            !> in Wannier basis
            j_spin_gamma_alpha= 0d0
            call mat_mul(Num_wann, spin_sigma, Vmn_wann(:, :, ialpha), j_spin_gamma_alpha(:, :))
            call mat_mul(Num_wann, Vmn_wann(:, :, ialpha), spin_sigma, mat_t)
            j_spin_gamma_alpha(:, :)= j_spin_gamma_alpha(:, :)+ mat_t
            j_spin_gamma_alpha= j_spin_gamma_alpha/2d0
        
            !> rotate to Hamiltonian basis
            mat_t= j_spin_gamma_alpha(:, :)
            call rotation_to_Ham_basis(UU, mat_t, j_spin_gamma_alpha(:, :))

            !> \Omega_spin^l_n^{\gamma}(k)=-2\sum_{m}*aimag(Im({js(\gamma),v(\alpha)}/2)_nm*v_beta_mn))/((w(n)-w(m))^2+eta_arc^2)
            do ibeta= 1, 3
                Omega_spin= 0d0
                do n= 1, Num_wann
                    do m= 1, Num_wann
                        if (abs(W(n)-W(m)) < band_degeneracy_threshold) cycle
                        deno_fac= -2d0/((W(n)-W(m))**2)
                        Omega_spin(n)= Omega_spin(n)+ &
                            aimag(j_spin_gamma_alpha(n, m)*Vmn_Ham(m, n, ibeta))*deno_fac*fermi(W(n)-E_arc, 1d0/Eta_Arc)
                    enddo
                enddo

                !> sum over all "spin" Berry curvature below chemical potential mu
                sigma_tensor_shc(igamma, ialpha, ibeta)= sum(Omega_spin(:))
            enddo ! ibeta  v
        enddo ! ialpha  j
    enddo ! igamma  spin

    !> in the latest version, we use the atomic unit
    !> in unit of ((hbar/e)(Ohm*m)^-1
    sigma_tensor_shc = sigma_tensor_shc * Echarge**2/hbar/Bohr_radius/2d0

    props      = 0d0
    props(1:3) = sigma_tensor_shc(1:3 ,3 ,1 )
   
    
end subroutine 