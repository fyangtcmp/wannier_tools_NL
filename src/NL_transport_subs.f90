!> subroutines based on the nonlinear_transport.mod
!> they are outside the module so that they can be directly called by the main program

subroutine sigma_TRAHC_static
    ! 对于对称性约化的带括号的项，因子3已经乘进去和1/3抵消了，计算电流时不需要再乘3了
    use wmpi
    use para
    use nonlinear_transport
    implicit none
    
    real(dp), parameter :: TRAHC_factor_tau1 = Echarge**4/hbar**2 *Bohr_radius /Hartree2J
    real(dp), parameter :: TRAHC_factor_tau3 = Echarge**4/hbar**4 *Bohr_radius *Hartree2J

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
    allocate( sk_list    (knv3))
    sk_list_mpi = 0
    sk_list = 0

    call now(time_start)
    do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0.and. mod(ik/num_cpu, 2000).eq.0) then
            call now(time_end)
            write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/knv3', &
                ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/2000d0/60d0
            time_start= time_end
        endif
        
        call ik_to_kpoint(ik,k)
        call sigma_TRAHC_k(k, sigma_k)
        sk_list_mpi(ik) = maxval(abs(sigma_k(:,1:4,:)))
    enddo ! ik

#if defined (MPI)
    call mpi_barrier(mpi_cmw,ierr)
    call mpi_allreduce(sk_list_mpi,sk_list,size(sk_list),mpi_dp,mpi_sum,mpi_cmw,ierr)
#endif

    call get_ik_adapt_list()
    deallocate(sk_list, sk_list_mpi)
    
    call now(time_start)
    do ik2= 1+ cpuid, Nk_adapt, num_cpu
        if (cpuid.eq.0.and. mod(ik2/num_cpu, 200).eq.0) then
            call now(time_end)
            write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/Nk_adapt', &
            ik2, Nk_adapt, '  time left', (Nk_adapt-ik2)*(time_end-time_start)/num_cpu/200d0/60d0
            time_start= time_end
        endif
        ik = ik_adapt_list(ik2)
        call ik_to_kpoint(ik,k)
 
        do ikfine=1, knv3_fine
            call sigma_TRAHC_k( k + k_fine_list(ikfine,:), sigma_k )
            sigma_tensor_mpi = sigma_tensor_mpi + sigma_k/dble(knv3_fine)
        enddo
        
    enddo ! ik


#if defined (MPI)
    call mpi_barrier(mpi_cmw,ierr)
    call mpi_allreduce(sigma_tensor_mpi,sigma_tensor,size(sigma_tensor),&
                    mpi_dp,mpi_sum,mpi_cmw,ierr)
#endif
 
    sigma_tensor(:,1:4,:)= sigma_tensor(:,1:4,:)/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume*TRAHC_factor_tau1 
    sigma_tensor(:,5:8,:)= sigma_tensor(:,5:8,:)/dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume*TRAHC_factor_tau3 
    
    if (Nk3==1) then
        sigma_tensor(:,1:4,:)= sigma_tensor(:,1:4,:)* Origin_cell%Ruc(3)/Ang2Bohr
        sigma_tensor(:,5:8,:)= sigma_tensor(:,5:8,:)* Origin_cell%Ruc(3)/Ang2Bohr
    endif

    outfileindex= outfileindex+ 1
    if (cpuid.eq.0) then
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


subroutine sigma_TRAHC ! dynamical mpi version
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


subroutine sigma_SOAHC_int

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
    !> 2023/11/06 Fan Yang, adaptive k-meshes methods to accelerate speed

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


subroutine sigma_NPHC_int_static
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

    knv3= int8(Nk1)*Nk2*Nk3

    call get_k_fine_list()

    allocate( sk_list_mpi(knv3))
    allocate( sk_list    (knv3))

    sk_list_mpi = 0
    sk_list = 0

    call now(time_start)
    do ik= 1+ cpuid, knv3, num_cpu
        if (cpuid.eq.0 .and. mod(ik/num_cpu, 2000).eq.0) then
            call now(time_end)
            write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/knv3', &
                ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/num_cpu/2000d0/60d0
            time_start= time_end
        endif

        call ik_to_kpoint(ik,k)

        call sigma_NPHC_int_single_k(k, Chi_xyyy_k, Chi_yxxx_k)

        max_tmp(1) = maxval(abs(Chi_xyyy_k))
        max_tmp(2) = maxval(abs(Chi_yxxx_k))
        sk_list_mpi(ik) = maxval( max_tmp ) 
    enddo ! ik

    !--------------------------------------------------------------------------
#if defined (MPI)
    call mpi_barrier(mpi_cmw,ierr)
    call mpi_allreduce(sk_list_mpi,sk_list,size(sk_list),mpi_dp,mpi_sum,mpi_cmw,ierr)
#endif

    call get_ik_adapt_list()
    deallocate(sk_list, sk_list_mpi)
    
    call now(time_start)
    do ik2= 1+ cpuid, Nk_adapt, num_cpu
        if (cpuid.eq.0.and. mod(ik2/num_cpu, 200).eq.0) then
            call now(time_end)
            write(stdout, '(a, i18, "/", i18, a, f10.2, "min")') 'ik/Nk_adapt', &
            ik2, Nk_adapt, '  time left', (Nk_adapt-ik2)*(time_end-time_start)/num_cpu/200d0/60d0
            time_start= time_end
        endif
        ik = ik_adapt_list(ik2)
        call ik_to_kpoint(ik,k)
 
        do ikfine=1, knv3_fine
            call sigma_NPHC_int_single_k(k + k_fine_list(ikfine,:), Chi_xyyy_k, Chi_yxxx_k)

            Chi_xyyy_tensor_mpi = Chi_xyyy_tensor_mpi + Chi_xyyy_k/dble(knv3_fine)
            Chi_yxxx_tensor_mpi = Chi_yxxx_tensor_mpi + Chi_yxxx_k/dble(knv3_fine)
        enddo
        
    enddo ! ik


#if defined (MPI)
    call mpi_barrier(mpi_cmw,ierr)
    call mpi_reduce(Chi_xyyy_tensor_mpi, Chi_xyyy_tensor, size(Chi_xyyy_tensor), mpi_dp,mpi_sum, 0, mpi_cmw,ierr)
    call mpi_reduce(Chi_yxxx_tensor_mpi, Chi_yxxx_tensor, size(Chi_yxxx_tensor), mpi_dp,mpi_sum, 0, mpi_cmw,ierr)
#endif

    Chi_xyyy_tensor = Chi_xyyy_tensor * NPHC_int_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume
    Chi_yxxx_tensor = Chi_yxxx_tensor * NPHC_int_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume

    if (Nk3==1) then
        Chi_xyyy_tensor = Chi_xyyy_tensor* Origin_cell%Ruc(3)/Ang2Bohr
        Chi_yxxx_tensor = Chi_yxxx_tensor* Origin_cell%Ruc(3)/Ang2Bohr
    endif

    outfileindex= outfileindex+ 1
    if (cpuid.eq.0) then
        do ieta=1, Eta_number
            write(Eta_name, '(f12.2)') Eta_array(ieta)*1000d0/eV2Hartree

            if (include_m_spin) then
                write(ahcfilename, '(7a)')'sigma_NPHC_int_S_eta', trim(adjustl(Eta_name)), 'meV.dat'
                open(unit=outfileindex, file=ahcfilename)
                write(outfileindex, '("#",a)')' Intrinsic nonlinear planar hall effect, in unit of A*V^-2*T^-1 for 3D cases, Ang*A*V^-2*T^-1 for 2D cases'
                write(outfileindex, '("#",a)')' Please refer to the Sec. III of the supplementary materials of 10.1103/PhysRevLett.130.126303, for the definition of term I and term II of the INPHE conductivities'

                write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', '\sigma_xyyy_I', '\sigma_xyyy_II', '\sigma_xyyy_tot', '\sigma_yxxx_I', '\sigma_yxxx_II','\sigma_yxxx_tot'
                do ie=1, OmegaNum
                    write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, &
                        Chi_xyyy_tensor(ie,ieta,1,1), Chi_xyyy_tensor(ie,ieta,1,2), Chi_xyyy_tensor(ie,ieta,1,1) + Chi_xyyy_tensor(ie,ieta,1,2),&
                        Chi_yxxx_tensor(ie,ieta,1,1), Chi_yxxx_tensor(ie,ieta,1,2), Chi_yxxx_tensor(ie,ieta,1,1) + Chi_yxxx_tensor(ie,ieta,1,2)
                enddo
                close(outfileindex)
            endif

            if (include_m_orb ) then
                write(ahcfilename, '(7a)')'sigma_NPHC_int_L_eta', trim(adjustl(Eta_name)), 'meV.dat'
                open(unit=outfileindex, file=ahcfilename)
                write(outfileindex, '("#",a)')' Intrinsic nonlinear planar hall effect, in unit of A*V^-2*T^-1 for 3D case, Ang*A*V^-2*T^-1 for 2D cases'
                write(outfileindex, '("#",a)')' Please refer to the Sec. III of the supplementary materials of 10.1103/PhysRevLett.130.126303, for the definition of term I and term II of the INPHE conductivities'

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


subroutine sigma_NPHC_int ! dynamical mpi version
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
            call sigma_NPHC_int_single_k(k, Chi_xyyy_k, Chi_yxxx_k)

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
                call sigma_NPHC_int_single_k(k + k_fine_list(ikfine,:), Chi_xyyy_k, Chi_yxxx_k)
    
                Chi_xyyy_tensor_mpi = Chi_xyyy_tensor_mpi + Chi_xyyy_k/dble(knv3_fine)
                Chi_yxxx_tensor_mpi = Chi_yxxx_tensor_mpi + Chi_yxxx_k/dble(knv3_fine)
            enddo
        enddo
    endif

    call mpi_reduce(Chi_xyyy_tensor_mpi, Chi_xyyy_tensor, size(Chi_xyyy_tensor), mpi_dp,mpi_sum, 0, mpi_cmw,ierr)
    call mpi_reduce(Chi_yxxx_tensor_mpi, Chi_yxxx_tensor, size(Chi_yxxx_tensor), mpi_dp,mpi_sum, 0, mpi_cmw,ierr)
#endif

    if (cpuid==0) then
        Chi_xyyy_tensor = Chi_xyyy_tensor * NPHC_int_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume
        Chi_yxxx_tensor = Chi_yxxx_tensor * NPHC_int_unit_factor /dble(knv3)/Origin_cell%CellVolume*kCubeVolume/Origin_cell%ReciprocalCellVolume

        if (Nk3==1) then
            Chi_xyyy_tensor = Chi_xyyy_tensor* Origin_cell%Ruc(3)/Ang2Bohr
            Chi_yxxx_tensor = Chi_yxxx_tensor* Origin_cell%Ruc(3)/Ang2Bohr
        endif

        outfileindex= outfileindex+ 1
        do ieta=1, Eta_number
            write(Eta_name, '(f12.2)') Eta_array(ieta)*1000d0/eV2Hartree

            if (include_m_spin) then
                write(ahcfilename, '(7a)')'sigma_NPHC_int_S_eta', trim(adjustl(Eta_name)), 'meV.dat'
                open(unit=outfileindex, file=ahcfilename)
                write(outfileindex, '("#",a)')' Intrinsic nonlinear planar hall effect, in unit of A*V^-2*T^-1 for 3D case, Ang*A*V^-2*T^-1 for 2D cases'
                write(outfileindex, '("#",a)')' Please refer to the Sec. III of the supplementary materials of 10.1103/PhysRevLett.130.126303, for the definition of term I and term II of the INPHE conductivities'

                write(outfileindex, '("#",a13, 20a16)')' Energy (eV)', '\sigma_xyyy_I', '\sigma_xyyy_II', '\sigma_xyyy_tot', '\sigma_yxxx_I', '\sigma_yxxx_II','\sigma_yxxx_tot'
                do ie=1, OmegaNum
                    write(outfileindex, '(200E16.8)')energy(ie)/eV2Hartree, &
                        Chi_xyyy_tensor(ie,ieta,1,1), Chi_xyyy_tensor(ie,ieta,1,2), Chi_xyyy_tensor(ie,ieta,1,1) + Chi_xyyy_tensor(ie,ieta,1,2),&
                        Chi_yxxx_tensor(ie,ieta,1,1), Chi_yxxx_tensor(ie,ieta,1,2), Chi_yxxx_tensor(ie,ieta,1,1) + Chi_yxxx_tensor(ie,ieta,1,2)
                enddo
                close(outfileindex)
            endif

            if (include_m_orb ) then
                write(ahcfilename, '(7a)')'sigma_NPHC_int_L_eta', trim(adjustl(Eta_name)), 'meV.dat'
                open(unit=outfileindex, file=ahcfilename)
                write(outfileindex, '("#",a)')' Intrinsic nonlinear planar hall effect, in unit of A*V^-2*T^-1 for 3D case, Ang*A*V^-2*T^-1 for 2D cases'
                write(outfileindex, '("#",a)')' Please refer to the Sec. III of the supplementary materials of 10.1103/PhysRevLett.130.126303, for the definition of term I and term II of the INPHE conductivities'

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


subroutine sigma_NPHC_tau2 ! dynamical mpi version
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
