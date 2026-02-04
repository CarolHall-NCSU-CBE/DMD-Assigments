module global
        implicit none
        save
        integer, parameter  :: n=108
        real, parameter  :: pi = 3.14159265359
        !real, parameter  :: box = 1.0, rlower = 0.0, rupper = box/2.0 
        real, parameter  :: timbig = 1.0e10                         ! Declare and intiialize max time for minimum collision time
        real, parameter :: tol = 1.0E-4					        	! tolerance for overlaps, no more than 0.0001
        real, parameter :: maxbin = 600     !radial distribution: number of histogram bins
        real, parameter :: delr = 0.001     !radial distribution: bin width, the delta r (size of the RDF shell)
        
        real :: series(n), rx(n), ry(n), rz(n), vx(n), vy(n), vz(n)
        real :: coltim(n+1)                                         !Declare collision time array
        integer :: partner(n)								        !Declare partner array
        real :: sigma									        	!Declare sigma = sphere diameter
        real :: total_energy                                        !total kinetic energy
        real :: energy                                              !
        !real :: enkt                                               ! total energy/ kT
        real :: time                                                ! accumulated system time
        real :: w                                                   ! virial from one collision
        real :: acw = 0                                             ! cumulative virial used for calculation Z=RV/NkT
        real :: Z                                                   !*** compressibility factor ***
        real :: density                                             !reduced density
        real :: destemp                                             !reduced temperature
        
        integer :: ncoll                                            !total collision number
        integer :: current_i, current_j                             !current colliding particles i and j
        integer :: counter                                          ! collision counter
        integer :: bin                                              ! Index for the gr histogram bins 
        integer :: hist(maxbin)                                     ! histogram for radial distribution
        integer :: grcounter                                        ! radial distribution counter
        real :: delv                                        ! radial distribution differential volume
        real :: gr_density                                  ! density at r
        real :: gr(maxbin)                                  ! radial distribution function with 500 bins
        real :: r_list(maxbin)                              ! used to store all radia data
            
        logical :: overlap
        
end module global
    
program hardsphere
    use global
    implicit none
    integer, parameter :: xyz_file_code = 11
    integer :: i, j            ! loop value
    real :: r1, r2              ! small and big radius for radial distribution differential volume
    character(len=8) :: date 
    character(len=20) :: world_time
    character(len=6) :: world_time_trim
    character(10) :: temp_s, n_s, density_s  ! open coordinate file 12.23
    character(40) :: gr_filename
    REAL :: e, en, enkt, pvnkt1, temp, tbc   ! energies
    
    !*********** Enter parameters **********************************************************
    write (*,*) "Enter reduced temperature kT/epsilon (suggested 1.5)"
    read (*,*) destemp                                      ! Input temperature
    write (*,*) "Enter reduced density (n/v)*sigma**3"	
    read (*,*) density										! Input density
    write (*,*) "Enter number of collisions required"		
    read (*,*) ncoll										! Input total collision number
    
    !*********** file names ********
    call date_and_time(date, world_time)    ! Get the date and time only for record output files
    world_time_trim = world_time(1:6)
    write(temp_s, '(F5.2)') destemp
    write(n_s, '(I4)') n
    write(density_s, '(F5.2)') density
    gr_filename = 'HS_gr_p'//trim(adjustl(trim(density_s)))//'_T'//trim(adjustl(trim(temp_s)))//'_n'//trim(adjustl(trim(n_s)))//'.csv'

    !initialize parameters **********************************************************
    time = 0          ! System time
    counter = 0       ! Collision number counter      
    hist = 0          ! initialize histogram for countering particles in certain bin (distance)
    grcounter = 0     ! Number of radial distribution recorded
    
    call coordinate_and_vel()
    sigma = (density/REAL(n))**(1.0/3.0)                    ! Calculate the sphere diameter in reduce unit
       
    total_energy = 0
    do i = 1, n
        total_energy = total_energy + vx(i)**2 + vy(i)**2 + vz(i)**2
    end do
    total_energy = 0.5*total_energy                         ! Calculate the total energy
    
    !en = e/REAL(n)													! Calculate the energy that will be used for the temperature calculation
    temp = 2.0/3.0 * total_energy/real(n)							! Calculate the temperature
    write (*, "(A,F10.5)") "Calculated Reduced Temperature:", temp  ! Recalculate the initial temperature, which is different than input
    
    call check_position()                                   ! Check initial configuration overlaps

    !start collisions **********************************************************
    call calculate_all_pairs()                              ! Do an initial collision
    call find_tmin_and_update()
    
    !open(unit=xyz_file_code, file='HS_xyz'//date//'_'//world_time_trim// '.xyz', status='replace', action='write')
    !call write_xyz_file(xyz_file_code)
    
    do while (counter < ncoll)                              ! Main loop
        
        call update_pairs()
        call find_tmin_and_update()
        
        
        if (mod(counter, 20) == 0 ) then    !record radial distribution every 20 collisions
            call update_radial_bin()
            grcounter = grcounter + 1
        end if
    
        !call write_xyz_file(xyz_file_code)
    end do
    !end collisions **********************************************************
    
    call check_position()
    
    Z = 1.0 + (acw / (real(n) * 3.0 * time * temp))                ! calculate compressibility factor
    
    
    total_energy = 0
    do i = 1, n
        total_energy = total_energy + vx(i)**2 + vy(i)**2 + vz(i)**2
    end do
    total_energy = 0.5*total_energy                                 !calcualte total energy again
    enkt = total_energy/(REAL(n) * temp)                            !Calculate E/kT

    
    write(*,"(A,F13.8)") "E/nkT = ", enkt
    write(*,"(A,F13.8)") "Z = PV/NkT = ", Z
    write(*,"(A,E15.7E2)") "Passed time =", time
    
    !calculation radial distribution **********************************************************
    do i = 1, maxbin
        r1 = (i - 1) * delr 
        r2 = i * delr
        
        delv = 4.0 * pi / 3.0 * (r2**3 - r1**3)        
        gr_density = real(hist(i)) / (real(grcounter) * real(n) * delv)           ! the g(r) is defined as (density at r)/(uniform density)
        gr(i) = gr_density / (density/sigma**3)                 ! reduce density/sigma**3 == real density = n/V = n here(V= 1^3) 
        r_list(i) = r1 + 0.5 * delr
        
    end do
    
    !export radial distribution data **********************************************************
    open(unit=10, file=gr_filename, status='replace', action='write')   ! output radial distribution data
    
    write(10,"(4A)") "Date,", date, ", Time,", world_time_trim
    write(10,"(A,I10)") "total particles = ,", n
    write(10,"(A,F10.2)") "reduced T = ,", temp
    write(10,"(A,F10.5)") "reduced density = ,", density
    write(10,"(A,I10)") "total collisions = ,", ncoll
    write(10,"(A,F13.8)") "E/nkT = ,", enkt
    write(10,"(A,F13.8)") "Z = PV/NkT = ,", Z
    write(10,"(A,E15.7E2)") "Passed time = ,", time
    
    do i = 1, maxbin
        write(10, *) i, ",", r_list(i), ",", gr(i), ","
        !bincounter = bincounter + hist(i)
    end do
    
end program
    
    

subroutine coordinate_and_vel()
    use global
    implicit none

    integer nc, nc2, i, j, k, m      !number of cell on one side
    real cell_len, cell_len2    
    real*8 :: sumx, sumy, sumz
    real :: rtemp, gauss, dummy
    real*8 :: vxsum, vysum, vzsum
	
	call init_random_seed() !initialize randam seed, only need call once 
    
! Producing number of 4*(3^3) spheres ***********************************
    nc = (n/4)**(1.0/3.0) 
    nc2 = nc**2
    ! Generate FCC lattice
	cell_len = 1.0/real(nc) ! lattice length
	cell_len2 = cell_len/2.0
 
	do i = 1, nc
		do j = 1, nc
			do k = 1, nc
			
				rx(((i-1)*nc2 + (j-1)*nc + k-1)*4 + 1) = cell_len*real(k-1)
				ry(((i-1)*nc2 + (j-1)*nc + k-1)*4 + 1) = cell_len*real(j-1)
				rz(((i-1)*nc2 + (j-1)*nc + k-1)*4 + 1) = cell_len*real(i-1)
				
				rx(((i-1)*nc2 + (j-1)*nc + k-1)*4 + 2) = cell_len*real(k-1) + cell_len2
				ry(((i-1)*nc2 + (j-1)*nc + k-1)*4 + 2) = cell_len*real(j-1) + cell_len2
				rz(((i-1)*nc2 + (j-1)*nc + k-1)*4 + 2) = cell_len*real(i-1)
			
				rx(((i-1)*nc2 + (j-1)*nc + k-1)*4 + 3) = cell_len*real(k-1) + cell_len2
				ry(((i-1)*nc2 + (j-1)*nc + k-1)*4 + 3) = cell_len*real(j-1)
				rz(((i-1)*nc2 + (j-1)*nc + k-1)*4 + 3) = cell_len*real(i-1) + cell_len2
				
				rx(((i-1)*nc2 + (j-1)*nc + k-1)*4 + 4) = cell_len*real(k-1)
				ry(((i-1)*nc2 + (j-1)*nc + k-1)*4 + 4) = cell_len*real(j-1) + cell_len2
				rz(((i-1)*nc2 + (j-1)*nc + k-1)*4 + 4) = cell_len*real(i-1) + cell_len2
				
			end do
		end do
    end do

! shift coordinate to the center
    do i = 1, n
        rx(i) = rx(i) - 0.5
        ry(i) = ry(i) - 0.5
        rz(i) = rz(i) - 0.5
    end do
! generate initial velocity***********************************
    dummy = 0 ! dummy vaiables for picking number in gaussian distribution
    rtemp = sqrt (destemp)    

    do i = 1, n   
           vx(i) = rtemp * gauss (dummy)
           vy(i) = rtemp * gauss (dummy)
           vz(i) = rtemp * gauss (dummy)

        end do
        
!       ** remove net momentum from system **
        sumx = 0.0
        sumy = 0.0
        sumz = 0.0
        vxsum = 0.0
        vysum = 0.0
        vzsum = 0.0
        
        do i = 1, n
           sumx = sumx + vx(i)
           sumy = sumy + vy(i)
           sumz = sumz + vz(i)
        end do
        
        sumx = sumx / real(n)
        sumy = sumy / real(n)
        sumz = sumz / real(n)
        
        do i = 1, n-1
           vx(i) = vx(i) - sumx
           vy(i) = vy(i) - sumy
           vz(i) = vz(i) - sumz
           vxsum = vxsum + vx(i)
           vysum = vysum + vy(i)
           vzsum = vzsum + vz(i)
           
        end do
        
        vx(n) = -vxsum
        vy(n) = -vysum
        vz(n) = -vzsum
        
        
        sumx = 0.0
        sumy = 0.0
        sumz = 0.0
        do i = 1, n
            sumx = sumx + vx(i)
            sumy = sumy + vy(i)
            sumz = sumz + vz(i)
        end do
        
        write(*,*) sumx, sumy, sumz
            
        
end subroutine coordinate_and_vel   
    
! Gaussian distribution***********************************    
real function gauss (dummy)
        real, parameter :: a1 = 3.949846138
        real, parameter :: a3 = 0.252408784
        real, parameter :: a5 = 0.076542912
        real, parameter :: a7 = 0.008355968
        real, parameter :: a9 = 0.029899776
        real    ::  rsum, r, r2, dummy, ran
        integer ::  i
        
        rsum = 0.0

        do i = 1, 12            
            call random_number(ran)
            !write(*,*) ran           
            rsum = rsum + ran
        end do
        
        r = (rsum - 6.0) / 4.0
        r2 = r * r
        
        gauss = (((( a9 * r2 + a7 ) * r2 + a5 ) * r2 + a3 ) * r2 + a1 ) * r
        return
    end function gauss
    
! Initialize Random seeds*********************************** 
subroutine init_random_seed()
    integer :: i, k, clock
    integer, dimension(:), allocatable :: seed
   
    call random_seed(size = k)
    allocate(seed(k))
   
    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, k) /)
    call random_seed(put = seed)
   
    deallocate(seed)
end subroutine init_random_seed
    
! Check overlaps ************************************************************
subroutine check_position() 
    use global
    implicit none
    
    real :: sigma2
    real :: rxi, ryi, rzi, rxj, ryj, rzj, rxij, ryij, rzij, rij2, rij
    
    integer :: i, j                                                    ! loop values
    
    sigma2 = sigma**2												   !Square dat sigma
    overlap= .false.
    
    do i=1, n-1   
        do j = i+1, n

            rxij = rx(i) - rx(j)
            ryij = ry(i) - ry(j)
            rzij = rz(i) - rz(j)
            rxij = rxij - anint(rxij)
            ryij = ryij - anint(ryij)
            rzij = rzij - anint(rzij)
            
            rij2 = rxij**2 + ryij**2 + rzij**2
            
            if (rij2 < sigma2) then ! check if distance between a pair is smaller than sphere diameter
                write (*,*) "rij < diameter spheres", i, j
                
                if ( (1-(rij2**0.5/sigma) ) > tol ) then
                    overlap = .true.
                    write (*,*) "overlap spheres", i, j
                end if
            end if
        end do
    end do
    
    if (overlap) then
        write (*,'(A,i8,A)') "The system has overlaps at ", counter, " collision." // CHAR(10)
        stop
    else if (overlap == .false.) then
        write (*, '(A,i8,A)') "The system has no overlap at ", counter, " collision."// CHAR(10)
    end if
end subroutine check_position                    


subroutine calculate_all_pairs()
    use global
    implicit none

    integer :: i                            ! loop values
    
    do i = 1, n                             ! assign maximum collision time
        coltim(i) = timbig
    end do
    
    do i = 1, n
        call uplist(i)
    end do
    
end subroutine calculate_all_pairs    
    
! find the minimal tij and colliding pair ij, update all collision time, position, and velocity
subroutine find_tmin_and_update()
    use global 
    implicit none
    
    real :: tij
    integer :: k
    integer :: i, j
    real :: rxij, ryij, rzij, vxij, vyij, vzij
    real :: delvx, delvy, delvz
    real :: bij
    real :: factor
    
    tij = timbig        ! set tij to a large value
    
    do k = 1, n         ! find the minimal tij and colliding pair i j
        
        if (coltim(k) < tij) then
            tij = coltim(k)
            i = k
        end if
        
    end do

    if (tij < 0) write(*,*) "tij < 0 at step:", counter, "bead:", i
    
    j = partner (i)
    current_i = i
    current_j = j
    
    time = time + tij
    
    do k = 1, n
        
        coltim(k) = coltim(k) - tij     !每个粒子的碰撞时间减小tij，并更新所有粒子tij后的坐标
        
        rx(k) = rx(k) + vx(k) * tij
        ry(k) = ry(k) + vy(k) * tij
        rz(k) = rz(k) + vz(k) * tij
        rx(k) = rx(k) - anint(rx(k))    
        ry(k) = ry(k) - anint(ry(k))
        rz(k) = rz(k) - anint(rz(k))
        
    end do
    
    rxij = rx(i) - rx(j)
    ryij = ry(i) - ry(j)
    rzij = rz(i) - rz(j)
    rxij = rxij - anint(rxij)
    ryij = ryij - anint(ryij)
    rzij = rzij - anint(rzij)
            
    vxij = vx(i) - vx(j)
    vyij = vy(i) - vy(j)
    vzij = vz(i) - vz(j)
    
    bij = rxij * vxij + ryij * vyij + rzij * vzij
    factor = -bij/sigma**2
    
    delvx = factor * rxij
    delvy = factor * ryij
    delvz = factor * rzij
    
    vx(i) = vx(i) + delvx
    vy(i) = vy(i) + delvy
    vz(i) = vz(i) + delvz
    
    vx(j) = vx(j) - delvx
    vy(j) = vy(j) - delvy
    vz(j) = vz(j) - delvz
    
    w = delvx * rxij + delvy * ryij + delvz * rzij
    
    acw = acw + w

    !write (*,"(A20, I4, I4, A10, E15.7E2)") "collision pairs:", current_i, current_j, "time:", tij
    !write (*,"(A30, E15.7E2, A20, E15.7E2)") "virial for this collision:", w, "cumulative virial:",acw 
    !
    !do i = 1, n
    !    write (*,*) i, rx(i), ry(i), rz(i), vx(i), vy(i), vz(i)
    !end do
    
    counter = counter + 1 ! increase collision counter by 1
    !write (*,*) counter
    !write (*,*)
end subroutine find_tmin_and_update

! only update the pairs that are affected by the last collision pairs
subroutine update_pairs()
    use global
    implicit none

    integer :: i, j, ii, jj, k
      
    ii = current_i
    jj = current_j
    
    do k = 1, n
        if ((k == ii) .or. (partner(k) == ii) .or. (k == jj) .or. (partner(k) == jj)) then
            call uplist(k)
        end if
    end do
    
    call dnlist(ii)
    call dnlist(jj)
            
    if (mod(counter, 100) == 0 ) then    ! periodically call whole uplist
        do i = 1, n
            call uplist(i)
        end do
    end if

end subroutine update_pairs

! calculate radial distribution total particles in every bin*********************
subroutine update_radial_bin()
    use global
    implicit none
    
    integer :: i, j
    
    real :: rxij, ryij, rzij, rij2, rij
    real :: sigma2
    !integer :: bincounter
    !bincounter = 0
    
    sigma2 = sigma**2
    
    do i = 1, n-1
        do j = i+1, n
            rxij = rx(i) - rx(j)
            ryij = ry(i) - ry(j)
            rzij = rz(i) - rz(j)
            rxij = rxij - anint(rxij)
            ryij = ryij - anint(ryij)
            rzij = rzij - anint(rzij)
    
            rij2 = rxij**2 + ryij**2 + rzij**2          ! Calculate distance between all pairs
            rij = sqrt(rij2)
            bin = int(rij/delr) + 1                     ! Determine these spheres belong to which bin

            if(bin.LE.maxbin) then						! If the current bin number is less than the maximum number of bins allowed...
	            hist(bin) = hist(bin) + 2				! Update the histogram by 2 because a pair of spheres
            end if 		

        end do
    end do
end subroutine update_radial_bin

subroutine uplist(i)
    use global
    implicit none

    integer :: i    ! enternal variable
    integer :: j    ! loop values
    real :: rxi, ryi, rzi
    real :: vxi, vyi, vzi
    real :: rxij, ryij, rzij, rij2, rij
    real :: vxij, vyij, vzij, vij2, vij
    real :: tij
    real :: bij      ! vector r_ij * vector v_ij
    real :: discr
    real :: sigma2
    
    if (i == n) return
    sigma2 = sigma**2
    ! assign maximum collision time
    coltim(i) = timbig
    
    rxi = rx(i)
    ryi = ry(i)
    rzi = rz(i)
    vxi = vx(i)
    vyi = vy(i)
    vzi = vz(i)
    
    do j = i + 1, n
        rxij = rxi - rx(j)
        ryij = ryi - ry(j)
        rzij = rzi - rz(j)
        rxij = rxij - anint(rxij)
        ryij = ryij - anint(ryij)
        rzij = rzij - anint(rzij)
            
        vxij = vxi - vx(j)
        vyij = vyi - vy(j)
        vzij = vzi - vz(j)
        
        bij = rxij * vxij + ryij * vyij + rzij * vzij
            
        if (bij < 0.0) then
            rij2 = rxij**2 + ryij**2 + rzij**2
            vij2 = vxij**2 + vyij**2 + vzij**2
            discr = bij**2 - vij2*(rij2 - sigma2) ! qudratic determinant, without coefficient 4
            
            if (discr > 0) then     ! calculate collision t_ij，compare with a current smallest collision time, if it is smaller, update
                tij = (-bij - sqrt(discr))/vij2
                    
                if (tij < coltim(i)) then
                    coltim(i) = tij
                    partner(i) = j 
                end if
                    
                !if (tij < coltim(j)) then
                !    coltim(j)  = tij
                !    partner(j) = i
                !end if 
                    
            end if
        end if
    end do
    
    return
end subroutine uplist

    
subroutine dnlist(j)
    use global
    implicit none
	! calculate collision time between j 和 i = 1~ (j-1) 
    integer :: i    
    integer :: j    ! external variable
    real :: rxj, ryj, rzj
    real :: vxj, vyj, vzj
    real :: rxij, ryij, rzij, rij2, rij
    real :: vxij, vyij, vzij, vij2, vij
    real :: tij
    real :: bij                     ! vector r_ij * vector v_ij
    real :: discr
    real :: sigma2
    
    if (j == 1) return
    sigma2 = sigma**2
                                    ! assign maximum collision time
    
    rxj = rx(j)
    ryj = ry(j)
    rzj = rz(j)
    vxj = vx(j)
    vyj = vy(j)
    vzj = vz(j)
    
    do i = 1, j-1
        rxij = rx(i) - rxj
        ryij = ry(i) - ryj
        rzij = rz(i) - rzj
        rxij = rxij - anint(rxij)
        ryij = ryij - anint(ryij)
        rzij = rzij - anint(rzij)
            
        vxij = vx(i) - vxj
        vyij = vy(i) - vyj
        vzij = vz(i) - vzj
        
        bij = rxij * vxij + ryij * vyij + rzij * vzij
            
        if (bij < 0.0) then
            rij2 = rxij**2 + ryij**2 + rzij**2
            vij2 = vxij**2 + vyij**2 + vzij**2
            discr = bij**2 - vij2*(rij2 - sigma2) 
            
            if (discr > 0) then
                tij = (-bij - sqrt(discr))/vij2
                    
                if (tij < coltim(i)) then
                    coltim(i) = tij
                    partner(i) = j 
                end if
                    
                !if (tij < coltim(j)) then
                !    coltim(j)  = tij
                !    partner(j) = i
                !end if 
                    
            end if
        end if
    end do
    
    return
end subroutine dnlist
    
subroutine write_xyz_file(file_code)
    use global
    implicit none

    integer :: i, file_code

    write(file_code, *) n
    write(file_code, *) " "

    do i = 1, n
        write(file_code, *) "Ar", rx(i), ry(i), rz(i)
    end do
    
return
end subroutine write_xyz_file