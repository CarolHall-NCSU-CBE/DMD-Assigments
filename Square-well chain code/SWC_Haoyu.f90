!  SWC.f90 
!
!  FUNCTIONS:
!  SWC - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: SWC
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

module global
        implicit none
        save
        
        integer, parameter  :: numchn = 64  ! ***Number of chains
        integer, parameter :: chainlen = 4  ! ***chain length (contains this number of spheres)
        integer, parameter :: n = 256       ! total number of spheres = chain_lenth * total chains.
        
        real(8), parameter  :: pi = 3.14159265359
        real, parameter :: lambda = 1.5     ! square-well diameter multiple
        real, parameter :: eps = 1.0        ! square-well potential parameter

        real, parameter :: timbig = 1.0e10  !Declare and intiialize max time for minimum collision time
        real, parameter :: tol = 1.0E-4                         
        real, parameter :: delbond = 0.01   ! covalent bond width
        real, parameter :: maxbin = 1000    ! radial distribution: number of bins
        real, parameter :: delr = 0.0005    ! radial distribution: bin width
        real, parameter :: ghostcoeff = 0.01
                
        real(8) :: series(n), rx(n), ry(n), rz(n), vx(n), vy(n), vz(n)
        real :: coltim(n+1)                                         ! Declare collision time array
        
        real(8) :: sigma1, sigma2                                   ! *** sigma1 = hardsphere diameter, sigma2 = square-well diameter
        real(8) :: d1, d2                                           ! *** 2024/1/5 max bond length, min bond length
        real(8) :: time, w, Z, Z2, density                          ! past time, virial per collision, compressibility factor,reduced density, 
        real :: destemp                                             ! reduced temperature
        real :: acw = 0                                             ! cumulative virial used for calculation Z=RV/NkT
        real :: acw2 = 0                                            ! 2024/1/22 don't considering the intrachain forces
        
        integer :: partner(n), coltype(n)                           ! partner array, collision type array
        integer :: ncoll, current_i, current_j, counter             ! total collision number, current pair ij, colliison number counter
        integer :: intra_counter                                    ! *** 2024/1/5 intrachain collision counter
        integer :: bonded_spheres(2, n)                             ! *** 2024/1/5 2D array for recording the covalent bond parter
        integer :: ngcoll_counter                                   ! counter records times when tij < 0
        integer :: bin, grcounter                                   ! bin index, number of gr recorded
        integer :: hist(maxbin)                                     ! histogram for radial distribution

        real :: delv                                          ! radial distribution differential volume
        real :: gr_density                                  ! density at r
        real :: gr(maxbin), r_list(maxbin)      ! *gr function*, * gr axis* 
            
        logical :: overlap
        
end module global
    
program swchain
    use global
    implicit none
    integer, parameter :: pdb_file_code = 30
    integer, parameter :: xyz_file_code = 11
    integer, parameter :: rv_file_code = 12
    integer, parameter :: gr_file_code = 13
    integer, parameter :: temp_file_code = 15
    real, parameter :: burninT0 = 500.0
    real, parameter :: burninT1 = 100.0
    real, parameter :: burninT2 = 10.0

    integer, parameter :: burnin1 = 3000000
    integer, parameter :: burnin2 = 6000000
    integer, parameter :: burnin3 = 9000000
    integer, parameter :: burnin4 = 15000000
    
    real, parameter  :: inittime = 0.001
    
    integer :: i, j, k, status      ! loop value
    integer :: numghost, ghost 
    
    real(8) :: tij, rij2, sig2sq, bij
    real(8) :: rxij, ryij, rzij, vxij, vyij, vzij
    real :: ran1, ran2, ran3           ! random number
    real :: ke, pe               ! kinetic energy, potential energy
    real :: ghosttime, burnin_time
    real :: r1, r2               ! small and big radius for radial distribution differential volume
    real :: e, en, enkt1, enkt2, enkt3, pvnkt1, temp, tbc
    character(10) :: temp_s, n_s, density_s, chainlen_s  ! open coordinate file 12.23
    character(40) :: rv_filename, pdb_filename, temp_filename, gr_filename ! open coordinate file 12.23
    character(len=8) :: date1, date2 
    character(len=20) :: world_time1, world_time2
    character(len=6) :: world_time1_trim, world_time2_trim
    
    call init_random_seed()
    
    !********enter parameters **********************************************************
    write (*,"(A)") "Square-Well Chains simulation. For running different chain lengths, please modify the code"
    write (*,"(A)") "n = chain-length * total chains"
    write (*,"(A)") "Enter reduced temperature kT/epsilon (suggested 1.5)"
    read (*,*) destemp                                      !Input temperature
    write (*,"(A)") "Enter reduced density (n/v)*sigma**3"	
    read (*,*) density										!Input density
    write (*,"(A)") "Enter number of collisions required"		
    read (*,*) ncoll										!Input total collision number
    
    !********file name strings **********************************************************
    call date_and_time(date1, world_time1)
    world_time1_trim = world_time1(1:6)
    write(temp_s, '(F5.2)') destemp
    write(n_s, '(I4)') n
    write(density_s, '(F5.2)') density
    write(chainlen_s, '(I3)') chainlen
    
    !********initialize parameters **********************************************************
    sigma1 = (density/REAL(n))**(1.0/3.0)                    
    sigma2 = lambda * sigma1
    d1 = sigma1 * (1 - delbond) ! minimal bond length
    d2 = sigma1 * (1 + delbond) ! maximum bond length
    bonded_spheres = -1         ! Initialize bond partner
    
    numghost = 0
    time = 0.0
    counter = 0
    intra_counter = 0
    hist = 0                     ! Initialize histogram for countering particles in certain bin (distance)
    grcounter = 0
    ngcoll_counter = 0
    call random_number(ran1)
    call random_number(ran2)
    ghosttime = ghostcoeff * abs(sqrt(-2 * alog(ran1)) * cos(2 * pi * ran2))
    coltim (n+1) = ghosttime 

    
    
    !*********CREATE coordinates and initial velocities **********************************************************
    call create_coordinate()
    call create_vel(burninT0) ! assign T* = 500 to burn it
    
    !*********READ coordinates and initial velocities **********************************************************
    !write(temp_string , '(F4.2)') destemp
    !write(n_string, '(I3)') n
    !rv_filename = "rv/rv_"//trim(temp_string)//"_"//trim(n_string)//".txt"
    !open(rv_file_code, file=rv_filename, status='old', action='read', iostat=status)
    !
    !if (status /= 0) then
    !    write(*,*) 'Error opening the file!'
    !    stop
    !end if
    !
    !read(rv_file_code,*)
    !do i = 1, n
    !    read(rv_file_code, '(F20.15,F20.15,F20.15,F20.15,F20.15,F20.15)') rx(i), ry(i), rz(i), vx(i), vy(i), vz(i)
    !    !write(*,'(F10.5,F10.5,F10.5,F20.15,F20.15,F20.15)') rx(i), ry(i), rz(i), vx(i), vy(i), vz(i)
    !end do
    !close(rv_file_code)
    
    !********Open output files **********************************************************
    pdb_filename = 'SWchain_pdb_p'//trim(adjustl(trim(density_s)))//'_T'//trim(adjustl(trim(temp_s)))//'_n'//trim(adjustl(trim(n_s)))//'_'//trim(adjustl(trim(chainlen_s)))//'mer.pdb'
    temp_filename = 'SWchain_temp_p'//trim(adjustl(trim(density_s)))//'_T'//trim(adjustl(trim(temp_s)))//'_n'//trim(adjustl(trim(n_s)))//'_'//trim(adjustl(trim(chainlen_s)))//'mer.csv'
    !open(unit=pdb_file_code, file=pdb_filename, status='replace', action='write')
    open(unit=temp_file_code, file=temp_filename, status='replace', action='write')
    
    !*******recalculate temperature **********************************************************
    call kecal(ke)              !Calculate the kinetic energy
    call pecal(pe)              !Calculate the potential energy
    !en = e/REAL(n)
    temp = 2.0/3.0 * ke/real(n)	!Calculate the temperature
    													
    									
    write (*, "(A,F10.5)") "Calculated Reduced Temperature:", temp
    
    call check_position()

    !***************do initial collision **********************************************************

    do i = 1, n                          
        coltim(i) = timbig              ! assign maximum collision time
        coltype(i) = 5
        partner(i) = n
    end do
    do i = 1, n - 1
        call uplist(i)
    end do
    counter = counter + 1
    ! * write in pdb file *
    
    !call write_pdb_file(pdb_file_code, pdb_filename)
    
    !***************************************
    !**************Main Loop****************
    !***************************************
    do while (counter < ncoll) ! main loop
        !**************Burnin time and Update full list every 100 time ****************
        ! Or annealing
        if (counter == burnin1) then
            call create_vel(burninT1) 
            do k = 1, n
                call uplist(k)
            end do
            write(*,"(A,F8.4)") "Phase 1 T* = 500, burnin finished at", time 
            
        else if (counter == burnin2) then
            call create_vel(burninT2) 
            do k = 1, n
                call uplist(k)
            end do
            write(*,"(A,F8.4)") "Phase 2 T* = 100, burnin finished at", time 
            
        else if (counter == burnin3) then
            call create_vel(destemp) 
            do k = 1, n
                call uplist(k)
            end do
            write(*,"(A,F8.4)") "Phase 3 T* = 10, burnin finished at", time
            
        else if (counter == burnin4) then
            acw = 0
            burnin_time = time
            write(*,"(A,F6.2,A,F8.4)") "Phase 4 T* =", destemp, ", burnin finished at", time
            write(*,"(A)") "Now entering equilibrium phase"
            
        else if (mod(counter, 200) == 0) then
            do k = 1, n
                call uplist(k)
            end do 
        
        end if
    
        !!***************** Write Temperature vs Coll ********************
        if (mod(counter, 600) == 0) then
            call pecal(pe)
            call kecal (ke)
            temp = 2.0 * ke / 3.0 / real (n)
            
            write(temp_file_code, "(I10,A,F15.10,A,F15.8,A,F15.8,A,F15.8)") counter, ",", time, ",", temp, ",", pe, ",", ke
        endif
    
        !**************Find t_min and negative collision times****************
        tij = timbig
        do k = 1, n + 1
            if (coltim(k) < tij) then
                tij = coltim(k)
                i = k
            end if
        end do
        
        if (tij < 0) then   
            ngcoll_counter = ngcoll_counter + 1
        end if
        !**************Update Partners****************
        if (i <= n) then
            j = partner(i)
        endif

        current_i = i
        current_j = j
        time = time + tij
        
        if (counter >= burnin4) then
            coltim(n+1) = coltim(n+1) - tij !after burnin , ghost collision time reduce by tij，update position for all particles
        end if

        !**************Move particles by t_min****************
        do k = 1, n
            coltim(k) = coltim(k) - tij
            rx(k) = rx(k) + vx(k) * tij
            ry(k) = ry(k) + vy(k) * tij
            rz(k) = rz(k) + vz(k) * tij
            rx(k) = rx(k) - anint(rx(k))    
            ry(k) = ry(k) - anint(ry(k))
            rz(k) = rz(k) - anint(rz(k))
        end do
        
        if (mod(counter, 100) == 0) then
            call check_position()
        end if
        
        ! ********* Update velocities ***********
        if (i <= n) then
            call bump()
        
        ! ********* Update list for affected particles ********
            do k = 1, n
                if ((k == i) .or. (partner(k) == i) .or. (k == j) .or. (partner(k) == j)) then
                    call uplist(k)
                end if
            end do
    
            call dnlist(i)
            call dnlist(j)
   
        else ! i == n + 1
            call random_number(ran1)
            call random_number(ran2)
            call random_number(ran3)
            ghost = int(ran3 * real(n)) + 1
            call ghostcoll(ghost)
            ghosttime = ghostcoeff * abs(sqrt(-2 * alog(ran1)) * cos(2 * pi * ran2))
            coltim (n+1) = ghosttime    !reset ghost time
            
            do k = 1, n
                 if ((k == ghost) .or. (partner(k) == ghost)) then
                    call uplist(k)
                 endif
            end do
            call dnlist (ghost)
            numghost = numghost + 1
        end if
        
        if ((mod(counter, 20) == 0 ) .AND. (counter >= burnin4)) then
        !if (mod(counter, 20) == 0 ) then
            call update_radial_bin()
            grcounter = grcounter + 1
        end if
        
        !call write_pdb_file(pdb_file_code, pdb_filename)
        counter = counter + 1
    end do
    !****************************************
    !***********end collisions **************
    !****************************************
    
    call check_position()
    
    Z = 1.0 + (acw / (real(n) * 3.0 * (time - burnin_time) * temp)) ! 计算压缩因子从平衡态的开始计算,现在是只计算一个粒子的Z
    Z = Z * chainlen ! 转换为分子的Z
    
    Z2 = 1.0 + (acw2 / (real(n) * 3.0 * (time - burnin_time) * temp))
    Z2 = Z * chainlen 
    call kecal(ke)
    call pecal (pe)
    enkt1 = ke/(REAL(n) * temp)              !Calculate E_kinetic/kT
    enkt2 = pe/(REAL(n) * temp)              !Calculate E_potential/kT
    enkt3 = enkt1 + enkt2


    write(*,"(A,F13.8)") "E_kinetic/nkT = ", enkt1
    write(*,"(A,F13.8)") "E_potential/nkT = ", enkt2
    write(*,"(A,F13.8)") "Z = PV/NkT = (considering intrachain virial)", Z
    write(*,"(A,F13.8)") "Z = PV/NkT = (on intrachain virial)", Z2
    write(*,"(A,E15.7E2)") "Passed time =", time
    write(*,"(A,I5)") "Number of ghost =", numghost
    
    !calculation radial distribution **********************************************************
    do i = 1, maxbin
        r1 = (i - 1) * delr 
        r2 = i * delr
        
        delv = 4.0 * pi / 3.0 * (r2**3 - r1**3)        
        gr_density = real(hist(i)) / (real(grcounter) * real(n) * delv)     ! the g(r) is defined as (density at r)/(uniform density)
        gr(i) = gr_density / (density/sigma1**3)                            ! reduce density/sigma**3 == real density = n/V = n here(V= 1^3) 
        r_list(i) = r1 + 0.5 * delr
        
    end do
    
    !export radial distribution data **********************************************************
    gr_filename = 'SWchain_gr_p'//trim(adjustl(trim(density_s)))//'_T'//trim(adjustl(trim(temp_s)))//'_n'//trim(adjustl(trim(n_s)))//'_'//trim(adjustl(trim(chainlen_s)))//'mer.csv'
    
    
    open(unit=gr_file_code, file=gr_filename, status='replace', action='write')   ! output radial distribution data
    call date_and_time(date2, world_time2)
    world_time2_trim = world_time2(1:6)
    
    write(gr_file_code,"(5A)") "Before simulation,", "Date,", date1, ", Time,", world_time1_trim
    write(gr_file_code,"(5A)") "After simulation,", "Date,", date2, ", Time,", world_time2_trim
    write(gr_file_code,"(A,I8)") "Number of chains = ,", numchn
    write(gr_file_code,"(A,I8)") "chain length = ,", chainlen
    write(gr_file_code,"(A,I10)") "total particles = ,", n
    write(gr_file_code,"(A,F10.2)") "reduced T = ,", temp
    write(gr_file_code,"(A,F10.5)") "reduced density = ,", density
    write(gr_file_code,"(A,I10)") "total collisions = ,", ncoll
    write(gr_file_code,"(A,F13.8)") "E_kinetic/nkT = ,", enkt1
    write(gr_file_code,"(A,F13.8)") "E_potential/nkT = ,", enkt2
    write(gr_file_code,"(A,F13.8)") "E_total/nkT = ,", enkt3
    write(gr_file_code,"(A,F13.8)") "Z = PV/NkT = ,", Z
    write(gr_file_code,"(A,F13.8)") "Z (no intrachain virial) = PV/NkT = ,", Z2
    write(gr_file_code,"(A,E15.7E2)") "Passed time = ,", time
    write(gr_file_code,"(A,I10)")"Number of ghosts,", numghost
    
    do i = 1, maxbin
        write(gr_file_code, *) i, ",", r_list(i), ",", gr(i), ","
        !bincounter = bincounter + hist(i)
    end do
    
    close(pdb_file_code)
    close(gr_file_code)
    
end program swchain
    
    
!Generate flat layer of spheres
subroutine create_coordinate()
    use global
    implicit none

    integer nc, nc2, i, j, k, m, ini_sphere 
    

    real :: chaindel, halfchaindel, nextx, nexty, nextz
    real :: xdir, ydir
    real :: new_rx, new_ry, new_rz
    real*8 :: sumx, sumy, sumz
    real :: rtemp, gauss, dummy
    real*8 :: vxsum, vysum, vzsum

	call init_random_seed()
    
    chaindel = 0.5 * (sigma1 + d2) ! half of maximum bond length 0.5sigma1*del
    halfchaindel = 0.5 * chaindel

    
    do i = 1, numchn 
        ini_sphere = 1 + chainlen * (i - 1) ! inital bead in a chain

        
        do j = 1, chainlen - 1 !do not set bond partner for last bead
            do k = 1, 2
                if ((bonded_spheres(k, ini_sphere + j - 1)) == -1) then ! (ini_sphere + j - 1) = current bad
                    bonded_spheres(k, ini_sphere + j - 1) = ini_sphere + j !set bond partner as next bead
                    EXIT
                end if
            end do
            
            do k = 1, 2
                if ((bonded_spheres(k, ini_sphere + j)) == -1) then ! (ini_sphere + j - 1) = current bead
                    bonded_spheres(k, ini_sphere + j) = ini_sphere + j - 1  !set bond partner as next bead
                    EXIT
                end if
            end do
        end do
    end do
    
    new_rx = halfchaindel - 0.5
    new_ry = 0.5 - halfchaindel
    new_rz = halfchaindel - 0.5
    
    xdir = 1.0
    ydir = -1.0
    
    rx(1) = new_rx
    ry(1) = new_ry
    rz(1) = new_rz    
    
    do i = 1, n-1
        
        new_rx = rx(i) + chaindel * xdir
        if ((new_rx * xdir) <= 0.5 - halfchaindel) then
            rx(i + 1) = new_rx
            ry(i + 1) = ry(i)
            rz(i + 1) = rz(i)
                
        else 
            new_ry = ry(i) + chaindel * ydir
            if ((new_ry * ydir) <= 0.5 - halfchaindel) then
                xdir = xdir * -1.0
                rx(i + 1) = rx(i)
                ry(i + 1) = new_ry
                rz(i + 1) = rz(i)
                    
            else
                xdir = xdir * -1.0
                ydir = ydir * -1.0                    
                nextz = (0.5 * chaindel ** 2) ** 0.5
                nextx = 0.5 * chaindel * xdir
                nexty = 0.5 * chaindel * ydir

                rx(i + 1) = rx(i) + nextx
                ry(i + 1) = ry(i) + nexty
                rz(i + 1) = rz(i) + nextz
                if (rz(i) > 0.5 - halfchaindel) then
                    write(*, *) "The system is too dense at chain",  int(i/chainlen)
                end if
            end if 
        end if

    end do
        
end subroutine create_coordinate 

subroutine create_vel(temp_temp)
    use global
    implicit none
    
    integer :: i
    real*8 :: sumx, sumy, sumz
    real :: temp_temp, rtemp, gauss, dummy
    real*8 :: vxsum, vysum, vzsum

	call init_random_seed()
    dummy = 0
    rtemp = sqrt (temp_temp)    

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
return
end subroutine create_vel

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

!Total kinetic energy
subroutine kecal(ke)
    use global
    implicit none
    integer :: i
    real :: ke
    ke = 0
    do i = 1, n
        ke = ke + vx(i)**2 + vy(i)**2 + vz(i)**2
    end do
    ke = 0.5*ke   
    return
end subroutine kecal
    
!Total potential energy
subroutine pecal (pe)
    use global
    implicit none
    
    integer :: i, j
    real :: rxi, ryi, rzi, rxij, ryij, rzij
    real :: pe, sig2sq, rijsq
    pe = 0
    sig2sq = sigma2 ** 2
    
    do i = 1, n-1
        rxi = rx(i)
        ryi = ry(i)
        rzi = rz(i)
        
        do j = i + 1, n
            rxij = rxi - rx(j)
            ryij = ryi - ry(j)
            rzij = rzi - rz(j)
            rxij = rxij - anint(rxij)
            ryij = ryij - anint(ryij)
            rzij = rzij - anint(rzij)
            rijsq = rxij**2 + ryij**2 + rzij**2
            
            if ((rijsq - sig2sq) < 0.0) then
                pe = pe - eps 
            end if
        end do
    end do
    return
end subroutine pecal
    
!check overlaps
subroutine check_position() 
    use global
    implicit none
    integer :: bonded
    real(8) :: sig1sq
    real(8) :: rxi, ryi, rzi, rxj, ryj, rzj, rxij, ryij, rzij, rijsq, rij
    
    integer :: i, j, k      ! loop values
    

    sig1sq = sigma1**2      !Square dat sigma
    overlap= .false.
    
    do i=1, n-1   
        rxi = rx(i)
        ryi = ry(i)
        rzi = rz(i)
        
        do j = i+1, n
            bonded = 0
            do k = 1, 2
                if ((bonded_spheres(k, i)) == j) then
                    bonded = 1
                    EXIT
                end if
            end do
                
            rxij = rxi - rx(j)
            ryij = ryi - ry(j)
            rzij = rzi - rz(j)
            rxij = rxij - anint(rxij)
            ryij = ryij - anint(ryij)
            rzij = rzij - anint(rzij)
            
            rijsq = rxij**2 + ryij**2 + rzij**2
            
            if (bonded == 1) then
                if ( (1-(rijsq**0.5/d1) ) > tol ) then
                    overlap = .true.
                    write (*,*) "overlap spheres in chain", rijsq**0.5, d1 
                    write (*,*) "overlap spheres in chain", i, j, coltype(i), partner(i), coltype(j), partner(j)
                else if ( ( (rijsq**0.5/d2) - 1) > tol ) then
                    overlap = .true.
                    write (*,*) "overlap spheres in chain", rijsq**0.5, d2
                    write (*,*) "chains comes apart in", i, j, coltype(i), partner(i), coltype(j), partner(j)
                end if
                    
            else 
                if (rijsq < sig1sq) then
                !write (*,*) "rij < diameter spheres", i, j
                
                    if ( (1-(rijsq**0.5/sigma1) ) > tol ) then
                        overlap = .true.
                        write (*,*) "overlap spheres", i, j
                    end if
                end if
            end if
        end do
    end do
    
    if (overlap) then
        write (*,'(A,i10,A)') "The system has overlaps at ", counter, " collision." // CHAR(10)
        stop
    !else if (overlap == .false.) then
    !    write (*, '(A,i10,A)') "The system has no overlap at ", counter, " collision."// CHAR(10)
    end if
end subroutine check_position                    
    

! calculate radial distribution total particles in every bin*********************
subroutine update_radial_bin()
    use global
    implicit none
    
    integer :: i, j, k, l, x, ini_sphere_i, ini_sphere_j
    integer :: chnid(n)
    
    real :: rxij, ryij, rzij, rijsq, rij, rij_temp
    real :: sig1sq
    !integer :: bincounter
    !bincounter = 0
    
    sig1sq = sigma1 ** 2
    
    do i = 1, n
        chnid(i) = int((i-1)/chainlen) + 1
    end do
    
    do i = 1, n-1
        do j = i+1, n
            if (chnid(i) /= chnid(j)) then
                rxij = rx(i) - rx(j)
                ryij = ry(i) - ry(j)
                rzij = rz(i) - rz(j)
                rxij = rxij - anint(rxij)
                ryij = ryij - anint(ryij)
                rzij = rzij - anint(rzij)
                rijsq = rxij**2 + ryij**2 + rzij**2
                rij = sqrt(rijsq)
                bin = int(rij/delr) + 1

                if(bin.LE.maxbin) then	!If the current bin number is less than the maximum number of bins allowed...
	                hist(bin) = hist(bin) + 2	!Update the histogram by 2 because reasons
                end if 		
            end if
        end do
    end do
    
    
    !do i = 1, maxbin
    !    write(*,*) i, hist(i)
    !    bincounter = bincounter + hist(i)
    !end do
    !
    !write(*,*) bincounter
    
end subroutine update_radial_bin
    
    
subroutine uplist(i)
    use global
    implicit none

    integer :: i, j, k   ! i 外部变量, loop values
    integer :: bonded
    real(8) :: rxi, ryi, rzi
    real(8) :: vxi, vyi, vzi
    real(8) :: rxij, ryij, rzij, rijsq, rij
    real(8) :: vxij, vyij, vzij, vijsq, vij
    real(8) :: tij, bij, bijsq, c1ij, c2ij ! bij = vector r_ij * vector v_ij, cij = rijsq - sig2sq
    real(8) :: discr1, discr2
    real(8) :: vsigma1, vsigma2
    real(8) :: sig1sq, sig2sq
    
    if (i == n) return
                                                   
    coltim(i) = timbig  ! assign maximum collision time
    tij = timbig
    coltype(i) = 5
    
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
        bijsq = bij ** 2
        
        rijsq = rxij**2 + ryij**2 + rzij**2
        vijsq = vxij**2 + vyij**2 + vzij**2
        
        bonded = 0
        do k = 1, 2
            if (bonded_spheres(k, i) == j) then
                bonded = 1
                EXIT
            end if
        end do        
        

        if (bonded == 1) then
            vsigma1 = d1
            vsigma2 = d2
            sig1sq = vsigma1 ** 2
            sig2sq = vsigma2 ** 2        
        
            c1ij = rijsq - sig1sq
            c2ij = rijsq - sig2sq
            
            if (bij < 0.0) then ! center approaching
                discr1 = bijsq - vijsq * c1ij
                  ! j inside the sigma 2
                if (discr1 > 0) then     
                    tij = (-bij - sqrt(discr1))/vijsq
                    
                    if (tij < coltim(i)) then
                        coltim(i) = tij
                        partner(i) = j
                        coltype(i) = 6                  ! bonded core collision
                    end if
                
                else                                    ! discr1 <= 0 core does not collide, and (c2ij < 0) (discr2 must <0) 
                    discr2 = bijsq - vijsq * c2ij 
                    tij = (-bij + sqrt(discr2))/vijsq
                    
                    if (tij < coltim(i)) then           
                        coltim(i) = tij
                        partner(i) = j
                        coltype(i) = 7                  ! center approaching bonded bounce
                    end if
                    
                end if
                
            end if
            
            if (bij >= 0.0) then ! center recedes AND j inside the d2n 
                discr2 = bijsq - vijsq * c2ij 
                tij = (-bij + sqrt(discr2))/vijsq
                    
                if (tij < coltim(i)) then           
                    coltim(i) = tij
                    partner(i) = j
                    coltype(i) = 7                  ! center receding bonded bounce
                end if
            end if
        
        else 
            vsigma1 = sigma1
            vsigma2 = sigma2
            sig1sq = vsigma1 ** 2
            sig2sq = vsigma2 ** 2        
        
            c1ij = rijsq - sig1sq
            c2ij = rijsq - sig2sq
        
            if (bij < 0.0) then ! center approaching
                discr1 = bijsq - vijsq * c1ij
                if (c2ij <= 0) then  ! j inside the sigma 2
                    if (discr1 > 0) then     
                        tij = (-bij - sqrt(discr1))/vijsq
                    
                        if (tij < coltim(i)) then
                            coltim(i) = tij
                            partner(i) = j
                            coltype(i) = 1                  ! core collision
                        end if
                
                    else                                    ! discr1 <= 0 core does not collide, and (c2ij < 0) (discr2 must <0) 
                        discr2 = bijsq - vijsq * c2ij 
                        tij = (-bij + sqrt(discr2))/vijsq
                    
                        if (tij < coltim(i)) then           
                            coltim(i) = tij
                            partner(i) = j
                            coltype(i) = 2                  ! bounce or dissociate
                        end if
                    
                    end if
            
                else                                        ! (c2ij > 0) outside the sigma2
                    discr2 = bijsq - vijsq * c2ij
                    if (discr2 >= 0) then
                        tij = (-bij - sqrt(discr2))/vijsq
                        if (tij < coltim(i)) then           
                            coltim(i) = tij
                            partner(i) = j
                            coltype(i) = 4                  ! capture
                        end if
                    
                    end if                                  ! (discr2 <= 0) no attractive collision
                end if    
            end if
            
            if ((bij > 0.0) .AND. (c2ij < 0)) then ! center recedes
                if (c2ij < 0) then ! j inside the sigma 2
                    discr2 = bijsq - vijsq * c2ij 
                    tij = (-bij + sqrt(discr2))/vijsq
                    
                    if (tij < coltim(i)) then           
                        coltim(i) = tij
                        partner(i) = j
                        coltype(i) = 2                  ! bounce or dissociate
                    end if
                end if
            end if
        
        end if
        
    end do
    
    return
end subroutine uplist

    
subroutine dnlist(j)
    use global
    implicit none

    integer :: i, j, k   ! i 外部变量, loop values
    integer :: bonded
    real(8) :: rxj, ryj, rzj
    real(8) :: vxj, vyj, vzj
    real(8) :: rxij, ryij, rzij, rijsq, rij
    real(8) :: vxij, vyij, vzij, vijsq, vij
    real(8) :: tij, bij, bijsq, c1ij, c2ij ! bij = vector r_ij * vector v_ij, cij = rijsq - sig2sq
    real(8) :: discr1, discr2
    real(8) :: vsigma1, vsigma2
    real(8) :: sig1sq, sig2sq
    
    if (j == 1) return
    
    rxj = rx(j)
    ryj = ry(j)
    rzj = rz(j)
    vxj = vx(j)
    vyj = vy(j)
    vzj = vz(j)
    
    do i = 1, j - 1
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
        bijsq = bij ** 2
        
        rijsq = rxij**2 + ryij**2 + rzij**2
        vijsq = vxij**2 + vyij**2 + vzij**2
        
        bonded = 0
        do k = 1, 2
            if (bonded_spheres(k, j) == i) then
                bonded = 1
                EXIT
            end if
        end do        
        
        if (bonded == 1) then
            vsigma1 = d1
            vsigma2 = d2
            sig1sq = vsigma1 ** 2
            sig2sq = vsigma2 ** 2        
        
            c1ij = rijsq - sig1sq
            c2ij = rijsq - sig2sq
            
            if (bij < 0.0) then ! center approaching
                discr1 = bijsq - vijsq * c1ij 
 
                if (discr1 > 0) then     
                    tij = (-bij - sqrt(discr1))/vijsq
                    
                    if (tij < coltim(i)) then
                        coltim(i) = tij
                        partner(i) = j
                        coltype(i) = 6                  ! bonded core collision
                    end if
                
                else                                ! discr1 <= 0 core does not collide, and (c2ij < 0) (discr2 must <0) 
                    discr2 = bijsq - vijsq * c2ij 
                    tij = (-bij + sqrt(discr2))/vijsq
                    
                    if (tij < coltim(i)) then           
                        coltim(i) = tij
                        partner(i) = j
                        coltype(i) = 7                  ! center approaching bonded bounce
                    end if
                    
                end if
                
            end if
            
            if (bij >= 0.0) then ! center recedes AND j inside the d2n 
                discr2 = bijsq - vijsq * c2ij 
                tij = (-bij + sqrt(discr2))/vijsq
                    
                if (tij < coltim(i)) then           
                    coltim(i) = tij
                    partner(i) = j
                    coltype(i) = 7                  ! center receding bonded bounce
                end if
            end if
            

        else 
            vsigma1 = sigma1
            vsigma2 = sigma2
            sig1sq = vsigma1 ** 2
            sig2sq = vsigma2 ** 2        
        
            c1ij = rijsq - sig1sq
            c2ij = rijsq - sig2sq
        
            if (bij < 0.0) then ! center approaching
                discr1 = bijsq - vijsq * c1ij 
                if (c2ij <= 0) then  ! j inside the sigma 2
                    if (discr1 > 0) then     
                        tij = (-bij - sqrt(discr1))/vijsq
                    
                        if (tij < coltim(i)) then
                            coltim(i) = tij
                            partner(i) = j
                            coltype(i) = 1                  ! core collision
                        end if
                
                    else                                    ! discr1 <= 0 core does not collide, and (c2ij < 0) (discr2 must <0) 
                        discr2 = bijsq - vijsq * c2ij 
                        tij = (-bij + sqrt(discr2))/vijsq
                    
                        if (tij < coltim(i)) then           
                            coltim(i) = tij
                            partner(i) = j
                            coltype(i) = 2                  ! bounce or dissociate
                        end if
                    
                    end if
            
                else                                        ! (c2ij > 0) outside the sigma2
                    discr2 = bijsq - vijsq * c2ij
                    if (discr2 >= 0) then
                        tij = (-bij - sqrt(discr2))/vijsq
                        if (tij < coltim(i)) then           
                            coltim(i) = tij
                            partner(i) = j
                            coltype(i) = 4                  ! capture
                        end if
                    
                    end if                                  ! (discr2 <= 0) no attractive collision
                end if    
            end if
            
            if ((bij > 0.0) .AND. (c2ij < 0)) then ! center recedes
                if (c2ij < 0) then ! j inside the sigma 2
                    discr2 = bijsq - vijsq * c2ij 
                    tij = (-bij + sqrt(discr2))/vijsq
                    
                    if (tij < coltim(i)) then           
                        coltim(i) = tij
                        partner(i) = j
                        coltype(i) = 2                  ! bounce or dissociate
                    end if
                end if
            end if
        
        end if
    
    end do
    
    return
    end subroutine dnlist

subroutine bump()
    use global 
    implicit none
    real(8), parameter :: smdist = 5e-12

    integer :: i, j, k, bonded
    real(8) :: rxij, ryij, rzij, vxij, vyij, vzij
    real(8) :: delvx, delvy, delvz
    real(8) :: vsigma1, vsigma2
    real(8) :: bij, bijsq, rijsq, sig1sq, sig2sq, discr_judge
    real(8) :: factor, smbump
    
    bonded = 0
    i = current_i
    j = current_j
    
    do k = 1, 2
        if (bonded_spheres(k, i) == j) then
            bonded = 1
            EXIT
        end if
    end do
    
    if (bonded == 1) then           !bond on the bonding info, set d1 d2 to sigma1 sigma2
        vsigma1 = d1
        vsigma2 = d2
    else
        vsigma1 = sigma1
        vsigma2 = sigma2
    end if
    
    sig1sq = vsigma1 ** 2
    sig2sq = vsigma2 ** 2
    
    rxij = rx(i) - rx(j)
    ryij = ry(i) - ry(j)
    rzij = rz(i) - rz(j)
    rxij = rxij - anint(rxij)
    ryij = ryij - anint(ryij)
    rzij = rzij - anint(rzij)
    rijsq = rxij**2 + ryij**2 + rzij**2
    
    vxij = vx(i) - vx(j)
    vyij = vy(i) - vy(j)
    vzij = vz(i) - vz(j)
    
    bij = rxij * vxij + ryij * vyij + rzij * vzij
    bijsq = bij ** 2
    discr_judge = bijsq - (4 * sig2sq * eps)
    
    smbump = 0.0
    

        select case (coltype(i))
            case (6) ! core collision
                factor = bij/sig1sq
                smbump = 0.0
            case (7) 
                factor = bij/sig2sq                                             ! bonded spheres bounce
                smbump = -1.0 
            case (1) ! core collision
                factor = bij/sig1sq
                smbump = 0.0
            case (2) ! bounce/dissociate
                if (discr_judge <= 0) then 
                    factor = bij/sig2sq                                          ! non-bonded spheres bounce
                    smbump = -1.0
                else ! (discr_judge >= 0)
                    factor = (-sqrt(- 4 * sig2sq * eps + bijsq) + bij)/(2 * sig2sq) ! dissociation
                    smbump = 1.0
                end if
            
            case (4) ! capture
                factor = (sqrt(4 * sig2sq * eps + bijsq) + bij)/(2 * sig2sq)
                smbump = -1.0
        end select      

    
    delvx = -factor * rxij
    delvy = -factor * ryij
    delvz = -factor * rzij
    
    vx(i) = vx(i) + delvx
    vy(i) = vy(i) + delvy
    vz(i) = vz(i) + delvz
    
    vx(j) = vx(j) - delvx
    vy(j) = vy(j) - delvy
    vz(j) = vz(j) - delvz
    
    w = delvx * rxij + delvy * ryij + delvz * rzij
    acw = acw + w
    
    if (bonded /= 1) then
        acw2 = acw2 + w
    end if

        
    if (smbump == -1.0) then
        do while (rijsq >= sig2sq) ! for bounce and capture, |rij| should < sigma2

            rx(i) = rx(i) + smbump * smdist * rxij * vsigma2
            ry(i) = ry(i) + smbump * smdist * ryij * vsigma2
            rz(i) = rz(i) + smbump * smdist * rzij * vsigma2
            rx(j) = rx(j) - smbump * smdist * rxij * vsigma2
            ry(j) = ry(j) - smbump * smdist * ryij * vsigma2
            rz(j) = rz(j) - smbump * smdist * rzij * vsigma2
            rxij = rx(i) - rx(j)
            ryij = ry(i) - ry(j)
            rzij = rz(i) - rz(j)   
            rxij = rxij - anint(rxij)
            ryij = ryij - anint(ryij)
            rzij = rzij - anint(rzij)
            rijsq = rxij**2 + ryij**2 + rzij**2
            !write(*,*) smbump, rijsq-sig2sq
            !write(*,"(7F15.10)") rx(i), ry(i), rz(i), rx(j), ry(j), rz(j), smbump * smdist * rxij * sigma2
        end do
    else if (smbump == 1.0) then 
        do while (rijsq <= sig2sq) ! for dissociation, |rij| should > sigma2

            rx(i) = rx(i) + smbump * smdist * rxij * vsigma2
            ry(i) = ry(i) + smbump * smdist * ryij * vsigma2
            rz(i) = rz(i) + smbump * smdist * rzij * vsigma2
            rx(j) = rx(j) - smbump * smdist * rxij * vsigma2
            ry(j) = ry(j) - smbump * smdist * ryij * vsigma2
            rz(j) = rz(j) - smbump * smdist * rzij * vsigma2
            rxij = rx(i) - rx(j)
            ryij = ry(i) - ry(j)
            rzij = rz(i) - rz(j)   
            rxij = rxij - anint(rxij)
            ryij = ryij - anint(ryij)
            rzij = rzij - anint(rzij)
            rijsq = rxij**2 + ryij**2 + rzij**2
            !write(*,"(7F15.10)") rx(i), ry(i), rz(i), rx(j), ry(j), rz(j), smbump * smdist * rxij * sigma2
            !write(*,*) smbump, rijsq-sig2sq
            
        end do
    end if
    
end subroutine bump
    
subroutine write_pdb_file(file_code, pdbfile)
    use global
    implicit none
    integer :: file_code, i
    character :: pdbfile*(*)
    

    write(file_code,'(A, I5)') "MODEL", n
        do i = 1, n
           write(file_code,'(A6,I5,A3,16X,3F8.3)') 'ATOM  ',I, 'N', rx(i), ry(i), rz(i)
        end do
        
    write(file_code,'(A)') "ENDMDL"   
        
    return
end subroutine write_pdb_file       

subroutine ghostcoll(ghost)
    use global
    implicit none
    integer :: ghost
    real :: rtemp, gauss, dummy

	call init_random_seed()
    dummy = 0
    rtemp = sqrt (destemp)    
  
    vx(ghost) = rtemp * gauss (dummy)
    vy(ghost) = rtemp * gauss (dummy)
    vz(ghost) = rtemp * gauss (dummy)
    
    return
end subroutine ghostcoll
