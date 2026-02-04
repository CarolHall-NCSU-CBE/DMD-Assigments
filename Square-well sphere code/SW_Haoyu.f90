module global
        implicit none
        save
        integer, parameter  :: n=500
        
        real(8), parameter  :: pi = 3.14159265359
        real, parameter :: lambda = 1.5     ! 引力半径与真实半径倍数
        real, parameter :: eps = 1.0        ! square-well potential parameter
        !real, parameter  :: box = 1.0, rlower = 0.0, rupper = box/2.0 
        real, parameter  :: timbig = 1.0e10                         !Declare and intiialize max time for minimum collision time
        real, parameter :: tol = 1.0E-4					        	!应该是说明 重叠的区域不超过 0.0001
        real, parameter :: maxbin = 500     !radial distribution 图表的总bin数
        real, parameter :: delr = 0.001     !radial distribution 的每个bin的宽度，即step, the delta r (size of the RDF shell)
        real, parameter :: ghostcoeff = 0.01
                
        real(8) :: series(n), rx(n), ry(n), rz(n), vx(n), vy(n), vz(n)
        real :: coltim(n+2)                                         ! Declare collision time array
        
        real(8) :: sigma1, sigma2                                    !*** sigma1 = hardsphere diameter, sigma2 = square well diameter
        
        real(8) :: time, w, Z, density, destemp ! total time, virial per collision, compressibility factor,reduced density, reduced temperature
        real :: acw = 0                                             ! cumulative virial used for calculation Z=RV/NkT

        integer :: partner(n), coltype(n)                           ! partner array, collision type array
        integer :: ncoll, current_i, current_j, counter             ! 碰撞总数, 当前碰撞粒子 i j, 碰撞计数器
        integer :: ngcoll_counter                                   ! tij 为负数counter
        integer :: bin, grcounter                                   ! 柱状图柱子序号, gr 计数器
        integer :: hist(maxbin)                                     ! histogram for radial distribution

        real :: delv                                          ! radial distribution differential volume
        real :: gr_density                                  ! density at r
        real :: gr(maxbin), r_list(maxbin)      ! *gr function*, *gr axis* 
            
        logical :: overlap
        
end module global
    
program squarewell
    use global
    implicit none
    integer, parameter :: pdb_file_code = 30
    integer, parameter :: xyz_file_code = 11
    integer, parameter :: rv_file_code = 12
    integer, parameter :: gr_file_code = 13
    integer, parameter :: temp_file_code = 15
    real, parameter  :: inittime = 0.001
    
    integer :: i, j, k, status      ! loop value
    integer :: numghost, ghost      ! number of ghost picked, random pick a sphere as ghost
    
    real(8) :: tij, rij2, sig2sq, bij
    real(8) :: rxij, ryij, rzij, vxij, vyij, vzij
    real :: ran1, ran2, ran3        ! random number
    real :: ke, pe                  ! kinetic energy, potential energy
    real :: ghosttime               ! time to perform ghost collision
    real :: r1, r2                  ! small and big radius for radial distribution differential volume
    real :: e, en, enkt, pvnkt1, temp, tbc	 ! total energy, single particle energy
    character(10) :: temp_s, n_s, density_s  ! open coordinate file 12.23
    character(40) :: rv_filename, pdb_filename, temp_filename, gr_filename ! open coordinate file 12.23
    character(len=8) :: date1, date2 
    character(len=20) :: world_time1, world_time2
    character(len=6) :: world_time1_trim, world_time2_trim
    !REAL :: density, dij, tij, t, rate					!Declare a bunch of real/float type variables
    
    call init_random_seed()
    
    !********enter parameters **********************************************************
    write (*,*) "Enter reduced temperature kT/epsilon (suggested 1.5)"
    read (*,*) destemp                                      !Input temperature
    write (*,*) "Enter reduced density (n/v)*sigma**3"	
    read (*,*) density										!Input density
    write (*,*) "Enter number of collisions required"		
    read (*,*) ncoll										!Input total collision number
    
    !********file name strings **********************************************************
    call date_and_time(date1, world_time1)
    world_time1_trim = world_time1(1:6)
    write(temp_s, '(F5.2)') destemp
    write(n_s, '(I4)') n
    write(density_s, '(F5.2)') density
    
    !*********CREATE coordinates and initial velocities **********************************************************
    call coordinate_and_vel()
    
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

    !********initialize parameters **********************************************************
    sigma1 = (density/REAL(n))**(1.0/3.0)                    
    sigma2 = lambda * sigma1
    numghost = 0
    time = 0.0
    counter = 0 
    hist = 0                                        ! initialize histogram for countering particles in certain bin (distance)
    grcounter = 0
    ngcoll_counter = 0
    call random_number(ran1)
    call random_number(ran2)
    ghosttime = ghostcoeff * abs(sqrt(-2 * alog(ran1)) * cos(2 * pi * ran2))
    coltim (n+1) = ghosttime
    coltim (n+2) = inittime
    
    !********Open output files **********************************************************
    pdb_filename = 'SW_pdb_p'//adjustl(trim(density_s))//'_T'//adjustl(trim(temp_s))//'_n'//adjustl(trim(n_s))//'.pdb'
    temp_filename = 'SW_temp_p'//trim(adjustl(trim(density_s)))//'_T'//trim(adjustl(trim(temp_s)))//'_n'//trim(adjustl(trim(n_s)))//'.csv'
    !open(unit=pdb_file_code, file=pdb_filename, status='replace', action='write')
    open(unit=temp_file_code, file=temp_filename, status='replace', action='write')
    
    !*******recalculate temperature **********************************************************
    call kecal(ke)              !Calculate the kinetic energy
    call pecal(pe)              !Calculate the potential energy
    !en = e/REAL(n)
    temp = 2.0/3.0 * ke/real(n)	!Calculate the temperature

    write (*, "(A,F10.5)") "Calculated Reduced Temperature:", temp

    call check_position()     ! check overlaps

    !***************do initial collision **********************************************************

    do i = 1, n                          
        coltim(i) = timbig              ! assign maximum collision time
        coltype(i) = 5
        partner(i) = n
    end do
    do i = 1, n - 1
        call uplist(i)					! initial collision
    end do
    counter = counter + 1
    ! * write in pdb file *
    
    !call write_pdb_file(pdb_file_code, pdb_filename)
    
    !***************************************
    !**************Main Loop****************
    !***************************************
    do while (counter < ncoll)                              ! main loop
        !**************Update full list every 100 time ****************
        if (mod(counter, 100) == 0) then
            do k = 1, n
                call uplist(k)
            end do
        end if
    
        !!***************** Write Temperature vs Coll ********************
        if (mod(counter, 20) == 0) then
           call kecal (ke)
           temp = 2.0 * ke / 3.0 / real (n)
           write(temp_file_code, "(I10,A,F15.10,A,F15.8)") counter, ",", time, ",", temp
        endif
    
        !**************Find t_min and negative collision times****************
        tij = timbig           ! set tij to a large value
        do k = 1, n + 1        ! find i j with smallest tij
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
        time = time + tij   ! 总时间
        coltim(n+1) = coltim(n+1) - tij      !ghost particle reduce by tij，update positions for all spheres

        !**************Move particles by t_min****************
        do k = 1, n                         
            coltim(k) = coltim(k) - tij     ! reduce collision by tij，update positions for all spheres
            rx(k) = rx(k) + vx(k) * tij
            ry(k) = ry(k) + vy(k) * tij
            rz(k) = rz(k) + vz(k) * tij
            rx(k) = rx(k) - anint(rx(k))    
            ry(k) = ry(k) - anint(ry(k))
            rz(k) = rz(k) - anint(rz(k))
        end do
        
        if (mod(counter, 1000) == 0) then
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
            ghost = int(ran3 * real(n)) + 1 ! ghost particle is picked randomly from all spheres
            call ghostcoll(ghost)
            ghosttime = ghostcoeff * abs(sqrt(-2 * alog(ran1)) * cos(2 * pi * ran2))
            coltim (n+1) = ghosttime    ! reset ghost collision time
            
            do k = 1, n
                 if ((k == ghost) .or. (partner(k) == ghost)) then
                    call uplist(k)
                 endif
            end do
            call dnlist (ghost)
            numghost = numghost + 1
        end if
        
        if (mod(counter, 20) == 0 ) then    !record radial distribution every 20 collisions
            call update_radial_bin()
            grcounter = grcounter + 1
        end if
        
        !call write_xyz_file(xyz_file_code)
        !call write_pdb_file(pdb_file_code, pdb_filename)
        counter = counter + 1
    end do
    !****************************************
    !***********end collisions **************
    !****************************************
    
    call check_position()
    
    Z = 1.0 + (acw / (real(n) * 3.0 * time * temp))                ! 计算压缩因子
    
    call kecal(ke)
    enkt = ke/(REAL(n) * temp)                            !Calculate E/kT

    
    write(*,"(A,F13.8)") "E/nkT = ", enkt
    write(*,"(A,F13.8)") "Z = PV/NkT = ", Z
    write(*,"(A,E15.7E2)") "Passed time =", time
    write(*,"(A,I5)") "Number of ghost =", numghost
    
    !calculation radial distribution **********************************************************
    do i = 1, maxbin
        r1 = (i - 1) * delr 
        r2 = i * delr
        
        delv = 4.0 * pi / 3.0 * (r2**3 - r1**3)        
        gr_density = real(hist(i)) / (real(grcounter) * real(n) * delv)           ! the g(r) is defined as (density at r)/(uniform density)
        gr(i) = gr_density / (density/sigma1**3)                 ! reduce density/sigma**3 == real density = n/V = n here(V= 1^3) 
        r_list(i) = r1 + 0.5 * delr
        
    end do
    
    !export radial distribution data **********************************************************
    gr_filename = 'SW_gr_p'//trim(adjustl(trim(density_s)))//'_T'//trim(adjustl(trim(temp_s)))//'_n'//trim(adjustl(trim(n_s)))//'.csv'
    
    
    open(unit=gr_file_code, file=gr_filename, status='replace', action='write')   ! 导出radial distribution data
    call date_and_time(date2, world_time2)    !这里是 提取run之前的时间
    world_time2_trim = world_time2(1:6)
    
    write(gr_file_code,"(5A)") "Before simulation,", "Date,", date1, ", Time,", world_time1_trim
    write(gr_file_code,"(5A)") "After simulation,", "Date,", date2, ", Time,", world_time2_trim
    write(gr_file_code,"(A,I10)") "total particles = ,", n
    write(gr_file_code,"(A,F10.2)") "reduced T = ,", temp
    write(gr_file_code,"(A,F10.5)") "reduced density = ,", density
    write(gr_file_code,"(A,I10)") "total collisions = ,", ncoll
    write(gr_file_code,"(A,F13.8)") "E/nkT = ,", enkt
    write(gr_file_code,"(A,F13.8)") "Z = PV/NkT = ,", Z
    write(gr_file_code,"(A,E15.7E2)") "Passed time = ,", time
    
    do i = 1, maxbin
        write(gr_file_code, *) i, ",", r_list(i), ",", gr(i), ","
        !bincounter = bincounter + hist(i)
    end do
    
    close(pdb_file_code)
    close(gr_file_code)
    
end program squarewell
    
    
!用于生成ffc晶格
subroutine coordinate_and_vel()
    use global
    implicit none

    integer nc, nc2, i, j, k, m      !number of cell on one side
    real cell_len, cell_len2    
    real*8 :: sumx, sumy, sumz
    real :: rtemp, gauss, dummy
    real*8 :: vxsum, vysum, vzsum
	
	call init_random_seed() !随机数生成初始化 only need call once 
    
! 用于产生 4*(3^3) 一个晶格包含4个，3x3的ffc晶格 ***********************************
    nc = (n/4)**(1.0/3.0) 
    nc2 = nc**2
    
	cell_len = 1.0/real(nc) !晶格长度
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
! 用于产生初始速度***********************************
    dummy = 0 !傀儡变量，用于高斯分布函数
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
        
end subroutine coordinate_and_vel   
    
! 高斯分布函数，输入随机数，产生 平均数为0，偏差为1的随机分布值***********************************    
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
    
! 用于初始化随机数seed，使其完全不同*********************************** 
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

!Calculate kinetic energy**************
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
    
!Calculate potential energy**************
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
    
 !用于检测初始位置************************************************************
subroutine check_position() 
    use global
    implicit none
    
    real(8) :: sig1sq
    real(8) :: rxi, ryi, rzi, rxj, ryj, rzj, rxij, ryij, rzij, rijsq, rij
    
    integer :: i, j                                                    ! loop values
    
    sig1sq = sigma1**2												!Square dat sigma
    overlap= .false.
    
    do i=1, n-1   
        do j = i+1, n

            rxij = rx(i) - rx(j)
            ryij = ry(i) - ry(j)
            rzij = rz(i) - rz(j)
            rxij = rxij - anint(rxij)
            ryij = ryij - anint(ryij)
            rzij = rzij - anint(rzij)
            
            rijsq = rxij**2 + ryij**2 + rzij**2
            
            if (rijsq < sig1sq) then !先判断两个粒子之间的距离是不是比粒子直径小
                !write (*,*) "rij < diameter spheres", i, j
                
                if ( (1-(rijsq**0.5/sigma1) ) > tol ) then
                    overlap = .true.
                    write (*,*) "overlap spheres", i, j
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
    
    real :: rxij, ryij, rzij, rijsq, rij
    real :: sig1sq
    !integer :: bincounter
    !bincounter = 0
    
    sig1sq = sigma1 ** 2
    
    do i = 1, n-1
        do j = i+1, n
            rxij = rx(i) - rx(j)
            ryij = ry(i) - ry(j)
            rzij = rz(i) - rz(j)
            rxij = rxij - anint(rxij)
            ryij = ryij - anint(ryij)
            rzij = rzij - anint(rzij)
    
            rijsq = rxij**2 + ryij**2 + rzij**2         ! calculate distance between all pairs
            rij = sqrt(rijsq)
            bin = int(rij/delr) + 1                     ! determine the pair belongs to which bin

            if(bin.LE.maxbin) then                          
	            hist(bin) = hist(bin) + 2               !Update the histogram by 2, updates n times for one spheres
            end if                                      !later need to divide by "n"

        end do
    end do
end subroutine update_radial_bin
    
    
subroutine uplist(i)
    use global
    implicit none

    integer :: i, j   ! i external variable, loop values
    real(8) :: rxi, ryi, rzi
    real(8) :: vxi, vyi, vzi
    real(8) :: rxij, ryij, rzij, rijsq, rij
    real(8) :: vxij, vyij, vzij, vijsq, vij
    real(8) :: tij, bij, bijsq, c1ij, c2ij ! bij = vector r_ij * vector v_ij, cij = rijsq - sig2sq
    real(8) :: discr1, discr2
    real(8) :: sig1sq, sig2sq
    
    if (i == n) return
    sig1sq = sigma1 ** 2
    sig2sq = sigma2 ** 2
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
        
        c1ij = rijsq - sig1sq
        c2ij = rijsq - sig2sq
        
        !write(*,"(I5,I5, A, F15.8,A, F15.8 )") i, j,"c2ij", c2ij,"bij", bij
        
        
        if (bij < 0.0) then ! center approaching
            discr1 = bijsq - vijsq * c1ij ! quadratic determinant 1
            if (c2ij <= 0) then  ! j inside the sigma 2
                if (discr1 > 0) then     
                    tij = (-bij - sqrt(discr1))/vijsq
                    
                    if (tij < coltim(i)) then           ! find smallest collision time
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
            
        if (bij > 0.0) then ! center recedes
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
    
    end do
    
    return
end subroutine uplist

    
subroutine dnlist(j)
    use global
    implicit none

    integer :: i, j   ! i 外部变量, loop values
    real(8) :: rxj, ryj, rzj
    real(8) :: vxj, vyj, vzj
    real(8) :: rxij, ryij, rzij, rijsq, rij
    real(8) :: vxij, vyij, vzij, vijsq, vij
    real(8) :: tij, bij, bijsq, c1ij, c2ij ! bij = vector r_ij * vector v_ij, cij = rijsq - sig2sq
    real(8) :: discr1, discr2
    real(8) :: sig1sq, sig2sq
    
    if (j == 1) return
    sig1sq = sigma1 ** 2
    sig2sq = sigma2 ** 2


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
        
        c1ij = rijsq - sig1sq
        c2ij = rijsq - sig2sq
        
        if (bij < 0.0) then ! center approaching
            discr1 = bijsq - vijsq * c1ij
            if (c2ij <= 0) then  ! j inside the sigma 2
                if (discr1 > 0) then     
                    tij = (-bij - sqrt(discr1))/vijsq
                    
                    if (tij < coltim(i)) then           ! find smallest collision time
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
            
        if (bij > 0.0) then ! center recedes
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
    
    end do
    
    return
    end subroutine dnlist

subroutine bump()
    use global 
    implicit none
    real(8), parameter :: smdist = 5e-8

    integer :: i, j, k
    real(8) :: rxij, ryij, rzij, vxij, vyij, vzij
    real(8) :: delvx, delvy, delvz
    real(8) :: bij, bijsq, rijsq, sig1sq, sig2sq, discr_judge
    real(8) :: factor, smbump
    
    i = current_i
    j = current_j
    sig1sq = sigma1 ** 2
    sig2sq = sigma2 ** 2
    
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
        case (1) ! core collision
            factor = bij/sig1sq
            smbump = 0.0
        case (2:3) ! bounce/dissociate
            if (discr_judge <= 0) then 
                factor = bij/sig2sq                                             ! bounce
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
    
	! small bumps to move spheres apart or into the square-well
    if (smbump == -1.0) then
        do while (rijsq >= sig2sq) ! for bounce and capture, |rij| should < sigma2
        !do k = 1, 100
            rx(i) = rx(i) + smbump * smdist * rxij * sigma2
            ry(i) = ry(i) + smbump * smdist * ryij * sigma2
            rz(i) = rz(i) + smbump * smdist * rzij * sigma2
            rx(j) = rx(j) - smbump * smdist * rxij * sigma2
            ry(j) = ry(j) - smbump * smdist * ryij * sigma2
            rz(j) = rz(j) - smbump * smdist * rzij * sigma2
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
        !do k = 1, 100
            rx(i) = rx(i) + smbump * smdist * rxij * sigma2
            ry(i) = ry(i) + smbump * smdist * ryij * sigma2
            rz(i) = rz(i) + smbump * smdist * rzij * sigma2
            rx(j) = rx(j) - smbump * smdist * rxij * sigma2
            ry(j) = ry(j) - smbump * smdist * ryij * sigma2
            rz(j) = rz(j) - smbump * smdist * rzij * sigma2
            rxij = rx(i) - rx(j)
            ryij = ry(i) - ry(j)
            rzij = rz(i) - rz(j)   
            rxij = rxij - anint(rxij)
            ryij = ryij - anint(ryij)
            rzij = rzij - anint(rzij)
            rijsq = rxij**2 + ryij**2 + rzij**2

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