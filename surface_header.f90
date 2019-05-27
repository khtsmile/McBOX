module surface_header

    use constants!, only : INFINITY, TOOLONG, EPSILON
    use omp_lib
    
    implicit none 
    
    
        
    type Surface
        character(20) :: surf_id            !> surface id
        integer       :: surf_type            
        integer       :: bc                 !> 0 : normal
                                            !> 1 : vacuum
                                            !> 2 : reflective
        
        real(8), dimension(5) :: parmtrs    !> p   : 
                                            !> sqr : 
                                            !> cyl : 
                                            !> sph : 

        integer, allocatable :: &
             neighbor_pos(:), &               ! List of cells on positive side
             neighbor_neg(:)                  ! List of cells on negative side
            
    contains 
        !procedure :: find_idx => find_surf_idx
    
    end type
    type(Surface), allocatable, target :: Surfaces(:), Surfaces_temp(:)
    type(surface) :: surf_obj
    

    
    
    contains
    
    subroutine read_surf (surfobj, line)
        type(Surface) :: surfobj
        character(*) :: line
        character(30):: temp, option, surf_type
        integer :: i, j, k
    
        j = 0; k = 0 
        do while (j.lt.len(line)) 
            j = j+1; k = k+1 
            if (line(j:j).ne.' ') then 
                call process_line (line, j, temp ) 
                if (k.eq.1) read(temp, *) option
                if (k.eq.2) read(temp, *) surfobj%surf_id
                if (k.eq.3) read(temp, *) surf_type
                if (k.eq.4) read(temp, *) surfobj%parmtrs(1)
                if (k.eq.5) read(temp, *) surfobj%parmtrs(2)
                if (k.eq.6) read(temp, *) surfobj%parmtrs(3)
                if (k.eq.7) read(temp, *) surfobj%parmtrs(4)
                if (k.eq.8) read(temp, *) surfobj%parmtrs(5)
            endif 
        enddo
        surfobj%surf_type = surf_type_converter(surf_type)
        surfobj%bc          = 0 ! by default 
    end subroutine 
    
    
    subroutine read_bc (surflist, line)
        type(Surface) :: surflist(:)
        character(*) :: line
        character(30):: temp, option, surf_id
        integer :: i, j, k, idx
    
        j = 0; k = 0 
        do while (j.lt.len(line)) 
            j = j+1; k = k+1 
            if (line(j:j).ne.' ') then 
                call process_line (line, j, temp ) 
                
                if (k.eq.1) read(temp, *) option
                if (k.eq.2) then 
                    read(temp, *) surf_id
                    idx = find_surf_idx(surflist, surf_id)
                endif
                if (k.eq.3) read(temp, *) surflist(idx)%bc  
            endif 
        enddo
        
    end subroutine 
    
    
    
    function surf_type_converter(surf_type) result(output) 
        character(*), intent(in) :: surf_type
        integer :: output
        output = 0 
        
        if (surf_type.eq.'px')   output = 1 
        if (surf_type.eq.'py')   output = 2
        if (surf_type.eq.'pz')   output = 3 
        if (surf_type.eq.'sqcx') output = 4
        if (surf_type.eq.'sqcy') output = 5
        if (surf_type.eq.'sqcz') output = 6
        if (surf_type.eq.'cylx') output = 7
        if (surf_type.eq.'cyly') output = 8
        if (surf_type.eq.'cylz') output = 9
        if (surf_type.eq.'sph')  output = 10
        
        if (output == 0) then 
            print *, '******************************************'
            print *, 'ERROR: WRONG SURFACE TYPE - ', surf_type 
            print *, '******************************************'
            stop  
        endif 
        
        
    end function 
    
    function find_surf_idx(this, surf_id) result (idx)
        type(surface) :: this(:)
        integer ::  i, idx
        character(*):: surf_id
        
        do i = 1, size(this) 
            if (this(i)%surf_id == surf_id) then 
                idx = i 
                goto 1
            endif
        enddo 
        print *, "no such surface id : ", surf_id 
        stop 
1        continue
    end function 


    function surf_neg_or_pos(this, xyz) result(neg) 
        type(surface) :: this
        real(8) :: xyz(3), xyz_tr(3), val
        logical :: neg
        integer :: type
        
        neg = .false.
        
        type = this % surf_type
        select case (type)
        case(1) !> px
            if (xyz(1) < this%parmtrs(1)) neg = .true.
        case(2) !> py
            if (xyz(2) < this%parmtrs(1)) neg = .true.
        case(3) !> pz
            if (xyz(3) < this%parmtrs(1)) neg = .true.
            
        case(6) !> sqcz
            xyz_tr(1) = xyz(1) - this%parmtrs(1)
            xyz_tr(2) = xyz(2) - this%parmtrs(2)
            xyz_tr(3) = xyz(3)
            if ((abs(xyz_tr(1))) < this%parmtrs(3)&
                .and.abs(xyz_tr(2)) < this%parmtrs(3))  &
                neg = .true.
            
        case(9) !> cylz
            xyz_tr(1) = xyz(1) - this%parmtrs(1)
            xyz_tr(2) = xyz(2) - this%parmtrs(2)
            xyz_tr(3) = xyz(3)
            
            val = sqrt((xyz_tr(1))**2 + (xyz_tr(2))**2 ) 
            if (val < this%parmtrs(3)) neg = .true.
            
        case(10) !> sph
            xyz_tr(1) = xyz(1) - this%parmtrs(1)
            xyz_tr(2) = xyz(2) - this%parmtrs(2)
            xyz_tr(3) = xyz(3) - this%parmtrs(3)
            
            val = sqrt((xyz_tr(1))**2 + (xyz_tr(2))**2 +(xyz_tr(3))**2 ) 
            if (val < this%parmtrs(4)) neg = .true.
            
        end select
        
        
    end function 

    !> distance to surface boundary
    function surf_gp(A,B,C,D,xyz,uvw) result(dist) 
        real(8),intent(in) :: A,B,C,D, xyz(3), uvw(3)
        real(8) :: dist
        
        dist = (D-(A*xyz(1)+B*xyz(2)+C*xyz(3)))/(A*uvw(1)+B*uvw(2)+C*uvw(3))
        
    end function  
    
    function surf_px(surf,xyz,uvw) result(dist)
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: dist 
        
        if (surf%surf_type /= px) print *, "ERROR : WRONG SURFACE" 
        
        if (uvw(1) == 0) then 
            dist = INFINITY
        else 
            dist = (surf%parmtrs(1)-xyz(1))/uvw(1)
        endif
        
        if (dist<0) dist = INFINITY
        
    end function
    function surf_py(surf,xyz,uvw) result(dist)
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: dist 
        
        if (surf%surf_type /= py) print *, "ERROR : WRONG SURFACE" 
        
        if (uvw(2) == 0) then 
            dist = INFINITY
        else 
            dist = (surf%parmtrs(1)-xyz(2))/uvw(2)
        endif
        
        if (dist<0) dist = INFINITY
        
    end function
    function surf_pz(surf,xyz,uvw) result(dist)
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: dist 
        
        if (surf%surf_type /= pz) print *, "ERROR : WRONG SURFACE" 
        
        if (uvw(3) == 0) then 
            dist = INFINITY
        else 
            dist = (surf%parmtrs(1)-xyz(3))/uvw(3)
        endif
        
        if (dist<0) dist = INFINITY
        
    end function
    
    
    
    function surf_cylz(surf,xyz,uvw) result(dist)
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: dist 
        real(8) :: a, k, c, xyz_(3), val 
        integer :: i 
        
        if (surf%surf_type /= cylz) print *, "ERROR : WRONG SURFACE" 
        xyz_(1:2) = xyz(1:2) - surf%parmtrs(1:2)
        
        a = uvw(1)**2 + uvw(2)**2
        c = xyz_(1)**2 + xyz_(2)**2 - surf%parmtrs(3)**2
        k = xyz_(1)*uvw(1) + xyz_(2)*uvw(2)
        val = k**2 - a*c
        
        if ((a == 0).or.(val < 0)) then 
            dist = INFINITY
        elseif (c < 0) then
            dist = (-k + sqrt(val))/a
            
        elseif (c > 0) then
            dist = (-k - sqrt(val))/a
        endif
        
        if (dist<0) dist = INFINITY
        
    end function
    
    function surf_sqcz(surf,xyz,uvw) result(dist) 
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: d(4), xyz_(3), temp, r
        real(8) :: dist 
        integer :: i
        
        if (surf%surf_type /= sqcz) print *, "ERROR : WRONG SURFACE" 
        
        xyz_(1:2) = xyz(1:2) - surf%parmtrs(1:2) 
        xyz_(3)   = xyz(3)  
        
        r = surf%parmtrs(3) 
        d(1) = (-r-xyz_(2))/uvw(2)
        d(2) = (-r-xyz_(1))/uvw(1)
        d(3) = ( r-xyz_(2))/uvw(2)
        d(4) = ( r-xyz_(1))/uvw(1)
        
        temp = xyz_(1)+d(1)*uvw(1)
        if ((temp < -r).or.(temp > r)) d(1) = INFINITY
        
        temp = xyz_(2)+d(2)*uvw(2)
        if ((temp < -r).or.(temp > r)) d(2) = INFINITY
        
        temp = xyz_(1)+d(3)*uvw(1)
        if ((temp < -r).or.(temp > r)) d(3) = INFINITY
        
        temp = xyz_(2)+d(4)*uvw(2)
        if ((temp < -r).or.(temp > r)) d(4) = INFINITY
        
        do i = 1, 4 
            if (d(i) < 0) d(i) = INFINITY 
        enddo 
        
        
        dist = minval(d(:))
        
    end function
    
    function surf_sph(surf,xyz,uvw) result(dist) 
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: xyz_(3), k, c
        real(8) :: dist, temp
        integer :: i 
        
        if (surf%surf_type /= sph) print *, "ERROR : WRONG SURFACE" 
        
        k = 0; c = 0 
        do i = 1, 3 
            xyz_(i) = xyz(i) - surf%parmtrs(i) 
            k = k + xyz_(i) * uvw(i)
            c = c + xyz_(i)**2 
        enddo 
        c = c - surf%parmtrs(4)**2
        temp = k**2 - c
        if (temp < 0) then
            dist = INFINITY
        elseif (c < 0) then 
            dist = -k + sqrt(temp)
        elseif (c > 0) then 
            dist = -k - sqrt(temp)
        endif 
        
        if (dist < 0) dist = INFINITY 
        
        
    end function
    
    
    subroutine surf_select(surf,xyz,uvw, dist) 
        type(surface) :: surf
        real(8), intent(in) :: xyz(3), uvw(3)
        real(8) :: dist
        
        
        if (surf%surf_type == px)    dist = surf_px(surf,xyz,uvw)
        if (surf%surf_type == py)    dist = surf_py(surf,xyz,uvw)
        if (surf%surf_type == pz)   dist = surf_pz(surf,xyz,uvw)
        if (surf%surf_type == cylz) dist = surf_cylz(surf,xyz,uvw)
        if (surf%surf_type == sqcz) dist = surf_sqcz(surf,xyz,uvw)
        if (surf%surf_type == sph)  dist = surf_sph(surf,xyz,uvw)
        
    end subroutine
    
    
end module 



