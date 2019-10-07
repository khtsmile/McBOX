subroutine process_line (line, ids, output ) 
    
    implicit none 
    
    character(*), intent(in) :: line 
    integer :: i, ids, j, ide 
    character(30) :: output 
    
    j = ids 
    do while (.true.) 
        j = j+1 
        if (line(j:j).eq.' ') then 
            ide = j-1
            goto 10 
        endif 
    enddo 
    
10    output = line(ids:ide) 

    i = 0 
    do
        if (line(j+i:j+i).eq.' ' .and. (j+i) < len(line)) then 
            i = i+1
        else
            exit
        endif 
    enddo 
    ids = j-1+i

end subroutine 
