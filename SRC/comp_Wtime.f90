!!!!!!!!!!!
  function comp_Wtime()

    implicit none
	
    double precision :: comp_Wtime
    integer :: time0(8)
    call date_and_time(values=time0)
    comp_Wtime = time0(5) * 3600 + time0(6)*60 + time0(7) +0.001 * time0(8)
    
  end function comp_Wtime


