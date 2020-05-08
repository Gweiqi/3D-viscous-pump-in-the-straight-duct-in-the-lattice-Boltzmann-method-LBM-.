

      Include 'Datatype.f90' 
      Include 'Geo.f90' 
      Include 'LBM.f90'
      Include 'Fluidflow.f90'
	  Include 'currentproblem.f90'
	  !===========================================================


	  !=================================================================================	     
      program prog
      use main
      implicit none


      call First_Calc
      !    must use for every simulation
      call ini_lattice(lattice)
      call ini_Geo
      Call Geo_maker(lattice)
      Call Geo_Analysis(lattice)
      !   must use for every distribution
      call distribution_allocation(f,f_now,f_next,f_eq,delta_x,f_delta_t,f_tau,f_concept,lattice%Q,f_WallBC)
      call distribution_display(f) 
      call distribution_initialaize(f_now,f_next,lattice,1.0_RK)

      !--------------------------------- Ready to go
      call cpu_time(stime)
      iteration=0
      do while (5000>iteration)
		    call LB_Stream(f_now,f_next,lattice)
		 !   call Bounce_Back(f_now,f_next,lattice)
           	call Flow_BC(f_next,lattice,f,f_WallBC)  
		    call F_Eq_calc(lattice,f,f_eq,f_next)          
		    call LB_Collision(f_now,f_next,f_eq,lattice%Q,f)
          !  call accelerate(lattice,f_now)
		    if (mod(iteration,100)==0) Print*, iteration
            iteration=iteration+1
      end do 
	  call cpu_time(rtime) 
	  print*,' Runtime: ',rtime-stime
	  call Flow_output(lattice,f,f_next)
          
    end program prog
    
	!================================
