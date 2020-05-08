	  module Main
		  use Geo 
		  use LBM
		  use Fluidflow
		  implicit none
		  !------------------  LB parameters
		  Real (RK):: stime,rtime
		  real (RK):: delta_x 
		  integer  :: iteration
		  !------------------ Geometry
		  type(Latice_D3Q15) :: lattice 
		  !------------------ distribution functions
		  type(DistroFunction) ::f
		  real (RK), allocatable :: f_now(:,:,:,:),f_next(:,:,:,:),f_eq(:,:,:,:)  
		  real (RK):: f_delta_t,f_tau   
		  Character (len=30) :: f_concept 
          type(Wall_info), allocatable :: F_WallBC(:,:,:)


		contains
			subroutine First_Calc
			implicit none
        
			open(11,file='Input_Parameters.Input')	
				read(11,*) ! nx
				read(11,*) nx
				read(11,*) ! ny
				read(11,*) ny
				read(11,*) ! nz
				read(11,*) nz                    
				read(11,*) ! Distribution Concept
				read(11,*) f_concept                       
				read(11,*) ! Tau	
				read(11,*) f_tau
			close(11)

			delta_x=1.0_RK
			f_delta_t=1.0_RK


			Print*,' Calculation Results for current program:'
			print*,' Delta x: ',delta_x
			print*,' Delta t: ',f_delta_t    
			print*;
		end subroutine
      end module Main