    module Fluidflow
		use DataType
		use Geo
		implicit none

		interface Flow_output
			module procedure Flow_output_3D
		end interface 
    
		interface F_eq_calc
			module procedure F_eq_calc_D3Q15
        end interface
        
     !   interface accelerate
        !    module procedure accelerate_D3Q15
     !   end interface
    interface Flow_BC
			module procedure Flow_BC_D3Q15
		end interface    
    
    
    contains

        !****************************    
    subroutine Flow_wallBC(distribution_WallBC)
		implicit none
		type(Wall_info), intent(INOUT)::distribution_WallBC(0:nx-1,0:ny-1,0:nz-1)
		integer::x,y,z
		Real(RK)::x_c,y_c,teta,u_r
		do x=0,nx-1
			do y=0,ny-1
				do z=0,nz-1
					distribution_WallBC(x,y,z)%ctr=0
					distribution_WallBC(x,y,z)%VAL1=0
					distribution_WallBC(x,y,z)%VAL2=0
					distribution_WallBC(x,y,z)%VAL3=0
				end do
			end do
		end do
		!--
		x_c=real(nx-1)/2.0
		y_c=real(ny-1)/3.0
		u_r=0.05
		!z=0
		do x=5,nx-5
		do y=5,ny-5
        do z=5,nz-5
			if (wall(x,y,z)) then
				teta=atan2(real(y)-y_c,real(x)-x_c)
				distribution_WallBC(x,y,z)%ctr=0
				distribution_WallBC(x,y,z)%VAL1=u_r*sin(teta)
				distribution_WallBC(x,y,z)%VAL2=-u_r*cos(teta)
                distribution_WallBC(x,y,z)%VAL3=0
			end if 
		end do
		end do
        end do
		!--
		x=nx-1;z=0
		do y=1,ny-2
			distribution_WallBC(x,y,z)%ctr=1
		end do
		!--	
		x=0;z=0
		do y=1,ny-2
			distribution_WallBC(x,y,z)%ctr=1
		end do
		!--			
		return
		end subroutine
    !---------------------
        subroutine local_flow_D3Q15(x,y,z,rho,u,lattice_in,info,distribution_next)
	    implicit none
        !--- In
        type(Latice_D3Q15), intent(IN) :: lattice_in
        type(DistroFunction), Intent(IN) ::info
	    integer,intent(in):: x,y,z
        real (RK), intent(IN) :: distribution_next(0:nx-1,0:ny-1,0:nz-1,0:lattice_in%Q)
        !--- IN/OUT
	    real (RK), intent(inout):: rho,u(0:2)        
        !--- local
        integer:: k              
        
        rho=0.d0;	    u=0.d0
	    do k=0,lattice_in%Q
            rho=rho+distribution_next(x,y,z,k)
		    u(0)=u(0)+distribution_next(x,y,z,k)*lattice_in%v_x(k)*info%C
		    u(1)=u(1)+distribution_next(x,y,z,k)*lattice_in%v_y(k)*info%C
            u(2)=u(2)+distribution_next(x,y,z,k)*lattice_in%v_z(k)*info%C
	    end do
	    u(0)=u(0)/rho
	    u(1)=u(1)/rho
        u(2)=u(2)/rho
	    return
        end subroutine 
        !****************************
	    subroutine Flow_output_3D(lattice_in,info,distribution_next)
        implicit none 
        !---- In
        type(Latice_D3Q15), intent(IN) :: lattice_in
        type(DistroFunction), Intent(IN) ::info
		real (RK), intent(IN) :: distribution_next(0:nx-1,0:ny-1,0:nz-1,0:lattice_in%Q)
        !--- IN/OUT
        !--- local
	    integer x,y,z
	    real(RK) u(0:2),rho,press,rhoave
	    character (len=300) :: filename
        
        rhoave=1.0
	    filename='FlowData_3D.dat'
	    open(2,file=filename)
	    write(2,*) 'VARIABLES = X, Y, Z, VX, VY, VZ, PRESS' 
	    write(2,*) 'ZONE I=', nx, ', J=', ny, ', K=', nz, ', F=POINT'
     do z=0,nz-1
	    do y=0,ny-1
		    do x=0,nx-1
                if (useless(x,y,z)) then
                    u=0.0
                    rho=1.0
                else
		            rho=0.0;u=0.0
		            call local_flow_D3Q15(x,y,z,rho,u,lattice_in,info,distribution_next)
                end if
			    press=(rho-rhoave)*info%Cs2
			    write(2,*) real(x), real(y), real(z), u(0), u(1), u(2), press
		    end do
	    end do
     end do
	    close(2)
	    return
        end subroutine   
        !*******************
        Subroutine F_eq_calc_D3Q15(lattice_in,info,eq_distribution,distribution_next)
        implicit none
        ! --- In
        type(Latice_D3Q15),intent(In) :: lattice_in 
        type(DistroFunction), Intent(IN) ::info  
        real (RK), intent(IN) :: distribution_next(0:nx-1,0:ny-1,0:nz-1,0:lattice_in%Q)
        !--- IN/OUT
        real (RK), intent(INOUT) :: eq_distribution(0:nx-1,0:ny-1,0:nz-1,0:lattice_in%Q) 
        ! --- local       
        real (RK):: rho,u(0:2),p
        integer  :: x,y,z,k
        Real (RK):: fact1,fact2,fact3
        
        fact1=3.0_RK/info%C2
        fact2=4.5_RK/info%C4
        fact3=1.5_Rk/info%C2

        !print*,lattice_in%Q
        !pause ('feq')
        
        do x=0,nx-1
			do y=0,ny-1
               do z=0,nz-1 
				if (.not. solid(x,y,0)) then
					rho=0.0;    u=0.0
					do k=0,lattice_in%Q
						rho=rho+distribution_next(x,y,z,k)
						u(0)=u(0)+distribution_next(x,y,z,k)*lattice_in%v_x(k)*info%C
						u(1)=u(1)+distribution_next(x,y,z,k)*lattice_in%v_y(k)*info%C
                        u(2)=u(2)+distribution_next(x,y,z,k)*lattice_in%v_z(k)*info%C
					end do
					u(0)=u(0)/rho
					u(1)=u(1)/rho
                    u(2)=u(2)/rho
					do k=0,lattice_in%Q
						p=info%C*lattice_in%v_x(k)*u(0)+info%C*lattice_in%v_y(k)*u(1)+info%C*lattice_in%v_z(k)*u(2)
						eq_distribution(x,y,z,k)=lattice_in%w_grid(k)*rho*(1.0_RK+fact1*p+fact2*p*p-fact3*(u(0)*u(0)+u(1)*u(1)+u(2) *u(2)) )
					end do
				end if
            end do
        end do
     end do
		return   
        end subroutine 
        !************************
  !      Subroutine accelerate_D3Q15(lattice_in,distribution_now)
  !      implicit none
  !      ! --- In
  !      type(Latice_D3Q15),intent(In) :: lattice_in
  !      !--- IN/OUT
  !      real (RK), intent(INOUT) :: distribution_now(0:nx-1,0:ny-1,0:nz-1,0:lattice_in%Q)
  !      ! --- local       
  !      integer  :: x,y,z
  !      do x=0,nx-1
		!	do y=0,ny-1
  !             do z=0,nz-1
		!		if (.not. solid(x,y,z)) then
		!				distribution_now(x,y,z,1)=distribution_now(x,y,z,1)+lattice_in%w_grid(1)*(0.005)
  !                      distribution_now(x,y,z,7)=distribution_now(x,y,z,7)+lattice_in%w_grid(7)*(0.005)
  !                      distribution_now(x,y,z,8)=distribution_now(x,y,z,8)+lattice_in%w_grid(8)*(0.005)
  !                      distribution_now(x,y,z,9)=distribution_now(x,y,z,9)+lattice_in%w_grid(9)*(0.005)
  !                      distribution_now(x,y,z,10)=distribution_now(x,y,z,10)+lattice_in%w_grid(10)*(0.005)
  !                      distribution_now(x,y,z,3)=distribution_now(x,y,z,3)-lattice_in%w_grid(3)*(0.005)
  !                      distribution_now(x,y,z,11)=distribution_now(x,y,z,11)-lattice_in%w_grid(11)*(0.005)
  !                      distribution_now(x,y,z,12)=distribution_now(x,y,z,12)-lattice_in%w_grid(12)*(0.005)
  !                      distribution_now(x,y,z,13)=distribution_now(x,y,z,13)-lattice_in%w_grid(13)*(0.005)
  !                      distribution_now(x,y,z,14)=distribution_now(x,y,z,14)-lattice_in%w_grid(14)*(0.005)
		!		end if
  !          end do
  !      end do
  !   end do
		!return   
  !      end subroutine 
        !*************************
        	subroutine Flow_BC_D3Q15(distribution_next,lattice_in,info,distribution_wallBC)
		implicit none
		!--- In
		type(Latice_D3Q15),intent(in):: lattice_in
		type(DistroFunction), Intent(IN) ::info
		type(Wall_info),intent(IN)::distribution_wallBC(0:nx-1,0:ny-1,0:nz-1)
        !--- InOut
        Real (RK), intent(INOUT):: distribution_next(0:nx-1,0:ny-1,0:nz-1,0:lattice_in%Q)  
        !--- Local
		integer::x,y,z,k,x_k,y_k,z_k,x_invk,y_invk,z_invk
		integer::x_bef,y_bef,z_bef,x_next,y_next,z_next
		logical::find
		logical::xk_out(0:lattice_in%Q),yk_out(0:lattice_in%Q),zk_out(0:lattice_in%Q),xinvk_out(0:lattice_in%Q),yinvk_out(0:lattice_in%Q),zinvk_out(0:lattice_in%Q)
		logical::k_out(0:lattice_in%Q),invk_out(0:lattice_in%Q)
		real(RK)::rho_next,f_neq_next(0:lattice_in%Q),eq_next,ux_next,uy_next,uz_next,p_next
		real(RK)::p_wall,ux_wall,uy_wall,uz_wall,f_eq_wall
		Real(RK):: fact1,fact2,fact3
        
        fact1=3.0_RK/info%C2
        fact2=4.5_RK/info%C4
        fact3=1.5_Rk/info%C2

		do x=0,nx-1
			do y=0,ny-1
               do z=0,nz-1
				if (wall(x,y,z)) then				
						!---------------------------- Next Dir
						find=.false.
						xk_out=.false.;yk_out=.false.;zk_out=.false.
						xinvk_out=.false.;yinvk_out=.false.;zinvk_out=.false.
						k_out=.false.;invk_out=.false.
						do k=1,6
							x_k=mod(x+int(lattice_in%v_x(k))+nx,nx)
							y_k=mod(y+int(lattice_in%v_y(k))+ny,ny)
                            z_k=mod(z+int(lattice_in%v_z(k))+nz,nz)
							if (abs(x_k-x)>1) xk_out(k)=.true.
							if (abs(y_k-y)>1) yk_out(k)=.true.
                            if (abs(z_k-z)>1) zk_out(k)=.true.
							if (xk_out(k) .or. yk_out(k) .or. zk_out(k)) k_out(k)=.true.
                         
							x_invk=mod(x+int(lattice_in%v_x(lattice_in%inv_v(k)))+nx,nx)
							y_invk=mod(y+int(lattice_in%v_y(lattice_in%inv_v(k)))+ny,ny)
                            z_invk=mod(z+int(lattice_in%v_z(lattice_in%inv_v(k)))+nz,nz)
							if (abs(x_invk-x)>1) xinvk_out(k)=.true.
							if (abs(y_invk-y)>1) yinvk_out(k)=.true.
                            if (abs(z_invk-z)>1) zinvk_out(k)=.true.
							if (xinvk_out(k) .or. yinvk_out(k) .or. zinvk_out(k)) invk_out(k)=.true.

							if ( (.not. solid(x_k,y_k,z_k)) .and. (.not. k_out(k)) .and. &
								(useless(x_invk,y_invk,z_invk) .or. invk_out(k)) ) then
								x_next=x_k
								y_next=y_k
                                z_next=z_k
								find=.true.
								exit
							end if
						end do
						if (.not. find) then
							do k=7,14
								x_k=mod(x+int(lattice_in%v_x(k))+nx,nx)
								y_k=mod(y+int(lattice_in%v_y(k))+ny,ny)
                                z_k=mod(z+int(lattice_in%v_z(k))+nz,nz)
								if (abs(x_k-x)>1) xk_out(k)=.true.
								if (abs(y_k-y)>1) yk_out(k)=.true.
                                if (abs(z_k-z)>1) zk_out(k)=.true.
								if (xk_out(k) .or. yk_out(k) .or. zk_out(k)) k_out(k)=.true.
                         
								x_invk=mod(x+int(lattice_in%v_x(lattice_in%inv_v(k)))+nx,nx)
								y_invk=mod(y+int(lattice_in%v_y(lattice_in%inv_v(k)))+ny,ny)
                                z_invk=mod(z+int(lattice_in%v_z(lattice_in%inv_v(k)))+nz,nz)
								if (abs(x_invk-x)>1) xinvk_out(k)=.true.
								if (abs(y_invk-y)>1) yinvk_out(k)=.true.
                                if (abs(z_invk-z)>1) zinvk_out(k)=.true.
								if (xinvk_out(k) .or. yinvk_out(k) .or. zinvk_out(k)) invk_out(k)=.true.

								if ( (.not. solid(x_k,y_k,z_k)) .and. (.not. k_out(k)) .and. &
									(useless(x_invk,y_invk,z_invk) .or. invk_out(k)) ) then
									x_next=x_k
									y_next=y_k
                                    z_next=z_k
									find=.true.
									exit
								end if
							end do                        
						end if
						!if (.not. find) then
						!	print*,'there is a problem in next direction detection!'
						!	print*,x,y
						!	stop
						!end if
					!---------------------------------------- calc non-eq part for the next node
					rho_next=0.0; ux_next=0.0; uy_next=0.0; uz_next=0.0
					do k=0,lattice_in%Q
						rho_next=rho_next+distribution_next(x_next,y_next,z_next,k)	
						ux_next=ux_next+distribution_next(x_next,y_next,z_next,k)*lattice_in%v_x(k)*info%C
						uy_next=uy_next+distribution_next(x_next,y_next,z_next,k)*lattice_in%v_y(k)*info%C
                        uz_next=uz_next+distribution_next(x_next,y_next,z_next,k)*lattice_in%v_z(k)*info%C						
					end do
					ux_next=ux_next/rho_next
					uy_next=uy_next/rho_next
                    uz_next=uz_next/rho_next
					do k=0,lattice_in%Q
						p_next=info%C*lattice_in%v_x(k)*ux_next+info%C*lattice_in%v_y(k)*uy_next+info%C*lattice_in%v_z(k)*uz_next
						eq_next=lattice_in%w_grid(k)*rho_next*(1.0_RK+fact1*p_next+fact2*p_next*p_next-fact3*(ux_next*ux_next+uy_next*uy_next+uz_next*uz_next))
						f_neq_next(k)=distribution_next(x_next,y_next,z_next,k)-eq_next
					end do
					!-------------------------- find unknown directions
					do k=0,lattice_in%Q
						x_bef=mod(x+int(lattice_in%v_x(lattice_in%inv_v(k)))+nx,nx)
						y_bef=mod(y+int(lattice_in%v_y(lattice_in%inv_v(k)))+ny,ny)
                        z_bef=mod(y+int(lattice_in%v_z(lattice_in%inv_v(k)))+nz,nz)
						if (useless(x_bef,y_bef,z_bef) .or. (abs(x_bef-x)>2) .or. (abs(y_bef-y)>2) .or. (abs(z_bef-z)>2)) then
							!----- calc eq part
							if (distribution_Wallbc(x,y,z)%ctr==0) then
								ux_wall=distribution_Wallbc(x,y,z)%Val1
								uy_wall=distribution_Wallbc(x,y,z)%Val2	
                                uz_wall=distribution_Wallbc(x,y,z)%Val3							
							else
								ux_wall=ux_next
								uy_wall=uy_next
                                uz_wall=uz_next
							end if
							p_wall=ux_wall*info%C*lattice_in%v_x(k)+uy_wall*info%C*lattice_in%v_y(k)+uz_wall*info%C*lattice_in%v_z(k)
							f_eq_wall=lattice_in%w_grid(k)*rho_next*(1.0_RK+fact1*p_wall+fact2*p_wall*p_wall-fact3*(ux_wall*ux_wall+uy_wall*uy_wall+uz_wall*uz_wall))
							!----------
							distribution_next(x,y,z,k)=f_eq_wall+f_neq_next(k)
						end if
					end do				 
				end if
			end do
		end do
     end do
		return
		end subroutine
        !-----------------
    end module Fluidflow