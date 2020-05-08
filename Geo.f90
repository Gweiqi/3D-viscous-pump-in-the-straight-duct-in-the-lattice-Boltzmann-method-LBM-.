    Module Geo
    use DataType   
	implicit none
	integer:: nx,ny,nz
	integer:: Wallnum,Fluidnum,uselessnum
	logical, allocatable :: solid(:,:,:),wall(:,:,:),useless(:,:,:)


	interface ini_lattice
          module procedure ini_lattice_D3Q15
    end interface  

	interface Geo_maker
          module procedure Geo_maker_3D
    end interface 

	interface Geo_Analysis
		module procedure Geo_Analysis_D3Q15
	end interface

    contains
        !----------------------------------
        subroutine ini_lattice_D3Q15(lattice_in)
        implicit none
        type(Latice_D3Q15),intent(inout):: lattice_in
		    lattice_in%v_x(0)=0.0_RK; lattice_in%v_x(1)=1.0_RK; lattice_in%v_x(2)=0.0_RK; lattice_in%v_x(3)=-1.0_RK; lattice_in%v_x(4)=0.0_RK
		    lattice_in%v_y(0)=0.0_RK; lattice_in%v_y(1)=0.0_RK; lattice_in%v_y(2)=1.0_RK; lattice_in%v_y(3)=0.0_RK ; lattice_in%v_y(4)=-1.0_RK
		    lattice_in%v_x(5)=0.0_RK; lattice_in%v_x(6)=0.0_RK; lattice_in%v_x(7)=1.0_RK; lattice_in%v_x(8)=1.0_RK
		    lattice_in%v_y(5)=0.0_RK; lattice_in%v_y(6)=0.0_RK ; lattice_in%v_y(7)=1.0_RK; lattice_in%v_y(8)=1.0_RK
            lattice_in%v_x(9)=1.0_RK; lattice_in%v_x(10)=1.0_RK; lattice_in%v_x(11)=-1.0_RK
            lattice_in%v_x(12)=-1.0_RK;lattice_in%v_x(13)=-1.0_RK;lattice_in%v_x(14)=-1.0_RK
            lattice_in%v_y(9)=-1.0_RK;lattice_in%v_x(10)=-1.0_RK;lattice_in%v_x(11)=1.0_RK
            lattice_in%v_x(12)=1.0_RK;lattice_in%v_x(13)=-1.0_RK;lattice_in%v_x(14)=-1.0_RK	
            lattice_in%v_z(0)=0.0_RK;lattice_in%v_z(1)=0.0_RK;lattice_in%v_z(2)=0.0_RK;lattice_in%v_z(3)=0.0_RK;lattice_in%v_z(4)=0.0_RK
            lattice_in%v_z(5)=1.0_RK;lattice_in%v_z(6)=-1.0_RK;lattice_in%v_z(7)=1.0_RK;lattice_in%v_z(8)=-1.0_RK;lattice_in%v_z(9)=-1.0_RK
            lattice_in%v_z(10)=1.0_RK;lattice_in%v_z(11)=-1.0_RK;lattice_in%v_z(12)=1.0_RK;lattice_in%v_z(13)=1.0_RK;lattice_in%v_z(14)=-1.0_RK
            !--------------------
		    lattice_in%inv_v(0)=0; lattice_in%inv_v(1)=3; lattice_in%inv_v(2)=4; lattice_in%inv_v(3)=1; lattice_in%inv_v(4)=2
		    lattice_in%inv_v(5)=6; lattice_in%inv_v(6)=5; lattice_in%inv_v(7)=14; lattice_in%inv_v(14)=7; lattice_in%inv_v(8)=13;lattice_in%inv_v(13)=8
            lattice_in%inv_v(9)=12;lattice_in%inv_v(12)=9;lattice_in%inv_v(10)=11;lattice_in%inv_v(11)=10
            !---------------
		    lattice_in%w_grid(0)=16.0_RK/72.0_RK
		    lattice_in%w_grid(1)=8.0_RK/72.0_RK;  lattice_in%w_grid(2)=8.0_RK/72.0_RK;  lattice_in%w_grid(3)=8.0_RK/72.0_RK;  lattice_in%w_grid(4)=8.0_RK/72.0_RK
		    lattice_in%w_grid(5)=8.0_RK/72.0_RK; lattice_in%w_grid(6)=8.0_RK/72.0_RK; lattice_in%w_grid(7)=1.0_RK/72.0_RK; lattice_in%w_grid(8)=1.0_RK/72.0_RK
            lattice_in%w_grid(9)=1.0_RK/72.0_RK;lattice_in%w_grid(10)=1.0_RK/72.0_RK;lattice_in%w_grid(11)=1.0_RK/72.0_RK
            lattice_in%w_grid(12)=1.0_RK/72.0_RK;lattice_in%w_grid(13)=1.0_RK/72.0_RK;lattice_in%w_grid(14)=1.0_RK/72.0_RK
            print*,' Lattice parameters are set for',lattice_in%Q+1,'directions'  
            print*;
        end subroutine     
        !----------------------------------
        subroutine ini_Geo
              implicit none
			  !-- Local
              integer :: AllocateStatus


			  if (nx==0 .or. ny==0 .or. nz==0) Stop (' Wrong Dimensions')

              allocate (solid(0:nx-1,0:ny-1,0:nz-1), STAT = AllocateStatus)
              IF (AllocateStatus /= 0) STOP "*** Not enough memory - Solid ***"
              solid=.false.
              allocate (wall(0:nx-1,0:ny-1,0:nz-1), STAT = AllocateStatus)
              IF (AllocateStatus /= 0) STOP "*** Not enough memory - Wall ***"
              wall=.false.
              allocate (useless(0:nx-1,0:ny-1,0:nz-1), STAT = AllocateStatus)
              IF (AllocateStatus /= 0) STOP "*** Not enough memory - Useless ***"
              useless=.false.

              Print*,'Solid, Wall, and Useless are initialaized'
        end subroutine     
        !--------------------------
        Subroutine Geo_maker_3D(lattice_in)
        implicit none
		!-- IN
        type(Latice_D3Q15),intent(inout):: lattice_in
		!-- Local
		integer:: x,y,z,geo_ctr
        geo_ctr=100
        do while (Geo_ctr>0) 
             print*, ' Geo Maker ...'
			 print*;
             Print*,' What do you want?   0 - STOP    1-Wall   2-Rectangle  3-Circle'
             read*,geo_ctr
                select case(geo_ctr)
                case default
                    Print*, 'Wrong choice! - geo_ctr'
                    stop
                case(0)
                    print*, ' It is Done!'
                case(1)
                    call Wall_3D
                case(2)
                    call Rectangle
                case(3)
                    call Circle
                end select
                call Bound_Chk_3D(Solid,Solid,Solid,lattice_in)
            end do 
        end subroutine        
        !-----------------
	    subroutine select_zone_3D(x_ld,y_ld,z_ld,x_ru,y_ru,z_ru,n_x,n_y,n_z)
		    implicit none
			!-- INOUT
			real (RK), intent(INOUT):: x_ld,x_ru,y_ld,y_ru,z_ld,z_ru
			!-- IN
			integer, intent(in):: n_x,n_y,n_z
			!-- Local
		    real (RK) a,b
		    integer zone_ctr
            
		    print*;
		    print*,' Select active zone: All(0) or Specified Zone(1)'
		    read*,zone_ctr
		    select case(zone_ctr)
			    case(1)
				    print*;
				    Print*,' Ref. L is Ny=',n_y,', Now input your value corresponding to this value';print*;
				    print*,' X-left & down = a/b*L : a,b?'
				    read*,a,b
				    x_ld=a*real(n_y-1)/b
				    print*,' Y-left & down = a/b*L : a,b?'
				    read*,a,b
				    y_ld=a*real(n_y-1)/b
                    print*,' Z-Front  = a/b*L : a,b?'
				    read*,a,b
				    z_ld=a*real(n_y-1)/b
				    print*,' X-Right & Up  = a/b*L : a,b?'
				    read*,a,b
				    x_ru=a*real(n_y-1)/b
				    print*,' Y-Right & Up  = a/b*L : a,b?'
				    read*,a,b
				    y_ru=a*real(n_y-1)/b
                    print*,' Z-Back  = a/b*L : a,b?'
				    read*,a,b
				    z_ru=a*real(n_y-1)/b
                    
			    case(0)
				    x_ld=0.0
				    y_ld=0.0
                    z_ld=0.0
				    x_ru=real(n_x-1)
				    y_ru=real(n_y-1)
                    z_ru=real(n_z-1)
                end select
                if (x_ld<0) x_ld=0.0
                if (y_ld<0) y_ld=0.0
                if (z_ld<0) z_ld=0.0  
                if (x_ru>real(n_x-1)) x_ru=real(n_x-1)
                if (y_ru>real(n_y-1)) y_ru=real(n_y-1)
                if (z_ru>real(n_z-1)) z_ru=real(n_z-1)                 
	    return
        end subroutine
        !---------------
        subroutine Rectangle
        implicit none
        integer:: fill_ctr,x,y,z
        real (RK):: a,b,x_ld,x_ru,y_ld,y_ru,z_ld,z_ru
        
        print*;
        print*,'                 Clear(0) or Fill (1) : ?'
        read*,fill_ctr
        print*,' Accurate corner cordinates :';print*;
        call select_zone_3D(x_ld,y_ld,z_ld,x_ru,y_ru,z_ru,nx,ny,nz)
        select case (fill_ctr)
            case (1)
                do x=max(0,int(x_ld)-1),min(int(x_ru)+1,nx-1)
                    do y=max(0,int(y_ld)-1),min(int(y_ru)+1,ny-1)
                      do z=max(0,int(z_ld)-1),min(int(z_ru)+1,nz-1)
                        if ((x>=x_ld) .and. (x<=x_ru) .and. (y>=y_ld) .and. (y<=y_ru) .and. (z>=z_ld) .and. (z<=z_ru)) then
                            solid(x,y,z)=.true.
                        end if
                    end do
                end do
             end do
            case (0)
                do x=max(0,int(x_ld)-1),min(int(x_ru)+1,nx-1)
                    do y=max(0,int(y_ld)-1),min(int(y_ru)+1,ny-1)
                       do z=max(0,int(z_ld)-1),min(int(z_ru)+1,nz-1)
                        if ((x>x_ld) .and. (x<x_ru) .and. (y>y_ld) .and. (y<y_ru) .and. (z>z_ld) .and. (z<z_ru)) then
                            solid(x,y,z)=.false.
                        end if
                    end do
                end do
             end do
            case default
				stop (' Wrong choice!')               
        end select
        end subroutine
        !-----------------------
        	subroutine Circle
        implicit none
        integer:: fill_ctr,x,y,z
        real (RK):: a,b,x_c,y_c,R_c
        
        print*;
        print*,'                 Clear(0) or Fill (1) : ?'
        read*,fill_ctr
        print*,' Accurate corner cordinates :';print*;
		print*,' X-Center = a/b*L : a,b?'
		read*,a,b
		x_c=a*real(ny-1)/b
		print*,' Y-Center = a/b*L : a,b?'
		read*,a,b
		y_c=a*real(ny-1)/b
		print*,' Radius = a/b*L : a,b?'
		read*,a,b
		R_c=a*real(ny-1)/b
        select case (fill_ctr)
            case (1)
                do x=0,nx-1
                    do y=0,ny-1
                       do z=0,nz-1
                        if ( (real(x)-x_c)**2+(real(y)-y_c)**2<=R_c**2 ) then
                            solid(x,y,z)=.true.
                        end if
                    end do
                end do
              end do
            case (0)
                do x=0,nx-1
                    do y=0,ny-1
                       do z=0,nz-1
                        if ( (real(x)-x_c)**2+(real(y)-y_c)**2<R_c**2 ) then
                            solid(x,y,z)=.false.
                        end if
                    end do
                end do
             end do
            case default
				stop (' Wrong choice!')               
        end select
        end subroutine
        !------------
	    subroutine Bound_Chk_3D(Solid_in,Wall_in,Useless_in,lattice_in)
		implicit none
		! --- Local
		integer:: S_val,W_val,Us_val,x,y,z
		! --- IN
		type(Latice_D3Q15), intent (in) :: lattice_in
        logical, intent (IN):: Solid_in(0:nx-1,0:ny-1,0:nz-1),Wall_in(0:nx-1,0:ny-1,0:nz-1),useless_in(0:nx-1,0:ny-1,0:nz-1)
        
		open(10,file='Bound_Chk.dat')
		write(10,*) 'VARIABLES = X, Y, Z, Solid, Wall, Useless' 
		write(10,*) 'ZONE I=', nx, ', J=', ny, ', K=', nz, ',F=POINT'
     do z=0,nz-1
		do y=0,ny-1
			do x=0,nx-1
				if (solid_in(x,y,z)) then
					s_val=1
				else
					s_val=0
                end if
				if (Wall_in(x,y,z)) then
					W_val=1
				else
					W_val=0
                end if
				if (useless_in(x,y,z)) then
					Us_val=1
				else
					Us_val=0
				end if                
				write(10,*) x, y, z, S_val, W_val, Us_val
			end do
		end do
      end do
		close(10)
		print*;
		print*,' RE : Bound_Chk.Dat is updated'
        end subroutine        
        !------------
        subroutine Wall_3D
        implicit none
        integer:: Wall_ctr,x,y,z

        Wall_ctr=10
        do while (Wall_ctr>0)
            print*;
            print*,' Where?  Right(6) Up(8) Left(4) Down(2) Front(1) Back(3)  Stop(0)'
            read*,Wall_ctr
			select case(Wall_ctr)
            case(6)    ! right
                    x=nx-1
					do y=0,ny-1
                       do z=0,nz-1
						solid(x,y,z)=.true.
					end do
                  end do
            case(8)    ! up
					y=ny-1
					do x=0,nx-1
                       do z=0,nz-1
						solid(x,y,z)=.true.
					end do
                  end do
            case(4)    ! left
					x=0
					do y=0,ny-1
                       do z=0,nz-1
						solid(x,y,z)=.true.
					end do
                  end do
            case(2)    ! down
					y=0
					do x=0,nx-1
                       do z=0,nz-1
						solid(x,y,z)=.true.
                    end do
                  end do
             case(1)    ! Front
					z=0
					do x=0,nx-1
                       do y=0,ny-1
						solid(x,y,z)=.true.
                    end do
                  end do 
             case(3)    ! Back
					z=nz-1
					do x=0,nx-1
                       do y=0,ny-1
						solid(x,y,z)=.true.
                    end do
                  end do        	
            case(0)
                print*,' It is Done!'
            case default
                stop (' Wrong choice!')  
            end select
            print*;
		end do
        end subroutine   
        !---------------
        subroutine Geo_Analysis_D3Q15(lattice_in)
        implicit none
		! --- In
		type(Latice_D3Q15), intent (in) :: lattice_in
		! -- Local
        integer::x,y,z,xp,yp,zp,xm,ym,zm

        wallnum=0
        fluidnum=0
        Uselessnum=0
        do x=0,nx-1
	        xp=min(x+1,nx-1)
	        xm=max(x-1,0)
	        do y=0,ny-1
		        yp=min(y+1,ny-1)
		        ym=max(y-1,0)
                  do z=0,nz-1
                    zp=min(z+1,nz-1)
		            zm=max(z-1,0)
		        if (solid(x,y,z)) then
			        if ( .not. (solid(xm,y,z) .and. solid(xp,y,z)  .and. solid(xm,yp,z) .and. solid(xp,yp,z) .and. &
				               solid(xm,ym,z) .and. solid(xp,ym,z) .and.  solid(x,yp,z) .and. solid(x,ym,z))) then
                             !  solid(x,ym,zm) .and. solid(x,ym,zp) .and. solid(x,yp,zm) .and. solid(x,yp,zp) .and. &
                             !  solid(x,y,zm) .and. solid(x,y,zp) .and. & 
                             !  solid(xm,y,zm) .and. solid(xm,y,zp) .and. solid(xp,y,zm) .and. solid(xp,y,zp)) then
				        wallnum=wallnum+1
				        wall(x,y,z)=.true.
			        else
				        Uselessnum=Uselessnum+1
				        Useless(x,y,z)=.true.	
			        end if
		        else
			        fluidnum=fluidnum+1
		        end if
	        end do
        end do
      end do
        call Bound_Chk_3D(Solid,Wall,Useless,lattice_in)
        print*;
        print*,' Total Useless nodes: ',uselessnum			 		
        print*,' Total fluid   nodes: ',fluidnum
        print*,' Total wall    nodes: ',wallnum  
        print*;
        end subroutine
        !---------------        
    End Module Geo