!-----------------------------------------------------------------------------!
!WENO Solver for inviscid Burger equation (one-dimensional)
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: July 2015
!-----------------------------------------------------------------------------!

program burger
implicit none
integer::nx,ns,nt
real*8 ::dt,tm,dx,ds,x,t
integer::i,k
real*8,allocatable::u(:,:)


!reading input file
open(7,file='input_sol.txt')
read(7,*)nx 	!number of intervals in x
read(7,*)ns 	!number of record in t
read(7,*)dt     !time step 
read(7,*)tm	    !maximum time 
close(7)

!spatial grid size
dx = 1.0d0/dfloat(nx)

!time interval for storing data
ds = tm/dfloat(ns)

!max number of time for numerical solution
nt = nint(tm/dt)

!numerical solution with RK3 + CRWENO schemes
allocate(u(1:nx,0:ns))
call numerical(nx,ns,nt,dx,dt,u)

!write solutions in Tecplot format

!write time variable solutions
open(3, file="burger_sol_time.plt")
write(3,*) 'variables ="x","u"'
  
do k=0,ns
  	t = dfloat(k)*ds
	write(3,'(a16,i8,a10,f10.4,a3)')'zone f=point i=',nx,',t="time',t,'"'
		
	do i=1,nx
		x = -0.5d0*dx + dfloat(i)*dx
        write(3,*) x, u(i,k)
    end do

end do
close(3)

!write 3D view solutions
open(4, file="burger_sol_3Dtx.plt")
write(4,*) 'variables ="t","x","u"'
write(4,*)'zone f=point i=',ns+1,',j=',nx

do i=1,nx
  
 		x = -0.5d0*dx + dfloat(i)*dx
        
		do k=0,ns
  		t = dfloat(k)*ds		

        write(4,*) t, x, u(i,k)
        end do

end do
close(4)

open(5, file="burger_sol_3Dxt.plt")
write(5,*) 'variables ="x","t","u"'
write(5,*)'zone f=point i=',nx,',j=',ns+1
     
do k=0,ns
  	t = dfloat(k)*ds
    
		do i=1,nx
 		x = -0.5d0*dx + dfloat(i)*dx
    
        write(5,*) x, t, u(i,k)
        end do

end do
close(5)


open(7, file="burger_sol_final.plt")
write(7,*) 'variables ="x","u"'
do i=1,nx
x = -0.5d0*dx + dfloat(i)*dx
write(7,*) x,u(i,ns)
end do
        
end


!-----------------------------------------------------------------------------!
!compute numerical solutions
!	* 3rd-order Runge-Kutta for temporal 
!	* 5th-order WENO Scheme for spatial
!-----------------------------------------------------------------------------!
subroutine numerical(nx,ns,nt,dx,dt,u)
implicit none
integer::nx,ns,nt,i,j,k,freq
real*8 ::dx,dt,x,pi
real*8 ::un(1:nx),ut(1:nx),r(1:nx)
real*8 ::u(1:nx,0:ns)

!u : stored solution 
!un: numerical solution at time n


!initial condition 
!a simple sin wave; developing shock at the center
pi = 4.0d0*datan(1.0d0)

do i=1,nx
x = -0.5d0*dx + dfloat(i)*dx
un(i) = dsin(2.0d0*pi*x)
end do
    

!initial time recording
do i=1,nx
u(i,0)=un(i)
end do


!time integration (with RK3)
k=0 !record index
freq = int(nt/ns) !record frequency
    
do j=1,nt

print*, j

	call rhs(nx,dx,un,r)
    
	do i=1,nx
    ut(i) = un(i) + dt*r(i)
    end do

	call rhs(nx,dx,ut,r)

	do i=1,nx
    ut(i) = 0.75d0*un(i) +0.25d0*ut(i) + 0.25d0*dt*r(i)
    end do

	call rhs(nx,dx,ut,r)

	do i=1,nx
    un(i) = 1.0d0/3.0d0*un(i) +2.0d0/3.0d0*ut(i) + 2.0d0/3.0d0*dt*r(i)
    end do

	!record data
    if(mod(j,freq).eq.0) then
    k=k+1
    do i=1,nx
    u(i,k)=un(i)
    end do
    end if
        
end do


return
end

!-----------------------------------------------------------------------------!
!Compute rhs for numerical solution of inviscid Burger equation
!  r = -u*u' 
!
!we are using conservative flux form: r = -(u*u/2)'
!-----------------------------------------------------------------------------!
subroutine rhs(nx,dx,u,r)
implicit none
integer::nx,i
real*8 ::dx
real*8 ::u(1:nx),r(1:nx)
real*8,allocatable ::uL(:),uR(:)
real*8,allocatable ::fL(:),fR(:)
real*8,allocatable ::cc(:),f(:)

allocate(uL(1:nx+1))
allocate(uR(1:nx+1))
allocate(fL(1:nx+1))
allocate(fR(1:nx+1))
allocate(cc(1:nx+1))
allocate(f(1:nx+1))

!compute upwind reconstruction for conserved variable u
call weno5L(nx,u,uL)

!compute downwind reconstruction for conserved variable u
call weno5R(nx,u,uR)

!compute fluxes
call fluxes(nx,uR,fR)
call fluxes(nx,uL,fL)

!compute wavespeed (Jacobian = df/du)
do i=2,nx
cc(i) = max(dabs(u(i)),dabs(u(i-1)))
end do
!periodic b.c.
cc(1) = max(dabs(u(1)),dabs(u(nx)))
cc(nx+1) = max(dabs(u(1)),dabs(u(nx)))

!Riemann solver
do i=1,nx+1
f(i) = 0.5d0*(fR(i)+fL(i)) - 0.5d0*cc(i)*(uR(i)-uL(i))
end do

!compute RHS
do i=1,nx	
	r(i) = -(f(i+1)-f(i))/dx
end do


return
end

!---------------------------------------------------------------------------!
!interface flux reconstruction formula
!---------------------------------------------------------------------------!

subroutine fluxes(nx,u,f)
implicit none
integer:: nx,i
real*8, dimension(1:nx+1) ::u,f  

do i=1,nx+1
f(i) = 0.5d0*u(i)*u(i)
end do

return
end


!-----------------------------------------------------------------------------!
!WENO5 reconstruction for upwind direction (positive & left to right)
!u(i): solution value at nodes i; i=1,2,...,N
!f(j): recontructed value at nodes j=i-1/2; j=1,...,N+1  
!periodic boundary condition
!-----------------------------------------------------------------------------!
subroutine weno5L(n,u,f)
implicit none
integer::n
real*8 ::u(1:n),f(1:n+1)
integer::i
real*8 ::a,b,c,d,e,w

i=0
  	a = u(i-2+n)
  	b = u(i-1+n)
  	c = u(i+n)
  	d = u(i+1)
  	e = u(i+2)
  	call weno5(a,b,c,d,e,w)
  	f(i+1) = w

i=1
  	a = u(i-2+n)
  	b = u(i-1+n)
  	c = u(i)
  	d = u(i+1)
  	e = u(i+2)
  	call weno5(a,b,c,d,e,w)
  	f(i+1) = w

i=2
  	a = u(i-2+n)
  	b = u(i-1)
  	c = u(i)
  	d = u(i+1)
  	e = u(i+2)
  	call weno5(a,b,c,d,e,w)
  	f(i+1) = w
       
do i=3,n-2
  	a = u(i-2)
  	b = u(i-1)
  	c = u(i)
  	d = u(i+1)
  	e = u(i+2)
  	call weno5(a,b,c,d,e,w)
  	f(i+1) = w
end do

i=n-1
  	a = u(i-2)
  	b = u(i-1)
  	c = u(i)
  	d = u(i+1)
  	e = u(i+2-n)
  	call weno5(a,b,c,d,e,w)
  	f(i+1) = w
    
i=n
  	a = u(i-2)
  	b = u(i-1)
  	c = u(i)
  	d = u(i+1-n)
  	e = u(i+2-n)
  	call weno5(a,b,c,d,e,w)
  	f(i+1) = w

return
end

!-----------------------------------------------------------------------------!
!WENO5 reconstruction for downwind direction (negative & right to left)
!u(i): solution value at nodes i; i=1,2,...,N
!f(j): recontructed value at nodes j=i-1/2; j=1,2,...,N+1 
!periodic boundary condition 
!-----------------------------------------------------------------------------!
subroutine weno5R(n,u,f)
implicit none
integer::n
real*8 ::u(1:n),f(1:n+1)
integer::i
real*8 ::a,b,c,d,e,w


i=1
  	a = u(i+2)
  	b = u(i+1)
  	c = u(i)
  	d = u(i-1+n)
  	e = u(i-2+n)
  	call weno5(a,b,c,d,e,w)
  	f(i) = w

i=2
  	a = u(i+2)
  	b = u(i+1)
  	c = u(i)
  	d = u(i-1)
  	e = u(i-2+n)
  	call weno5(a,b,c,d,e,w)
  	f(i) = w
        
do i=3,n-2
  	a = u(i+2)
  	b = u(i+1)
  	c = u(i)
  	d = u(i-1)
  	e = u(i-2)
  	call weno5(a,b,c,d,e,w)
  	f(i) = w
end do

i=n-1
  	a = u(i+2-n)
  	b = u(i+1)
  	c = u(i)
  	d = u(i-1)
  	e = u(i-2)
  	call weno5(a,b,c,d,e,w)
  	f(i) = w

i=n
  	a = u(i+2-n)
  	b = u(i+1-n)
  	c = u(i)
  	d = u(i-1)
  	e = u(i-2)
  	call weno5(a,b,c,d,e,w)
  	f(i) = w

i=n+1
  	a = u(i+2-n)
  	b = u(i+1-n)
  	c = u(i-n)
  	d = u(i-1)
  	e = u(i-2)
  	call weno5(a,b,c,d,e,w)
  	f(i) = w
    
return
end

!----------------------------------------------------------------------------------!
!WENO5 
!----------------------------------------------------------------------------------!
subroutine weno5(a,b,c,d,e,f)
implicit none
real*8 ::a,b,c,d,e,f
real*8 ::q1,q2,q3
real*8 ::s1,s2,s3
real*8 ::a1,a2,a3
real*8 ::eps

q1 = a/3.0d0 - 7.0d0/6.0d0*b + 11.0d0/6.0d0*c
q2 =-b/6.0d0 + 5.0d0/6.0d0*c + d/3.0d0
q3 = c/3.0d0 + 5.0d0/6.0d0*d - e/6.0d0

s1 = 13.0d0/12.0d0*(a-2.0d0*b+c)**2 + 0.25d0*(a-4.0d0*b+3.0d0*c)**2
s2 = 13.0d0/12.0d0*(b-2.0d0*c+d)**2 + 0.25d0*(d-b)**2
s3 = 13.0d0/12.0d0*(c-2.0d0*d+e)**2 + 0.25d0*(3.0d0*c-4.0d0*d+e)**2

!Jiang-Shu estimator
eps = 1.0d-6
a1 = 1.0d-1/(eps+s1)**2
a2 = 6.0d-1/(eps+s2)**2
a3 = 3.0d-1/(eps+s3)**2

!Shen-Zha estimator
!eps = 1.0d-20
!a1 = 1.0d-1*(1.0d0 + (dabs(s1-s3)/(eps+s1))**2)
!a2 = 6.0d-1*(1.0d0 + (dabs(s1-s3)/(eps+s2))**2)
!a3 = 3.0d-1*(1.0d0 + (dabs(s1-s3)/(eps+s3))**2)

f = (a1*q1 + a2*q2 + a3*q3)/(a1 + a2 + a3)

return
end subroutine



