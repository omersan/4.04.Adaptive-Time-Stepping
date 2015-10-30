!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Adaptive Runge Kutta Method for system of equations (Bryne and Hindmarsh)
!     Dormand-Prince method (ode45) for solving du/dt=rhs
!     
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) -- Ch.1
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Oct. 1, 2015
!-----------------------------------------------------------------------------!

program adaptiveRK
implicit none
integer,parameter::ne=2
real*8 ::a,b,hmin,hmax,sf,tol
real*8 ::u(ne)

call RKTABLE

!minimum and maximum available time steps
hmin= 1.0d-3
hmax= 1.0d2

!tolerance and safety factor
tol = 1.0d-6
sf  = 0.5d0

!time interval
a=0.0d0
b=7.0d5

!initial guess
u(1) =-1.0d0 
u(2) = 0.0d0 



open(20,file='num.plt')
write(20,*) 'variables ="t","u1","u2"'
write(20,*) a,u(1),u(2)

open(30,file='dt.plt')
write(30,*) 'variables ="t","dt"'

!Dormand-Prince method (ode45) for solving du/dt=f(u)
call DP45(a,b,hmin,hmax,sf,tol,ne,u)  !ode45

end


!-----------------------------------------------------------!
!Adaptive time stepping
!Dormand-Prince method (ode45) Method for solving du/dt=rhs
!-----------------------------------------------------------!
SUBROUTINE DP45(a,b,hmin,hmax,sf,tol,ne,u)
implicit none
integer::ne
real*8::a,b,hmin,hmax,sf,tol
real*8::u(ne)
integer::i,j,k
real*8::adp45(7,7),cdp45(7),bdp4(7),bdp5(7)
real*8::r1(ne),r2(ne),r3(ne),r4(ne),r5(ne),r6(ne),r7(ne),uu(ne),u4(ne),u5(ne)
real*8::dt,t,tt,est,ratio,err(ne),tiny

common /dp45table/ adp45,cdp45,bdp4,bdp5

tiny=1.0d-8

!initial step size (start optimistically)
dt=hmax

!t holds the current time 
t=a


j=0
i=0 !step counter
do while (t.le.b-tiny)

!compute r1,r2,r3,r4,r5,r6,r7
call rhs(ne,t,u,r1)
do k=1,ne
uu(k)=u(k)+dt*adp45(2,1)*r1(k)
end do
tt=t+cdp45(2)*dt
call rhs(ne,tt,uu,r2)
do k=1,ne
uu(k)=u(k)+dt*(adp45(3,1)*r1(k)+adp45(3,2)*r2(k))
end do
tt=t+cdp45(3)*dt
call rhs(ne,tt,uu,r3)
do k=1,ne
uu(k)=u(k)+dt*(adp45(4,1)*r1(k)+adp45(4,2)*r2(k)+adp45(4,3)*r3(k))
end do
tt=t+cdp45(4)*dt
call rhs(ne,tt,uu,r4)
do k=1,ne
uu(k)=u(k)+dt*(adp45(5,1)*r1(k)+adp45(5,2)*r2(k)+adp45(5,3)*r3(k)+adp45(5,4)*r4(k))
end do
tt=t+cdp45(5)*dt
call rhs(ne,tt,uu,r5)
do k=1,ne
uu(k)=u(k)+dt*(adp45(6,1)*r1(k)+adp45(6,2)*r2(k)+adp45(6,3)*r3(k)+adp45(6,4)*r4(k)+adp45(6,5)*r5(k))
end do
tt=t+cdp45(6)*dt
call rhs(ne,tt,uu,r6)
do k=1,ne
uu(k)=u(k)+dt*(adp45(7,1)*r1(k)+adp45(7,2)*r2(k)+adp45(7,3)*r3(k)+adp45(7,4)*r4(k)+adp45(7,5)*r5(k)+adp45(7,6)*r6(k))
end do
tt=t+cdp45(7)*dt
call rhs(ne,tt,uu,r7)

!compute forth-order solution
do k=1,ne
u4(k) = u(k) + dt*(bdp4(1)*r1(k)+bdp4(2)*r2(k)+bdp4(3)*r3(k)+bdp4(4)*r4(k)+bdp4(5)*r5(k)+bdp4(6)*r6(k)+bdp4(7)*r7(k))
end do
!compute fifth-order solution 
do k=1,ne
u5(k) = u(k) + dt*(bdp5(1)*r1(k)+bdp5(2)*r2(k)+bdp5(3)*r3(k)+bdp5(4)*r4(k)+bdp5(5)*r5(k)+bdp5(6)*r6(k)+bdp5(7)*r7(k))
end do
!estimated error per unit time, should be at most tol
do k=1,ne
err(k)=dabs(u5(k)-u4(k))/dt
end do

est =maxval(err)

!check
if (est.le.tol.or.dt.le.hmin) then !accept
t=t+dt
	do k=1,ne
	u(k)=u5(k) !accept the 5th order solution
	end do
i=i+1
print*,t
write(20,*) t,u(1),u(2)
write(30,*) t,dt
end if

ratio = (sf*tol/(est+tiny))**(1.0d0/4.0d0)
!adjust step size for next time step
dt=dt*ratio

if(dt.lt.hmin) dt=hmin
if(dt.gt.hmax) dt=hmax
if((t+dt).gt.b) dt=b-t

j=j+1
end do

print*,"number of time steps taken",i
print*,"total number of time steps taken",j
print*,"number of null steps",j-i

return
end




!---------------------------------------------------!
!Butcher Table
!---------------------------------------------------!
SUBROUTINE RKTABLE
implicit none
integer::i,j
real*8::adp45(7,7),cdp45(7),bdp4(7),bdp5(7)


common /dp45table/ adp45,cdp45,bdp4,bdp5

!Butcher table for Dormand-Prince method (ode45)
do i=1,7
do j=1,7
adp45(i,j)=0.0d0
end do
cdp45(i)=0.0d0
bdp4(i) =0.0d0
bdp5(i) =0.0d0
end do

adp45(2,1)=1.0d0/5.0d0
adp45(3,1)=3.0d0/40.0d0
adp45(3,2)=9.0d0/40.0d0
adp45(4,1)=44.0d0/45.0d0
adp45(4,2)=-56.0d0/15.0d0
adp45(4,3)=32.0d0/9.0d0
adp45(5,1)=19372.0d0/6561.0d0
adp45(5,2)=-25360.0d0/2187.0d0
adp45(5,3)=64448.0d0/6561.0d0
adp45(5,4)=-212.0d0/729.0d0
adp45(6,1)=9017.0d0/3168.0d0
adp45(6,2)=-355.0d0/33.0d0
adp45(6,3)=46732.0d0/5247.0d0
adp45(6,4)=49.0d0/176.0d0
adp45(6,5)=-5103.0d0/18656.0d0
adp45(7,1)=35.0d0/384.0d0
adp45(7,2)=0.0d0
adp45(7,3)=500.0d0/1113.0d0
adp45(7,4)=125.0d0/192.0d0
adp45(7,5)=-2187.0d0/6784.0d0
adp45(7,6)=11.0d0/84.0d0

cdp45(2)=1.0d0/5.0d0
cdp45(3)=3.0d0/10.0d0
cdp45(4)=4.0d0/5.0d0
cdp45(5)=8.0d0/9.0d0
cdp45(6)=1.0d0
cdp45(7)=1.0d0

bdp4(1)=5179.0d0/57600.0d0
bdp4(2)=0.0d0
bdp4(3)=7571.0d0/16695.0d0
bdp4(4)=393.0d0/640.0d0
bdp4(5)=-92097.0d0/339200.0d0
bdp4(6)=187.0d0/2100.0d0
bdp4(7)=1.0d0/40.0d0

bdp5(1)=35.0d0/384.0d0
bdp5(2)=0.0d0
bdp5(3)=500.0d0/1113.0d0
bdp5(4)=125.0d0/192.0d0
bdp5(5)=-2187.0d0/6784.0d0
bdp5(6)=11.0d0/84.0d0
bdp5(7)=0.0d0


Return 
End


!---------------------------------------------------!
!rhs function (Bryne and Hindmarsh=Laser oscillator)
!---------------------------------------------------!
subroutine rhs(ne,t,u,f)
integer::ne
real*8 ::t,u(ne),f(ne)
real*8 ::alpha,beta,gamma,rho,sigma,tau

alpha = 1.5d-18
beta  = 2.5d-6
gamma = 2.1d-6

rho = 0.6d0
sigma = 0.18d0
tau = 0.016d0

f(1) =-u(1)*(alpha*u(2)+beta) + gamma
f(2) = u(2)*(rho*u(1)-sigma) + tau*(1.0d0 + u(1))


return
end

