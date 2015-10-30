!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Adaptive Runge Kutta Method 
!     for stiff problems
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

program aRK
implicit none

real*8::a,b,hmin,hmax,sf,tol
real*8::u
real*8::ua
external ua


call RKTABLE


!minimum and maximum available time steps
hmin= 1.0d-8
hmax= 1.0d-1

!tolerance and safety factor
tol = 1.0d-8
sf = 0.5

!time interval
a=0.0d0
b=10.0d0

!initial guess
u = 1.0d0 

open(20,file='num.plt')
write(20,*) 'variables ="t","u","dt"'
write(20,*) a,u,hmax

!call RK4adaptive(a,b,hmin,hmax,sf,tol,u)
!call RKF45(a,b,hmin,hmax,sf,tol,u)
!call DP45(a,b,hmin,hmax,sf,tol,u)  !ode45
call BS23(a,b,hmin,hmax,sf,tol,u)  !ode23

write(*,*)"time t=",b
write(*,*)"analytical result=",ua(b)
write(*,*)"computed result  =",u
write(*,*)"absolute error=",dabs(ua(b)-u)


end

!---------------------------------------------------!
!Adaptive time stepping
!Runge-Kutta-Fehlberg-5 Method for solving du/dt=rhs
!---------------------------------------------------!
SUBROUTINE RKF45(a,b,hmin,hmax,sf,tol,u)
implicit none
real*8::a,b,hmin,hmax,sf,tol
real*8::u
integer::i,j
real*8::arkf45(6,6),crkf45(6),brkf4(6),brkf5(6)
real*8::r1,r2,r3,r4,r5,r6,uu
real*8::dt,t,tt,est

common /rkf45table/arkf45,crkf45,brkf4,brkf5

!initial step size (start optimistically)
dt=hmax

!t holds the current time 
t=a

!u holds the current u(t)
u=u

j=0
i=0 !step counter
do while (t.le.b-hmin/2.0d0)

!compute r1,r2,r3,r4
call rhs(t,u,r1)
uu=u+dt*arkf45(2,1)*r1
tt=t+crkf45(2)*dt
call rhs(tt,uu,r2)
uu=u+dt*(arkf45(3,1)*r1+arkf45(3,2)*r2)
tt=t+crkf45(3)*dt
call rhs(tt,uu,r3)
uu=u+dt*(arkf45(4,1)*r1+arkf45(4,2)*r2+arkf45(4,3)*r3)
tt=t+crkf45(4)*dt
call rhs(tt,uu,r4)
uu=u+dt*(arkf45(5,1)*r1+arkf45(5,2)*r2+arkf45(5,3)*r3+arkf45(5,4)*r4)
tt=t+crkf45(5)*dt
call rhs(tt,uu,r5)
uu=u+dt*(arkf45(6,1)*r1+arkf45(6,2)*r2+arkf45(6,3)*r3+arkf45(6,4)*r4+arkf45(6,5)*r5)
tt=t+crkf45(6)*dt
call rhs(tt,uu,r6)

!compute forth-order solution
uu = u + dt*(brkf4(1)*r1+brkf4(2)*r2+brkf4(3)*r3+brkf4(4)*r4+brkf4(5)*r5+brkf4(6)*r6)

!compute fifth-order solution 
r5 = u + dt*(brkf5(1)*r1+brkf5(2)*r2+brkf5(3)*r3+brkf5(4)*r4+brkf5(5)*r5+brkf5(6)*r6)

!estimated error per unit time, should be at most tol
est=dabs(r5-uu)/dt

!check
if (est.le.tol.or.dt.le.hmin) then !accept
t=t+dt
u=r5
i=i+1
write(20,*) t,u,dt
end if

!adjust step size for next time step
dt=dt*(sf*tol/est)**(1.0d0/4.0d0)

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

!-----------------------------------------------------------!
!Adaptive time stepping
!Dormand-Prince method (ode45) Method for solving du/dt=rhs
!-----------------------------------------------------------!
SUBROUTINE DP45(a,b,hmin,hmax,sf,tol,u)
implicit none
real*8::a,b,hmin,hmax,sf,tol
real*8::u
integer::i,j
real*8::adp45(7,7),cdp45(7),bdp4(7),bdp5(7)
real*8::r1,r2,r3,r4,r5,r6,r7,uu
real*8::dt,t,tt,est

common /dp45table/ adp45,cdp45,bdp4,bdp5


!initial step size (start optimistically)
dt=hmax

!t holds the current time 
t=a

!u holds the current u(t)
u=u

j=0
i=0 !step counter
do while (t.le.b-hmin/2.0d0)

!compute r1,r2,r3,r4
call rhs(t,u,r1)
uu=u+dt*adp45(2,1)*r1
tt=t+cdp45(2)*dt
call rhs(tt,uu,r2)
uu=u+dt*(adp45(3,1)*r1+adp45(3,2)*r2)
tt=t+cdp45(3)*dt
call rhs(tt,uu,r3)
uu=u+dt*(adp45(4,1)*r1+adp45(4,2)*r2+adp45(4,3)*r3)
tt=t+cdp45(4)*dt
call rhs(tt,uu,r4)
uu=u+dt*(adp45(5,1)*r1+adp45(5,2)*r2+adp45(5,3)*r3+adp45(5,4)*r4)
tt=t+cdp45(5)*dt
call rhs(tt,uu,r5)
uu=u+dt*(adp45(6,1)*r1+adp45(6,2)*r2+adp45(6,3)*r3+adp45(6,4)*r4+adp45(6,5)*r5)
tt=t+cdp45(6)*dt
call rhs(tt,uu,r6)
uu=u+dt*(adp45(7,1)*r1+adp45(7,2)*r2+adp45(7,3)*r3+adp45(7,4)*r4+adp45(7,5)*r5+adp45(7,6)*r6)
tt=t+cdp45(7)*dt
call rhs(tt,uu,r7)

!compute forth-order solution
uu = u + dt*(bdp4(1)*r1+bdp4(2)*r2+bdp4(3)*r3+bdp4(4)*r4+bdp4(5)*r5+bdp4(6)*r6+bdp4(7)*r7)

!compute fifth-order solution 
r5 = u + dt*(bdp5(1)*r1+bdp5(2)*r2+bdp5(3)*r3+bdp5(4)*r4+bdp5(5)*r5+bdp5(6)*r6+bdp5(7)*r7)

!estimated error per unit time, should be at most tol
est=dabs(r5-uu)/dt

!check
if (est.le.tol.or.dt.le.hmin) then !accept
t=t+dt
u=r5
i=i+1
write(20,*) t,u,dt
end if

!adjust step size for next time step
dt=dt*(sf*tol/est)**(1.0d0/4.0d0)

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

!-----------------------------------------------------------!
!Adaptive time stepping
!Bogacki-Shampine Method(ode23) for solving du/dt=rhs
!-----------------------------------------------------------!
SUBROUTINE BS23(a,b,hmin,hmax,sf,tol,u)
implicit none
real*8::a,b,hmin,hmax,sf,tol
real*8::u
integer::i,j
real*8::abs23(4,4),cbs23(4),bbs2(4),bbs3(4)
real*8::r1,r2,r3,r4,uu
real*8::dt,t,tt,est

common /bs23table/abs23,cbs23,bbs2,bbs3

!initial step size (start optimistically)
dt=hmax

!t holds the current time 
t=a

!u holds the current u(t)
u=u

j=0
i=0 !step counter
do while (t.le.b-hmin/2.0d0)

!compute r1,r2,r3,r4
call rhs(t,u,r1)
uu=u+dt*abs23(2,1)*r1
tt=t+cbs23(2)*dt
call rhs(tt,uu,r2)
uu=u+dt*(abs23(3,1)*r1+abs23(3,2)*r2)
tt=t+cbs23(3)*dt
call rhs(tt,uu,r3)
uu=u+dt*(abs23(4,1)*r1+abs23(4,2)*r2+abs23(4,3)*r3)
tt=t+cbs23(4)*dt
call rhs(tt,uu,r4)

!compute second-order solution
uu = u + dt*(bbs2(1)*r1+bbs2(2)*r2+bbs2(3)*r3+bbs2(4)*r4)

!compute third-order solution 
r3 = u + dt*(bbs3(1)*r1+bbs3(2)*r2+bbs3(3)*r3+bbs3(4)*r4)

!estimated error per unit time, should be at most tol
est=dabs(r3-uu)/dt

!check
if (est.le.tol.or.dt.le.hmin) then !accept
t=t+dt
u=r3
i=i+1
write(20,*) t,u,dt
end if

!adjust step size for next time step
dt=dt*(sf*tol/est)**(1.0d0/2.0d0)

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


!-------------------------------------------------------!
!Adaptive time stepping (Richardson type error estimate)
!Regular Runge-Kutta-4 Method for solving du/dt=rhs
!time interval; [a,b]
!tolerance; tol
!maximum and minimum allowable time step; hmax, hmin
!-------------------------------------------------------!
SUBROUTINE RK4adaptive(a,b,hmin,hmax,sf,tol,u)
implicit none
real*8::a,b,hmin,hmax,sf,tol
real*8::u
integer::i,j
real*8::dt,t,tnew,tmid,thalf,unew,umid,uhalf,est

!initial step size (start optimistically)
dt=hmax

!t holds the current time 
t=a

!u holds the current u(t)
u=u

j=0
i=0 !step counter
do while (t.le.b-hmin/2.0d0)
	
	!provisional new values
	tnew=t
	unew=u	
	call RK4(tnew,dt,unew)

	!now do two steps of lenght dt/2 for Richardson error estimate
	tmid=t
	umid=u	
	call RK4(tmid,dt/2.0d0,umid)
	thalf=tmid
	uhalf=umid	
	call RK4(thalf,dt/2.0d0,uhalf)

	!estimated error per unit time, should be at most tol
	est=dabs(uhalf-unew)/dt

	!check
	if (est.le.tol.or.dt.le.hmin) then !accept
	t=tnew
	u=unew
	i=i+1
	write(20,*) t,u,dt
	end if

	!adjust step size for next time step
	dt=dt*(sf*tol/est)**(1.0d0/4.0d0)

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
!Regular Runge-Kutta-4 Method for solving du/dt=rhs
!---------------------------------------------------!
SUBROUTINE RK4(t,h,u)
implicit none
real*8::t,h,u
real*8::ark4(4,4),crk4(4),brk4(4)
real*8::r1,r2,r3,r4,uu
real*8::tt

common /rk4table/ark4,crk4,brk4

!compute r1,r2,r3,r4
call rhs(t,u,r1)
uu=u+h*ark4(2,1)*r1
tt=t+crk4(2)*h
call rhs(tt,uu,r2)
uu=u+h*(ark4(3,1)*r1+ark4(3,2)*r2)
tt=t+crk4(3)*h
call rhs(tt,uu,r3)
uu=u+h*(ark4(4,1)*r1+ark4(4,2)*r2+ark4(4,3)*r3)
tt=t+crk4(4)*h
call rhs(tt,uu,r4)

!compute fourth-order solution and update
u = u + h*(brk4(1)*r1+brk4(2)*r2+brk4(3)*r3+brk4(4)*r4)
t = t + h

return
end



!---------------------------------------------------!
!Butcher Tables for Explicit Runge-Kutta Methods
!---------------------------------------------------!
SUBROUTINE RKTABLE
implicit none
integer::i,j
!Runge-Kutta Coefficients
real*8::abs23(4,4),cbs23(4),bbs2(4),bbs3(4)
real*8::ark4(4,4),crk4(4),brk4(4)
real*8::arkf45(6,6),crkf45(6),brkf4(6),brkf5(6)
real*8::adp45(7,7),cdp45(7),bdp4(7),bdp5(7)


common /bs23table/abs23,cbs23,bbs2,bbs3
common /rk4table/ark4,crk4,brk4
common /rkf45table/arkf45,crkf45,brkf4,brkf5
common /dp45table/ adp45,cdp45,bdp4,bdp5

!Ititialize Butcher table
do i=1,4
do j=1,4
abs23(i,j)=0.0d0
ark4(i,j) =0.0d0
end do
cbs23(i)=0.0d0
bbs2(i) =0.0d0
bbs3(i) =0.0d0
crk4(i) =0.0d0
brk4(i) =0.0d0
end do


do i=1,6
do j=1,6
arkf45(i,j)=0.0d0
end do
crkf45(i)=0.0d0
brkf4(i) =0.0d0
brkf5(i) =0.0d0
end do


do i=1,7
do j=1,7
adp45(i,j)=0.0d0
end do
cdp45(i)=0.0d0
bdp4(i) =0.0d0
bdp5(i) =0.0d0
end do

!Butcher table for runge-kutta 4
ark4(2,1)=1.0d0/2.0d0
ark4(3,1)=0.0d0
ark4(3,2)=1.0d0/2.0d0
ark4(4,1)=0.0d0
ark4(4,2)=0.0d0
ark4(4,3)=1.0d0

crk4(2)=1.0d0/2.0d0
crk4(3)=1.0d0/2.0d0
crk4(4)=1.0d0


brk4(1)=1.0d0/6.0d0
brk4(2)=2.0d0/6.0d0
brk4(3)=2.0d0/6.0d0
brk4(4)=1.0d0/6.0d0

!Butcher table for Bogacki-Shampine method (ode23)
abs23(2,1)=1.0d0/2.0d0
abs23(3,1)=0.0d0
abs23(3,2)=3.0d0/4.0d0
abs23(4,1)=2.0d0/9.0d0
abs23(4,2)=3.0d0/9.0d0
abs23(4,3)=4.0d0/9.0d0

cbs23(2)=1.0d0/2.0d0
cbs23(3)=3.0d0/4.0d0
cbs23(4)=1.0d0

bbs2(1)=7.0d0/24.0d0
bbs2(2)=1.0d0/4.0d0
bbs2(3)=1.0d0/3.0d0
bbs2(4)=1.0d0/8.0d0


bbs3(1)=2.0d0/9.0d0
bbs3(2)=3.0d0/9.0d0
bbs3(3)=4.0d0/9.0d0
bbs3(4)=0.0d0


!Butcher table for Runge-Kutta-Fehlberg method 
arkf45(2,1)=1.0d0/4.0d0
arkf45(3,1)=3.0d0/32.0d0
arkf45(3,2)=9.0d0/32.0d0
arkf45(4,1)=1932.0d0/2197.0d0
arkf45(4,2)=-7200.0d0/2197.0d0
arkf45(4,3)=7296.0d0/2197.0d0
arkf45(5,1)=439.0d0/216.0d0
arkf45(5,2)=-8.0d0
arkf45(5,3)=3680.0d0/513.0d0
arkf45(5,4)=-845.0d0/4104.0d0
arkf45(6,1)=-8.0d0/27.0d0
arkf45(6,2)=2.0d0
arkf45(6,3)=-3544.0d0/2565.0d0
arkf45(6,4)=1859.0d0/4104.0d0
arkf45(6,5)=-11.0d0/40.0d0

crkf45(2)=1.0d0/4.0d0
crkf45(3)=3.0d0/8.0d0
crkf45(4)=12.0d0/13.0d0
crkf45(5)=1.0d0
crkf45(6)=1.0d0/2.0d0

brkf4(1)=25.0d0/216.0d0
brkf4(2)=0.0d0
brkf4(3)=1408.0d0/2565.0d0
brkf4(4)=2197.0d0/4104.0d0
brkf4(5)=-1.0d0/5.0d0
brkf4(6)=0.0d0

brkf5(1)=16.0d0/135.0d0
brkf5(2)=0.0d0
brkf5(3)=6656.0d0/12825.0d0
brkf5(4)=28561.0d0/56430.0d0
brkf5(5)=-9.0d0/50.0d0
brkf5(6)=2.0d0/55.0d0


!Butcher table for Dormand-Prince method (ode45)
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






!--------------!
!rhs function
!--------------!
subroutine rhs(t,u,f)
real*8::t,u,f

!f = u**2.0d0
!f = 1.0d0-10.0d0*u
!f = -2.0d0*(t-1.0d0)*u
f = -u + dsin(t)

return
end

!--------------!
!exact solution
!--------------!
real*8 function ua(t)
implicit none
real*8::t
ua = (dsin(t)-dcos(t))/2.0d0 + 3.0d0/2.0d0 * dexp(-t)
return
end
