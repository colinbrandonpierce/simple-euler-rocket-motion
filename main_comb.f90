program main
implicit none
character(len=*), parameter :: FILE_NAME = 'data_comb.txt'

integer:: i, j!, fu
real :: g, Isp, m_r, m_f, f, b, h, tmin, tmid, tmax, R, g_s,  d, dd, csa
integer:: Na, Nb
real, dimension(:), allocatable :: p, v, m, t
real, dimension (:), allocatable :: pp, vv, mm, tt
i = 1 ! index
g_s = -9.81 ! parameters
Isp = 200
m_r = 5000
m_f = 1000
f = m_r-m_f
b = 60

R = 6.378*10**6
csa = 3*3.1415/4

tmin=0
tmid = 60
tmax=400

h = 0.1

Na = int((tmid-tmin)/h)+1+1
Nb = int((tmax-tmid)/h)+1+1

allocate(p(Na))
allocate(v(Na))
allocate(m(Na))
allocate(t(Na))
allocate(pp(Nb))
allocate(vv(Nb))
allocate(mm(Nb))
allocate(tt(Nb))

print*,Na, Nb
t(i)=0
p(i)=0
v(i)=0
m(i)=m_r

do while (i <= Na) ! simple euler
dd=1.225*exp(-0.1385*p(i))
d = (0.6*dd*csa*(v(i)**2))/2
g = g_s*(R**2/(R+p(i))**2)
t(i+1) = t(i) + h
p(i+1)=p(i)+ h*v(i)
v(i+1)=v(i)+ h*(g + -d + (Isp*g*(-f/b))/(m(i)) )
m(i+1)=m(i)+ h*(-f/b)
i = i+1
end do

j=1
tt(j)=t(i-1)
pp(j)=p(i-1)
vv(j)=v(i-1)
mm(j)=m(i-1)

do while (j < Nb)
dd=1.225*exp(-0.1385*pp(j))
d = (0.6*dd*csa*(vv(j)**2))/2
g = g_s*(R**2/(R+pp(i))**2)
tt(j+1)=tt(j)+ h
pp(j+1)=pp(j)+ h*vv(j)
vv(j+1)=vv(j)+ h*(g-d)
mm(j+1)=mm(j)
j = j+1
end do

print*, MAXVAL(p), MAXVAL(pp), MAXVAL(v), MAXVAL(vv)

!~ open (action='write', file=FILE_NAME, newunit=fu, status='replace')

!~ do i = 1, Na+Nb
!~ if (i<=Na) then
!~ write (fu, *) t(i), p(i), v(i), m(i)
!~ else
!~ write (fu, *) tt(i-Na), pp(i-Na), vv(i-Na), mm(i-Na)
!~ end if
!~ end do

!~ close (fu)

!~ call execute_command_line('gnuplot -p plot_g.plt')

end program main
