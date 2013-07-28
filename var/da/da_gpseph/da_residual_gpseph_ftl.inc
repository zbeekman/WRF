!SUBROUTINE splinx_ftl ( n9, x9, y9, y, c9, c, nn9, xx9, yy9, yy )
subroutine splinx_ftl( n, x, y, g_y, c, g_c, nn, xx, yy, g_yy )

implicit none
integer n
integer nn
real(kind=8) c(3,n)
real(kind=8) g_c(3,n)
real(kind=8) g_y(n)
real(kind=8) g_yy(n)
real(kind=8) x(n)
real(kind=8) xx(n)
real(kind=8) y(n)
real(kind=8) yy(n)
real(kind=8) d
real(kind=8) del
real(kind=8) der
real(kind=8) dx
real(kind=8) dxl
real(kind=8) dxr
real(kind=8) dy
real(kind=8) dyl
real(kind=8) dyr
real(kind=8) dyx
real(kind=8) g_der
real(kind=8) g_dy
real(kind=8) g_dyl
real(kind=8) g_dyr
real(kind=8) g_dyx
integer i
integer j
integer m
g_c(1,1) = 0.
c(1,1) = 0.
g_c(2,1) = 0.
c(2,1) = 0.
g_c(2,n) = 0.
c(2,n) = 0.
dxr = x(2)-x(1)
g_dyr = g_y(2)/dxr-g_y(1)/dxr
dyr = (y(2)-y(1))/dxr
do i = 2, n-1
  dxl = dxr
  dxr = x(i+1)-x(i)
  g_dyl = g_dyr
  dyl = dyr
  g_dyr = g_y(i+1)/dxr-g_y(i)/dxr
  dyr = (y(i+1)-y(i))/dxr
  dx = dxr+dxl
  g_dy = (-(g_dyl/dx))+g_dyr/dx
  dy = (dyr-dyl)/dx
  g_c(1,i) = g_c(1,i-1)*(dxr*dxl/((2.*dx+dxl*c(1,i-1))*(2.*dx+dxl*c(1,i-1))))
  c(1,i) = -(dxr/(2.*dx+dxl*c(1,i-1)))
  g_c(2,i) = (-(g_c(2,i-1)*(dxl/(2.*dx+dxl*c(1,i-1)))+g_c(1,i-1)  &
            *((6.*dx*dy-dxl*c(2,i-1))*dxl/((2.*dx+dxl*c(1,i-1))   &
            *(2.*dx+dxl*c(1,i-1))))))+g_dy*(6*dx/(2.*dx+dxl*c(1,i-1)))
  c(2,i) = (6.*dx*dy-dxl*c(2,i-1))/(2.*dx+dxl*c(1,i-1))
end do
do i = n-1, 2, -1
  g_c(2,i) = g_c(2,i+1)*c(1,i)+g_c(2,i)+g_c(1,i)*c(2,i+1)
  c(2,i) = c(1,i)*c(2,i+1)+c(2,i)
end do
do i = 1, n-1
  dx = x(i+1)-x(i)
  g_dy = g_y(i+1)-g_y(i)
  dy = y(i+1)-y(i)
  g_dyx = g_dy/dx
  dyx = dy/dx
  g_c(1,i) = (-(0.16666667*g_c(2,i+1)*dx+0.33333333*g_c(2,i)*dx))+g_dyx
  c(1,i) = dyx-dx*(c(2,i)/3.+c(2,i+1)/6.)
  g_c(2,i) = 0.5*g_c(2,i)
  c(2,i) = c(2,i)/2.
  g_c(3,i) = (-(g_c(2,i)*(dx**2/dx**3)+g_c(1,i)*(dx/dx**3)))+g_dy/dx**3
  c(3,i) = (dy-c(1,i)*dx-c(2,i)*dx**2)/dx**3
end do
do j = 1, nn
  if (xx(j) .le. x(1)) then
    g_yy(j) = g_c(1,1)*(xx(j)-x(1))+g_y(1)
    yy(j) = y(1)+c(1,1)*(xx(j)-x(1))
  else if (xx(j) .ge. x(n)) then
    del = x(n)-x(n-1)
    g_der = 3*g_c(3,n-1)*del**2+2*g_c(2,n-1)*del+g_c(1,n-1)
    der = c(1,n-1)+2.d0*c(2,n-1)*del+3.d0*c(3,n-1)*del**2
    g_yy(j) = g_der*(xx(j)-x(n))+g_y(n)
    yy(j) = y(n)+der*(xx(j)-x(n))
  else
    do i = 1, n
      if (x(i) .le. xx(j)) then
        m = i
      endif
    end do
    d = xx(j)-x(m)
    g_yy(j) = g_c(3,m)*d*d*d+g_c(2,m)*d*d+g_c(1,m)*d+g_y(m)
    yy(j) = y(m)+d*(c(1,m)+d*(c(2,m)+d*c(3,m)))
  endif
end do

end subroutine splinx_ftl

!-------------------------------------------------------------------

subroutine mod_nonlocal_ftl (num, xb, new_ob, mean_h, lowest_level, &
                             ref_mean_h_ftl , y)
      IMPLICIT NONE
      TYPE (xb_type), INTENT(IN)     :: xb       ! first guess state.
      TYPE (ob_in_mean_h), INTENT(IN)   :: new_ob
      TYPE (y_type), INTENT(INout)   :: y
!      real, dimension(its:ite,jts:jte,kts:kte),intent(in) :: ref_mean_h
      real, dimension(ids:ide,jds:jde,kds:kde),intent(in) :: ref_mean_h_ftl
      integer :: i, j, k, l, m, n, i1, i2, i3, nbot, ntop, num
      integer,intent(in) :: lowest_level
      integer :: j1
      integer :: is, ie, js, je, ks, ke
      integer :: i1l,i2l,i1r,i2r
      real :: dtr, rtd, step, dst, rox, glat, glon, h, rxabs, tmp_ref,tmp_ref9
      real :: xlatd, xlond, help_h
!zx   real,dimension(ids:ide,jds:jde) :: x2, y2
      real,dimension(kts:kte),intent(inout) :: mean_h
      real,dimension(kts:kte) :: refm
      real,dimension(3):: rp, rn, rr, rx
      real,dimension(2) :: w1
      real,dimension(2,2) :: w2
      integer,dimension(kts:kte) :: ilocal
      integer :: ip1,ip2,ip3,je2
      real :: refp, ray_length


!zx   dtr=pi/180.
!zx   rtd=180./pi
      step=5.
!zx   is = ids
!zx   ie = ide
!zx   js = jds
!zx   je = jde
      ks = kds
      ke = kde

!zx   do j=js, je
!zx    do i=is, ie
!zx    if (global_lon(i,j) >= 0.) then
!zx      x2(i,j)= global_lat(i,j)
!zx      y2(i,j)= global_lon(i,j)
!zx    else
!zx      x2(i,j)= global_lat(i,j)
!zx      y2(i,j)= (global_lon(i,j)+360.)
!zx    endif
!zx    enddo
!zx   enddo

      nbot=lowest_level+1
      nbot = gps_rays(num)%nbot
      ntop=ke-1
      ntop = gps_rays(num)%ntop

      mean_h = gps_rays(num)%mean_h

      do i=nbot,ntop
!zx   rp(1)=(new_ob%rfict+mean_h(i))*cos(dtr*new_ob%lat(i)) &
!zx                                 *cos(dtr*new_ob%lon(i))
!zx   rp(2)=(new_ob%rfict+mean_h(i))*cos(dtr*new_ob%lat(i)) &
!zx                                 *sin(dtr*new_ob%lon(i))
!zx   rp(3)=(new_ob%rfict+mean_h(i))*sin(dtr*new_ob%lat(i))
!zx   rn(1)=-(sin(dtr*new_ob%lat(i))*cos(dtr*new_ob%lon(i)))
!zx   rn(2)=-(sin(dtr*new_ob%lat(i))*sin(dtr*new_ob%lon(i)))
!zx   rn(3)=cos(dtr*new_ob%lat(i))
!zx   help_h = -(dtr*new_ob%azim(i))
!zx   CALL spin(rn,rp,help_h,rr)
!zx   CALL rnorm(rr)

!zx    rox=sqrt(rp(1)**2+rp(2)**2)
!zx    glat=rtd*atan(rp(3)/rox)
!zx    glon=rtd*atan2(rp(2),rp(1))
!zx    CALL absv(rp,rxabs)
!zx    h = rxabs - new_ob%rfict
       h = gps_rays(num)%ip123(i)%h(0)
!zx    ray_length=sqrt((mean_h(ke)-h)*(2*new_ob%rfict+mean_h(ke)+h))
!zx    je2=int(ray_length/step)+1
       je2 = gps_rays(num)%je2(i)
!zx    ip1=0
!zx    ip2=0
!zx    ip3=0
!zx    CALL find2(ie,je,rxabs,x2,y2,glat,glon,ip1,ip2,w2)
       ip1 = gps_rays(num)%ip123(i)%i1(0)
       ip2 = gps_rays(num)%ip123(i)%i2(0)
       w2 = gps_rays(num)%ip123(i)%w2(1:2,1:2,0)

!zx    CALL find1(ke,mean_h,h,ip3,w1)
       ip3 = gps_rays(num)%ip123(i)%i3(0)
       w1 = gps_rays(num)%ip123(i)%w1(1:2,0)

       refp=0.
        do l=1,2
        do m=1,2
        do n=1,2
        refp=refp+ref_mean_h_ftl(ip1+l-1,ip2+m-1,ip3+n-1)*w2(l,m)*w1(n)
        end do
        end do
        end do
       ilocal(i)=1

       y%gpseph(num)%eph(i)=0.0

!zx   i1l=0
!zx   i2l=0
!zx   i1r=0
!zx   i2r=0
!zx   i1= 0
!zx   i2= 0

       DO j=1,je2
       do k=-1, 1, 2
!zx    dst=step*j*k
!zx     do l=1,3
!zx     rx(l)=rp(l)+dst*rr(l)
!zx     end do
!zx    rox=sqrt(rx(1)*rx(1)+rx(2)*rx(2))
!zx    glat=rtd*atan(rx(3)/rox)
!zx    glon=rtd*atan2(rx(2),rx(1))
!zx    CALL absv(rx,rxabs)
!       h = rxabs - earth_radius
!zx    h = rxabs - new_ob%rfict
       h = gps_rays(num)%ip123(i)%h(k*j)
       IF (h <= mean_h(ke-1)) then
!       CALL find1(ke,mean_h,h,i3,w1)
!zx      if (k == -1)then
!zx      i1=i1l
!zx      i2=i2l
!zx      endif
!zx      if (k == 1) then
!zx      i1=i1r
!zx      i2=i2r
!zx      endif
!zx    CALL find2(ie,je,rxabs,x2,y2,glat,glon,i1,i2,w2)
       i1 = gps_rays(num)%ip123(i)%i1(k*j)
       i2 = gps_rays(num)%ip123(i)%i2(k*j)
       w2 = gps_rays(num)%ip123(i)%w2(1:2,1:2,k*j)


! avoid boundry effect
!zx     if (i1==1 .or. i1>=ie-1 .or. i2==1 .or. i2>=je-1) then
!zx     w2(:,:)=0.
!zx     endif
! avoid mountain effect
!zx    CALL find1(ke,mean_h,h,i3,w1)
!zx    if (h < (global_terr(i1,i2)/1000.)) then
!zx    w1(:)=0.
!zx    endif
       i3 = gps_rays(num)%ip123(i)%i3(k*j)
       w1 = gps_rays(num)%ip123(i)%w1(1:2,k*j)

!zx      if (k == -1)then
!zx      i1l=i1
!zx      i2l=i2
!zx      endif
!zx      if (k == 1) then
!zx      i1r=i1
!zx      i2r=i2
!zx      endif
       tmp_ref=0.0
        do l=1,2
        do m=1,2
        do n=1,2
!        tmp_ref=tmp_ref+xa%ref(i1+l-1,i2+m-1,i3+n-1)*w2(l,m)*w1(n)
        tmp_ref=tmp_ref+ref_mean_h_ftl(i1+l-1,i2+m-1,i3+n-1)*w2(l,m)*w1(n)
        end do
        end do
        end do
       y%gpseph(num)%eph(i)=y%gpseph(num)%eph(i)+step*tmp_ref
       END IF
       end do
       END DO
      if ( .not. gpseph_local ) then 
        if (y%gpseph(num)%eph(i) == 0.) then
        y%gpseph(num)%eph(i)=step*refp
        ilocal(i)=0
        endif
      else if ( gpseph_local ) then
        ilocal(i)=0
        y%gpseph(num)%eph(i)=step*refp
      endif

      end do
END SUBROUTINE mod_nonlocal_ftl
