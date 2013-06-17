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
      real, dimension(its:ite,jts:jte,kts:kte),intent(in) :: ref_mean_h_ftl
      integer :: i, j, k, l, m, n, i1, i2, i3, nbot, ntop, num
      integer,intent(in) :: lowest_level
      integer :: j1
      integer :: is, ie, js, je, ks, ke
      integer :: i1l,i2l,i1r,i2r
      real :: dtr, rtd, step, dst, rox, glat, glon, h, rxabs, tmp_ref,tmp_ref9
      real :: xlatd, xlond, help_h
      real,dimension(its:ite,jts:jte) :: x2, y2
      real,dimension(kts:kte),intent(in) :: mean_h
      real,dimension(kts:kte) :: refm
      real,dimension(3):: rp, rn, rr, rx
      real,dimension(2) :: w1
      real,dimension(2,2) :: w2
      integer,dimension(kts:kte) :: ilocal
      integer :: ip1,ip2,ip3,je2
      real :: refp, ray_length


      dtr=pi/180.
      rtd=180./pi
      step=5.
      is = its
      ie = ite
      js = jts
      je = jte
      ks = kts
      ke = kte

      do i=is, ie
       do j=js, je
       if (xb%lon(i,j) >= 0.) then
         x2(i,j)= xb%lat(i,j)
         y2(i,j)= xb%lon(i,j)
       else
         x2(i,j)= xb%lat(i,j)
         y2(i,j)= (xb%lon(i,j)+360.)
       endif
       enddo
      enddo

      nbot=lowest_level+1
      ntop=ke-1
      do i=nbot,ntop
      rp(1)=(new_ob%rfict+mean_h(i))*cos(dtr*new_ob%lat(i)) &
                                    *cos(dtr*new_ob%lon(i))
      rp(2)=(new_ob%rfict+mean_h(i))*cos(dtr*new_ob%lat(i)) &
                                    *sin(dtr*new_ob%lon(i))
      rp(3)=(new_ob%rfict+mean_h(i))*sin(dtr*new_ob%lat(i))
      rn(1)=-(sin(dtr*new_ob%lat(i))*cos(dtr*new_ob%lon(i)))
      rn(2)=-(sin(dtr*new_ob%lat(i))*sin(dtr*new_ob%lon(i)))
      rn(3)=cos(dtr*new_ob%lat(i))
      help_h = -(dtr*new_ob%azim(i))
      CALL spin(rn,rp,help_h,rr)
      CALL rnorm(rr)

       rox=sqrt(rp(1)**2+rp(2)**2)
       glat=rtd*atan(rp(3)/rox)
       glon=rtd*atan2(rp(2),rp(1))
       CALL absv(rp,rxabs)
       h = rxabs - new_ob%rfict
       ray_length=sqrt((mean_h(ke)-h)*(2*new_ob%rfict+mean_h(ke)+h))
       je2=int(ray_length/step)+1
       ip1=0
       ip2=0
       ip3=0
       CALL find2(ie,je,rxabs,x2,y2,glat,glon,ip1,ip2,w2)

       CALL find1(ke,mean_h,h,ip3,w1)
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

      i1l=0
      i2l=0
      i1r=0
      i2r=0
      i1= 0
      i2= 0

       DO j=1,je2
       do k=-1, 1, 2
       dst=step*j*k
        do l=1,3
        rx(l)=rp(l)+dst*rr(l)
        end do
       rox=sqrt(rx(1)*rx(1)+rx(2)*rx(2))
       glat=rtd*atan(rx(3)/rox)
       glon=rtd*atan2(rx(2),rx(1))
       CALL absv(rx,rxabs)
!       h = rxabs - earth_radius
       h = rxabs - new_ob%rfict
       IF (h <= mean_h(ke-1)) then
!       CALL find1(ke,mean_h,h,i3,w1)
         if (k == -1)then
         i1=i1l
         i2=i2l
         endif
         if (k == 1) then
         i1=i1r
         i2=i2r
         endif
       CALL find2(ie,je,rxabs,x2,y2,glat,glon,i1,i2,w2)

! avoid boundry effect
        if (i1==1 .or. i1>=ie-1 .or. i2==1 .or. i2>=je-1) then
        w2(:,:)=0.
        endif
! avoid mountain effect
       CALL find1(ke,mean_h,h,i3,w1)
       if (h < (xb%terr(i1,i2)/1000.)) then
       w1(:)=0.
       endif

         if (k == -1)then
         i1l=i1
         i2l=i2
         endif
         if (k == 1) then
         i1r=i1
         i2r=i2
         endif
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
