subroutine splinx_adj ( n, x, y, a_y, nn, xx, c, a_c, a_yy )
implicit none
integer :: n, nn
integer :: i, j, m, i1
real :: a_c(3,n), a_y(n), a_yy(n), c(3,n), x(n), xx(n), y(n)
real :: a_der, a_dy, a_dyl, a_dyr, a_dyx, dyx
real :: ch(3,n), ch_1(3,n), ch_2(3,n), ch_3(3,n), ch_4(3,n), ch_5(3,n), ch_6(3,n)
real :: d, del, dx, dxl, dxr, dy, dyl, dyr

!----------------------------------------------
! RESET LOCAL ADJOINT VARIABLES
!----------------------------------------------
a_der = 0.
a_dy = 0.
a_dyl = 0.
a_dyr = 0.
a_dyx = 0.

  ch(:,:) = c(:,:)

c(1,1) = 0.
c(2,1) = 0.
c(2,n) = 0.
dxr = x(2)-x(1)
dyr = (y(2)-y(1))/dxr

do i = 2, n-1
  dxl = dxr
  dxr = x(i+1)-x(i)
  dyl = dyr
  dyr = (y(i+1)-y(i))/dxr
  dx = dxr+dxl
  dy = (dyr-dyl)/dx
  c(1,i) = -(dxr/(2.*dx+dxl*c(1,i-1)))

  ch_2(1,i-1) = c(1,i-1)
  ch_6(2,i-1) = c(2,i-1)

  c(2,i) = (6.*dx*dy-dxl*c(2,i-1))/(2.*dx+dxl*c(1,i-1))
end do

  do i=n-1,2,-1

   ch_4(2,i+1)=c(2,i+1)
   ch_5(1,i)=c(1,i)

   c(2,i)=c(1,i)*c(2,i+1)+c(2,i)
  end do

do j = 1, nn
  a_der = 0.
  if (xx(j) .le. x(1)) then
    a_c(1,1) = a_c(1,1)+a_yy(j)*(xx(j)-x(1))
    a_y(1) = a_y(1)+a_yy(j)
    a_yy(j) = 0.
  else if (xx(j) .ge. x(n)) then
    del = x(n)-x(n-1)
    a_der = a_der+a_yy(j)*(xx(j)-x(n))
    a_y(n) = a_y(n)+a_yy(j)
    a_yy(j) = 0.
    a_c(3,n-1) = a_c(3,n-1)+3*a_der*del**2
    a_c(2,n-1) = a_c(2,n-1)+2*a_der*del
    a_c(1,n-1) = a_c(1,n-1)+a_der
    a_der = 0.
  else
    do i = 1, n
      if (x(i) .le. xx(j)) then
        m = i
      endif
    end do
    d = xx(j)-x(m)
    a_c(3,m) = a_c(3,m)+a_yy(j)*d*d*d
    a_c(2,m) = a_c(2,m)+a_yy(j)*d*d
    a_c(1,m) = a_c(1,m)+a_yy(j)*d
    a_y(m) = a_y(m)+a_yy(j)
    a_yy(j) = 0.
  endif
end do

do i = n-1,1,-1
  dx = x(i+1)-x(i)
  a_c(2,i) = a_c(2,i)-a_c(3,i)*(dx**2/dx**3)
  a_c(1,i) = a_c(1,i)-a_c(3,i)*(dx/dx**3)
  a_dy = a_dy+a_c(3,i)/dx**3
  a_c(3,i) = 0.
  a_c(2,i) = 0.5*a_c(2,i)
  a_c(2,i+1) = a_c(2,i+1)-0.16666667*a_c(1,i)*dx
  a_c(2,i) = a_c(2,i)-0.33333333*a_c(1,i)*dx
  a_dyx = a_dyx+a_c(1,i)
  a_c(1,i) = 0.
  a_dy = a_dy+a_dyx/dx
  a_dyx = 0.
  a_y(i+1) = a_y(i+1)+a_dy
  a_y(i) = a_y(i)-a_dy
  a_dy = 0.
end do

do i = 2, n-1
!do i = n-1,2,-1
   c(2,i+1)=ch_4(2,i+1)
   c(1,i)=ch_5(1,i)

  a_c(1,i) = a_c(1,i)+a_c(2,i)*c(2,i+1)
  a_c(2,i+1) = a_c(2,i+1)+a_c(2,i)*c(1,i)
end do

do i = n-1, 2, -1

  c(:,:) = ch(:,:)

  c(1,1) = 0.
  c(2,1) = 0.
  c(2,n) = 0.
  dxr = x(2)-x(1)
  dyr = (y(2)-y(1))/dxr
  do i1 = 2, i-1
    dxl = dxr
    dxr = x(i1+1)-x(i1)
    dyl = dyr
    dyr = (y(i1+1)-y(i1))/dxr
    dx = dxr+dxl
    dy = (dyr-dyl)/dx
    c(1,i1) = -(dxr/(2.*dx+dxl*c(1,i1-1)))
    c(2,i1) = (6.*dx*dy-dxl*c(2,i1-1))/(2.*dx+dxl*c(1,i1-1))
  end do
  ch_1(:,:) = c(:,:)

  dxl = dxr
  dxr = x(i+1)-x(i)
  dyl = dyr
  dyr = (y(i+1)-y(i))/dxr
  dx = dxr+dxl
  dy = (dyr-dyl)/dx
  c(1,i) = -(dxr/(2.*dx+dxl*c(1,i-1)))


  a_c(2,i-1) = a_c(2,i-1)-a_c(2,i)*(dxl/(2.*dx+dxl*c(1,i-1)))
  a_c(1,i-1) = a_c(1,i-1)-a_c(2,i)*((6.*dx*dy-dxl*c(2,i-1))*dxl/((2.*dx+dxl*c(1,i-1))*(2.*dx+dxl*c(1,i-1))))
  a_dy = a_dy+a_c(2,i)*(6*dx/(2.*dx+dxl*c(1,i-1)))
  a_c(2,i) = 0.
  c(:,:) = ch_1(:,:)
  a_c(1,i-1) = a_c(1,i-1)+a_c(1,i)*(dxr*dxl/((2.*dx+dxl*c(1,i-1))*(2.*dx+dxl*c(1,i-1))))
  a_c(1,i) = 0.
  a_dyl = a_dyl-a_dy/dx
  a_dyr = a_dyr+a_dy/dx
  a_dy = 0.
  a_y(i+1) = a_y(i+1)+a_dyr/dxr
  a_y(i) = a_y(i)-a_dyr/dxr
  a_dyr = 0.
  a_dyr = a_dyr+a_dyl
  a_dyl = 0.
end do

dxr = x(2)-x(1)
a_y(2) = a_y(2)+a_dyr/dxr
a_y(1) = a_y(1)-a_dyr/dxr
a_dyr = 0.
a_c(2,n) = 0.
a_c(2,1) = 0.
a_c(1,1) = 0.


end subroutine splinx_adj

!-----------------------------------------------------------------------------

subroutine mod_nonlocal_adj( num, jo_grad_y, xb, new_ob, mean_h, &
                             lowest_level, ref_mean_h_adj, ilocal )
      implicit none
      type (y_type), intent (in)     :: jo_grad_y
      type (xb_type), intent(in)        :: xb       ! first guess state.
      type (ob_in_mean_h), intent(in)   :: new_ob
      integer,intent(in) :: lowest_level
      integer, intent(in) :: num
      integer,dimension(kts:kte),intent(in) :: ilocal
      real,dimension(kts:kte),intent(in) :: mean_h
      real,dimension(ids:ide,jds:jde,kds:kde),intent(out):: ref_mean_h_adj

      integer :: i, j, k, l, m, n, i1, i2, i3, nbot, ntop
      integer :: is, ie, js, je, ks, ke
      integer :: i1l,i2l,i1r,i2r
      integer :: j1
      real :: dtr, rtd, step, dst, rox, glat, glon, h, rxabs, tmp_ref
      real :: xlatd, xlond, a_rfr
      real,dimension(ids:ide,jds:jde) :: x2, y2
      real,dimension(kts:kte) :: refm
      real,dimension(3):: rp, rn, rr, rx
      real,dimension(2) :: w1
      real,dimension(2,2) :: w2
      integer :: ip1,ip2,ip3,je2
      real :: refp, ray_length

      is = ids
      ie = ide
      js = jds
      je = jde
      ks = kds
      ke = kde

      dtr = pi/180.
      rtd = 180./pi
      step = 5.   !  km for ray integration

      a_rfr = 0.
      ref_mean_h_adj(:,:,:) = 0.0
     
       do j=js, je
       do i=is, ie
       if (global_lon(i,j) >= 0.) then
         x2(i,j)= global_lat(i,j)
         y2(i,j)= global_lon(i,j)
       else
         x2(i,j)= global_lat(i,j)
         y2(i,j)= (global_lon(i,j)+360.)
       endif
       enddo
       enddo
      
      nbot=lowest_level+1
      ntop=ke-1
      do i = nbot, ntop
!       rp(1)=(earth_radius+mean_h(i))*cos(dtr*new_ob%lat(i)) &
       rp(1)=(new_ob%rfict+mean_h(i))*cos(dtr*new_ob%lat(i)) &
                                     *cos(dtr*new_ob%lon(i))
!       rp(2)=(earth_radius+mean_h(i))*cos(dtr*new_ob%lat(i)) &
       rp(2)=(new_ob%rfict+mean_h(i))*cos(dtr*new_ob%lat(i)) &
                                     *sin(dtr*new_ob%lon(i))
!       rp(3)=(earth_radius+mean_h(i))*sin(dtr*new_ob%lat(i))
       rp(3)=(new_ob%rfict+mean_h(i))*sin(dtr*new_ob%lat(i))
       rn(1)=-sin(dtr*new_ob%lat(i))*cos(dtr*new_ob%lon(i))
       rn(2)=-sin(dtr*new_ob%lat(i))*sin(dtr*new_ob%lon(i))
       rn(3)=+cos(dtr*new_ob%lat(i))

       CALL spin(rn,rp,-dtr*new_ob%azim(i),rr)
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

      i1l=0
      i2l=0
      i1r=0
      i2r=0
      i1 = 0
      i2 = 0

      if (ilocal(i)==1)then   
      do j = 1, je2
       
       do k=-1,+1,2
        dst=step*j*k
 
        do l=1,3
        rx(l)=rp(l)+dst*rr(l)
        end do

       rox=sqrt(rx(1)**2+rx(2)**2)
       glat=rtd*atan(rx(3)/rox)
       glon=rtd*atan2(rx(2),rx(1))
       CALL absv(rx,rxabs)
!       h = rxabs - earth_radius
       h = rxabs - new_ob%rfict

       IF (h <= mean_h(ke-1))then

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
       if ( h < (global_terr(i1,i2)/1000.))then
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
       
       a_rfr = a_rfr + step * jo_grad_y%gpseph(num)%eph(i)
       
       do l=1,2
       do m=1,2
       do n=1,2
       ref_mean_h_adj(i1+l-1,i2+m-1,i3+n-1) =    &
       ref_mean_h_adj(i1+l-1,i2+m-1,i3+n-1) + w2(l,m)*w1(n)*a_rfr
       enddo
       enddo
       enddo
     
     a_rfr = 0.0

     ENDIF
      enddo
     ENDDO

     jo_grad_y%gpseph(num)%eph(i) = 0.
  
     else if (ilocal(i)==0) then
       a_rfr = a_rfr + step * jo_grad_y%gpseph(num)%eph(i)
       
       do l=1,2
       do m=1,2
       do n=1,2
       ref_mean_h_adj(ip1+l-1,ip2+m-1,ip3+n-1) =    &
       ref_mean_h_adj(ip1+l-1,ip2+m-1,ip3+n-1) + w2(l,m)*w1(n)*a_rfr
       enddo
       enddo
       enddo
     
     a_rfr = 0.0
     jo_grad_y%gpseph(num)%eph(i) = 0.

     endif 
     
     enddo

end subroutine mod_nonlocal_adj
      
