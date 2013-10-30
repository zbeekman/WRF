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
      real,dimension(kts:kte),intent(inout) :: mean_h
      real,dimension(ids:ide,jds:jde,kds:kde),intent(out):: ref_mean_h_adj

      integer :: i, j, k, l, m, n, i1, i2, i3, nbot, ntop
      integer :: is, ie, js, je, ks, ke
      integer :: i1l,i2l,i1r,i2r
      integer :: j1
      real :: dtr, rtd, step, dst, rox, glat, glon, h, rxabs, tmp_ref
      real :: xlatd, xlond, a_rfr
!zx   real,dimension(ids:ide,jds:jde) :: x2, y2
      real,dimension(kts:kte) :: refm
      real,dimension(3):: rp, rn, rr, rx
      real,dimension(2) :: w1
      real,dimension(2,2) :: w2
      integer :: ip1,ip2,ip3,je2
      real :: refp, ray_length

!zx   is = ids
!zx   ie = ide
!zx   js = jds
!zx   je = jde
      ks = kds
      ke = kde

!zx   dtr = pi/180.
!zx   rtd = 180./pi
      step = 5.   !  km for ray integration

      a_rfr = 0.
      ref_mean_h_adj(:,:,:) = 0.0
     
!zx    do j=js, je
!zx    do i=is, ie
!zx    if (global_lon(i,j) >= 0.) then
!zx      x2(i,j)= global_lat(i,j)
!zx      y2(i,j)= global_lon(i,j)
!zx    else
!zx      x2(i,j)= global_lat(i,j)
!zx      y2(i,j)= (global_lon(i,j)+360.)
!zx    endif
!zx    enddo
!zx    enddo
      
      nbot=lowest_level+1
      nbot = gps_rays(num)%nbot

      ntop=ke-1
      ntop = gps_rays(num)%ntop

      mean_h = gps_rays(num)%mean_h

      do i = nbot, ntop
!       rp(1)=(earth_radius+mean_h(i))*cos(dtr*new_ob%lat(i)) &
!zx    rp(1)=(new_ob%rfict+mean_h(i))*cos(dtr*new_ob%lat(i)) &
!zx                                  *cos(dtr*new_ob%lon(i))
!       rp(2)=(earth_radius+mean_h(i))*cos(dtr*new_ob%lat(i)) &
!zx    rp(2)=(new_ob%rfict+mean_h(i))*cos(dtr*new_ob%lat(i)) &
!zx                                  *sin(dtr*new_ob%lon(i))
!       rp(3)=(earth_radius+mean_h(i))*sin(dtr*new_ob%lat(i))
!zx    rp(3)=(new_ob%rfict+mean_h(i))*sin(dtr*new_ob%lat(i))
!zx    rn(1)=-sin(dtr*new_ob%lat(i))*cos(dtr*new_ob%lon(i))
!zx    rn(2)=-sin(dtr*new_ob%lat(i))*sin(dtr*new_ob%lon(i))
!zx    rn(3)=+cos(dtr*new_ob%lat(i))

!zx    CALL spin(rn,rp,-dtr*new_ob%azim(i),rr)
!zx    CALL rnorm(rr)   

       
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

!zx   i1l=0
!zx   i2l=0
!zx   i1r=0
!zx   i2r=0
!zx   i1 = 0
!zx   i2 = 0

      if (gps_rays(num)%ilocal(i)==1)then   
      do j = 1, je2
       
       do k=-1,+1,2
!zx     dst=step*j*k
 
!zx     do l=1,3
!zx     rx(l)=rp(l)+dst*rr(l)
!zx     end do

!zx    rox=sqrt(rx(1)**2+rx(2)**2)
!zx    glat=rtd*atan(rx(3)/rox)
!zx    glon=rtd*atan2(rx(2),rx(1))
!zx    CALL absv(rx,rxabs)
!       h = rxabs - earth_radius
!zx    h = rxabs - new_ob%rfict
       h = gps_rays(num)%ip123(i)%h(k*j)

       IF (h <= mean_h(ke-1))then

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

! avoid boundry effect
!zx     if (i1==1 .or. i1>=ie-1 .or. i2==1 .or. i2>=je-1) then
!zx     w2(:,:)=0.
!zx     endif
       i1 = gps_rays(num)%ip123(i)%i1(k*j)
       i2 = gps_rays(num)%ip123(i)%i2(k*j)
       w2 = gps_rays(num)%ip123(i)%w2(1:2,1:2,k*j)

! avoid mountain effect
!zx    CALL find1(ke,mean_h,h,i3,w1)
!zx    if ( h < (global_terr(i1,i2)/1000.))then
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
  
     else if (gps_rays(num)%ilocal(i)==0) then
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
      
