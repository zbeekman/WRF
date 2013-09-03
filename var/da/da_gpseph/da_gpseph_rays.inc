SUBROUTINE da_gpseph_rays ( n, xb, iv, mean_h, lowest_level, new_ob)

      implicit none

      type (xb_type), intent(in)     :: xb       ! first guess state.
      type (ob_in_mean_h), intent(inout)   :: new_ob
      type (iv_type), intent(in)      :: iv       ! O-B structure.
      integer, intent(in)             :: n
       
      real,dimension(kts:kte),intent(out) :: mean_h
      integer,intent(out) :: lowest_level

      integer :: i, j, k, l, m, i1, i2, i3, nbot, ntop
      integer :: is, ie, js, je, ks, ke
      integer :: i1l,i2l,i1r,i2r
      integer                         :: nlevel        ! the count of observational levels
      integer :: ip1,ip2,ip3,je2
      real :: dtr, rtd, step, dst, rox, glat, glon, h, rxabs
      real :: xlatd, xlond
      real,dimension(ids:ide,jds:jde) :: x2, y2
      real,dimension(kts:kte) :: refm
      real,dimension(3):: rp, rn, rr, rx
      real,dimension(2) :: w1
      real,dimension(2,2) :: w2  
      real :: refp, ray_length

      real                            :: nstart, nstop , temp_lat
      real, dimension(:), allocatable :: log_ref
      real, dimension(:), allocatable :: lat_inv
      real, dimension(:), allocatable :: lon_inv
      real, dimension(:), allocatable :: azim_inv
      real, dimension(:), allocatable :: tmp_h
      real, dimension(interpolate_level) :: tmp_z   ! obs. levels are interpolated
                                                    ! into the amount of interpolate_level
      real, dimension(interpolate_level) :: tmp_ref
      real, dimension(interpolate_level) :: tmp_lat
      real, dimension(interpolate_level) :: tmp_lon
      real, dimension(interpolate_level) :: tmp_azim

      real :: requator,rpolar
      requator=6378.1363
      rpolar=6356.7516

      temp_lat=(iv%info(gpseph)%lat(1,n)*pi)/180.
      new_ob%rfict=SQRT((requator**2*rpolar**2)/                    &
                        (requator**2*(sin(abs(temp_lat)))**2+  &
                         rpolar**2*(cos(abs(temp_lat)))**2))
      CALL model_mean_altitude(xb, new_ob, mean_h)          ! mean_h in km

      allocate( gps_rays(n)%mean_h(kts:kte) )
      gps_rays(n)%mean_h = mean_h

      allocate( log_ref (iv%info(gpseph)%levels(n)) )
      allocate( lat_inv (iv%info(gpseph)%levels(n)) )
      allocate( lon_inv (iv%info(gpseph)%levels(n)) )
      allocate( azim_inv(iv%info(gpseph)%levels(n)) )
      allocate( tmp_h   (iv%info(gpseph)%levels(n)) )

      nlevel = iv%info(gpseph)%levels(n)
      do k=1,nlevel
         tmp_h(k)    = iv%gpseph(n)%h(k)/1000.    ! observational height in km
      end do

      do k=kts,kte
         if (mean_h(k) > tmp_h(1)) goto 1
      end do

1     lowest_level = k

      do k=kts,kte
       new_ob%height(k)=mean_h(k)
      enddo

! interpolate_level =2000, tmp_z(k) increase with constant interval(0.01 km) in the vertical direction
      do k=1,interpolate_level
      tmp_z(k) = (iv%gpseph(n)%h(1))/1000.0+0.01*(k-1)  ! iv%gpseph(n)%h(1) is in meter
      end do

      do k=1,nlevel
      log_ref(k)  = log(iv%gpseph(n)%ref(k)%inv)
      lat_inv(k)  = iv%gpseph(n)%lat(k)%inv
      lon_inv(k)  = iv%gpseph(n)%lon(k)%inv
!     if ( use_nonlocal ) then
!     azim_inv(k) = iv%gpseph(n)%azim(k)%inv - 180.0
!     else
      azim_inv(k) = iv%gpseph(n)%azim(k)%inv
!     endif
!     tmp_h(k)    = iv%gpseph(n)%h(k)/1000.    ! observational height in km
      end do

!! interpolate variables and smooth them on the standard height grid
      CALL lintp( nlevel, tmp_h, log_ref,               &
                  interpolate_level, tmp_z, tmp_ref )
      CALL lintp( nlevel, tmp_h, lat_inv,               &
                  interpolate_level, tmp_z, tmp_lat )
      CALL lintp( nlevel, tmp_h, lon_inv,               &
                  interpolate_level, tmp_z, tmp_lon )
      CALL lintp( nlevel, tmp_h, azim_inv,              &
                  interpolate_level, tmp_z, tmp_azim )      

      do i = lowest_level+1, kte-1

       do j = 1, interpolate_level
       if (tmp_z(j) < (mean_h(i)+mean_h(i-1))/2.d0) nstart = j
       if (tmp_z(j) < (mean_h(i)+mean_h(i+1))/2.d0) nstop  = j
       end do 

      new_ob%ref(i) = 0.d0
      new_ob%lat(i) = 0.d0
      new_ob%lon(i) = 0.d0
      new_ob%azim(i)= 0.d0

       do j = nstart, nstop
       new_ob%ref(i)  = new_ob%ref(i)  + tmp_ref(j)
       new_ob%lat(i)  = new_ob%lat(i)  + tmp_lat(j)
       new_ob%lon(i)  = new_ob%lon(i)  + tmp_lon(j)
       new_ob%azim(i) = new_ob%azim(i) + tmp_azim(j) 
       end do

      new_ob%ref(i)  = new_ob%ref(i) / (nstop-nstart+1)
      new_ob%lat(i)  = new_ob%lat(i) / (nstop-nstart+1)
      new_ob%lon(i)  = new_ob%lon(i) / (nstop-nstart+1)
      new_ob%azim(i) = new_ob%azim(i)/ (nstop-nstart+1)
      end do

      do k = lowest_level+1, kte-1
      new_ob%ref(k) = exp(new_ob%ref(k))
      end do
 
      deallocate (log_ref, lat_inv, lon_inv, azim_inv, tmp_h ) 

       dtr=pi/180.
       rtd=180./pi
       step=5.
       is = ids
       ie = ide
       js = jds
       je = jde
       ks = kds
       ke = kde

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
       gps_rays(n)%nbot = nbot
       gps_rays(n)%ntop = ntop


       allocate ( gps_rays(n)%je2 (nbot:ntop) )
       allocate ( gps_rays(n)%ip123 (nbot:ntop) )
       allocate ( gps_rays(n)%ilocal (nbot:ntop) )

       do i=nbot,ntop
!
!* rp is a vector along tangent point direction
!* rn is a vector normal to rp along the meridian
!* transform coordiante from sphere to cartesian coordinate 
!* A is an angle from vector to x-y plane, 
!* B is an angle from x-axis to the project of the vector on x-y plane    
!* x=r*cos(A)*cos(B)
!* y=r*cos(A)*sin(B)
!* z=r*sin(A)
       rp(1)=(new_ob%rfict+mean_h(i))*cos(dtr*new_ob%lat(i)) & 
                                     *cos(dtr*new_ob%lon(i))
       rp(2)=(new_ob%rfict+mean_h(i))*cos(dtr*new_ob%lat(i)) &
                                     *sin(dtr*new_ob%lon(i))
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
       gps_rays(n)%je2(i)=int(ray_length/step)+1

       allocate ( gps_rays(n)%ip123(i)%i1(    -je2:je2) )
       allocate ( gps_rays(n)%ip123(i)%i2(    -je2:je2) )
       allocate ( gps_rays(n)%ip123(i)%i3(    -je2:je2) )
       allocate ( gps_rays(n)%ip123(i)%h (    -je2:je2) )
       allocate ( gps_rays(n)%ip123(i)%w1(2,  -je2:je2) )
       allocate ( gps_rays(n)%ip123(i)%w2(2,2,-je2:je2) )

       gps_rays(n)%ip123(i)%h(0) = h

       ip1=0
       ip2=0
       ip3=0
       CALL find2(ie,je,rxabs,x2,y2,glat,glon,ip1,ip2,w2)
       gps_rays(n)%ip123(i)%i1(0) = ip1 
       gps_rays(n)%ip123(i)%i2(0) = ip2 
       gps_rays(n)%ip123(i)%w2(1:2,1:2,0) = w2

       CALL find1(ke,mean_h,h,ip3,w1)
       gps_rays(n)%ip123(i)%i3(0) = ip3 
       gps_rays(n)%ip123(i)%w1(1:2,0) = w1

       i1l=0
       i2l=0
       i1r=0
       i2r=0
       i1= 0
       i2= 0
       
      DO j=1,je2
!* calculate S from TP point and integrate to different direction (WRF)
!* S is asymmetric
       do k=-1,+1,2
       dst=step*j*k
        do l=1,3
        rx(l)=rp(l)+dst*rr(l)
        end do
!* transform coordiante from cartesian to sphere coordinate
       rox=sqrt(rx(1)**2+rx(2)**2)
       glat=rtd*atan(rx(3)/rox)
       glon=rtd*atan2(rx(2),rx(1))
       CALL absv(rx,rxabs)
!       h = rxabs - earth_radius
       h = rxabs - new_ob%rfict
       gps_rays(n)%ip123(i)%h(k*j) = h
       
!        IF (h > mean_h(ke)) goto 2
       IF (h <= mean_h(ke-1))then
!       IF (h <= mean_h(27)) then
!       CALL find1(ke,mean_h,h,i3,w1)
        if (k == -1) then
        i1=i1l
        i2=i2l
        endif
        if (k == 1) then
        i1=i1r
        i2=i2r
        endif
       CALL find2(ie,je,rxabs,x2,y2,glat,glon,i1,i2,w2)
       gps_rays(n)%ip123(i)%i1(k*j) = i1 
       gps_rays(n)%ip123(i)%i2(k*j) = i2 
       gps_rays(n)%ip123(i)%w2(1:2,1:2,k*j) = w2

! avoid boundry effect 
        if (i1==1 .or. i1>=ie-1 .or. i2==1 .or. i2>=je-1) then
        w2(:,:)=0.
        gps_rays(n)%ip123(i)%w2(1:2,1:2,k*j) = w2
        endif
! avoid mountain effect
       CALL find1(ke,mean_h,h,i3,w1)
       gps_rays(n)%ip123(i)%i3(k*j) = i3 
       gps_rays(n)%ip123(i)%w1(1:2,k*j) = w1
       if (h < (global_terr(i1,i2)/1000.)) then
       !print*,'h<terr',i,j,k,i1,i2,xb%terr(i1,i2)/1000.
       w1(:)=0.
       gps_rays(n)%ip123(i)%w1(1:2,k*j) = w1
       endif

        if (k == -1) then
        i1l=i1
        i2l=i2
        endif
        if (k == 1) then
        i1r=i1
        i2r=i2
        endif


      END IF
       end do ! from left to right

      END DO !from the begining to the end

      end do

!
END SUBROUTINE da_gpseph_rays

! ................................................................

SUBROUTINE spin ( v1, vs, a, v2 )
! *****************************************************
!  S.V.Sokolovskiy
! ***************************************************
!  This subroutine rotates vector V1 around vector VS
!  at angle A. V2 is the vector after the rotation.
! ***************************************************
       implicit  none
       integer :: i, j, k
       real ::  vsabs, a1, a2, a3 
       real,intent(in)::a
       real, dimension(3),intent(in)::v1
       real, dimension(3),intent(in)::vs
       real, dimension(3)::vsn
       real, dimension(3),intent(out)::v2
       real, dimension(3,3) :: s

! Calculation of the unit vector around which
! the rotation should be done.

       vsabs=sqrt(vs(1)**2+vs(2)**2+vs(3)**2)
       do 1 i=1,3
1      vsn(i)=vs(i)/vsabs

!-------------------------------------------------------------------------------------
! Calculation of the rotation matrix.
! Any linear transformation of a vector A can be 
! represented via a matrix M:
!    -           -   
!   | a11 a12 a13 | 
! A=| a21 a22 a23 |
!   | a31 a32 a33 | 
!    -           -         
!              -     -                 -                 -             -          -
!             | 1 0 0 |               | c1**2 c1*c2 c1*c3 |           |  0 -c3  c2 |
!  =cos(delta)| 0 1 0 |+(1-cos(delta))| c2*c1 c2**2 c2*c3 |+sin(delta)|  c3  0 -c1 |
!             | 0 0 1 |               | c3*c1 c3*c2 c3**2 |           | -c2  c1  0 |
!              -     -                 -                 -             -          -
!
!--------------------------------------------------------------------------------------
       a1=cos(a)
       a2=1.-a1
       a3=sin(a)
       s(1,1)=a2*vsn(1)*vsn(1)+a1
       s(1,2)=a2*vsn(1)*vsn(2)-a3*vsn(3)
       s(1,3)=a2*vsn(1)*vsn(3)+a3*vsn(2)
       s(2,1)=a2*vsn(2)*vsn(1)+a3*vsn(3)
       s(2,2)=a2*vsn(2)*vsn(2)+a1
       s(2,3)=a2*vsn(2)*vsn(3)-a3*vsn(1)
       s(3,1)=a2*vsn(3)*vsn(1)-a3*vsn(2)
       s(3,2)=a2*vsn(3)*vsn(2)+a3*vsn(1)
       s(3,3)=a2*vsn(3)*vsn(3)+a1

! Calculation of the rotated vector.

       do 2 i=1,3
2      v2(i)=s(i,1)*v1(1)+s(i,2)*v1(2)+s(i,3)*v1(3)

END SUBROUTINE spin

! .............................................................

SUBROUTINE rnorm ( v )
       implicit none
       integer :: i, j, k
       real :: s
       real, dimension(3),intent(inout):: v

       s=sqrt(v(1)**2+v(2)**2+v(3)**2)
       do i=1,3
       v(i)=v(i)/s
       end do
END SUBROUTINE rnorm

! ...............................................................

SUBROUTINE absv ( r, rabs )
!*****************************************************************
!* The absolute value of a vector
!*****************************************
       implicit none
       real,intent(out)::rabs
       real,dimension(3),intent(in):: r
 
       rabs=sqrt(r(1)**2+r(2)**2+r(3)**2)

END SUBROUTINE absv

! ..............................................................

SUBROUTINE find1 ( n, x, xx, n1, w )
       implicit none
       integer :: i, j, k, n2
       integer,intent(in)::n
       integer,intent(out):: n1
       real :: dst, dst1, ww
       real,intent(in)::xx
       real, dimension(n),intent(in)::x
       real, dimension(2),intent(out)::w
!  n:  the number of total model levels
!=====================================================================
!       dst=1.0e+38
!       do i=1,n-1
!
!       dst1=(x(i)-xx)**2+(x(i+1)-xx)**2
!       if (dst1.lt.dst) then
!       dst=dst1
!       n1=i
!       end if
!       end do
!
!       w(1)=abs(x(n1+1)-xx)
!       w(2)=abs(x(n1)-xx)
!       ww=w(1)+w(2)
!       w(1)=w(1)/ww
!       w(2)=w(2)/ww
!=====================================================================
      do i=1,n-1
      if ((xx.ge.x(i)).and.(xx.le.x(i+1))) then
      n1=i
      n2=i+1
      ww=x(n2)-x(n1)
      w(1)=(x(n2)-xx)/ww
      w(2)=(xx-x(n1))/ww
      goto 1
      end if
      end do
1     continue
      if (xx.lt.x(1)) then
!  for xx outside the model bottom 
      n1= 1
      w(1)= 1.0
      w(2)= 0.0
      endif
      if (xx.gt.x(n)) then
!  for xx outside the model top 
      n1= n-1
      w(1)= 0.0   ! weighting on n1 level
      w(2)= 1.0   ! weighting on n1+1 level (i.e., the top level)
      endif
!=====================================================================

END SUBROUTINE find1
! ...........................................................................

SUBROUTINE find2 ( ide, jde, rfict, x2, y2, xp, yp, i1, j1, w )
      IMPLICIT NONE
      integer :: i, j, k, l, m, is, ie, js, je
      integer,intent(inout)::i1, j1
      integer :: i2,j2   ! for test finding points
      integer :: ide, jde
      real :: dtr, dst, dst1, dst2, tmp, ww
      real :: xpd,ypd,xps,yps,zps,x2d,y2d,x2s,y2s,z2s,xdis,ydis,zdis
      real,intent(in)::xp,yp,rfict
      real, dimension(ide,jde),intent(in):: x2, y2
      real, dimension(2,2),intent(out) :: w
!
       dtr=pi/180.

!  transfer model xlat, xlon to spherical x and y
!

      if((i1 > 0) .and. (j1 > 0))then   

       is= max0(i1-3,1)
       ie= min0(i1+3,ide-1)
       js= max0(j1-3,1)
       je= min0(j1+3,jde-1)

      else

       is=1
       ie=ide-1
       js=1
       je=jde-1

      end if

       dst=1.e29
       do i=is,ie
       do j=js,je

       dst1=0.0
       do l=i,i+1
       do m=j,j+1
       x2d= x2(l,m)-xp
       y2d= y2(l,m)-yp
       dst1=dst1+(x2d*x2d+y2d*y2d)
       end do
       end do
       if (dst1 < dst) then
       dst=dst1
       i1=i
       j1=j
       end if

       end do
       end do

!**********************************************************
! for test finding points
!
!       dst=1.e29
!       do i=1,xb%ide-1
!       do j=1,xb%jde-1
!
!       dst1=0.0
!       do l=i,i+1
!       do m=j,j+1
!       x2d= x2(l,m)-xp2
!       y2d= y2(l,m)-yp2
!       dst1=dst1+(x2d*x2d+y2d*y2d)
!       end do
!       end do
!       if (dst1 < dst) then
!       dst=dst1
!       i2=i
!       j2=j
!       end if
!
!       end do
!       end do
!
!       if ((i1 /= i2).or.(j1 /= j2))then
!       write(225,*) 'test for finding points'
!       write(225,*) 'i1,i2=',i1,i2
!       write(225,*) 'j1,j2=',j1,j2
!       endif
!**********************************************************
!  calculate weighting 
       xpd= xp*dtr
       ypd= yp*dtr
       xps= rfict*cos(xpd)*cos(ypd)
       yps= rfict*cos(xpd)*sin(ypd)
       zps= rfict*sin(xpd)
       do i=1,2
       do j=1,2
       x2d= x2(i1+i-1,j1+j-1)*dtr
       y2d= y2(i1+i-1,j1+j-1)*dtr
       x2s= rfict*cos(x2d)*cos(y2d)
       y2s= rfict*cos(x2d)*sin(y2d)
       z2s= rfict*sin(x2d)
       xdis= x2s-xps
       ydis= y2s-yps
       zdis= z2s-zps  
       tmp= xdis*xdis + ydis*ydis + zdis*zdis
       w(i,j)=1.0/sqrt(tmp)
       end do
       end do


       ww= 0.
       do i=1,2
       do j=1,2
       ww= ww+w(i,j)
       end do
       end do

       do i=1,2
       do j=1,2
       w(i,j)= w(i,j)/ww
       end do
       end do

END SUBROUTINE find2

!--------------------------------------------------------------------
!  DA_Residual_Gpseph.inc for all subroutines which are needed 
!  by nonlocal operator calculation
!                                   Shu-Ya Chen      Nov. 2005
!--------------------------------------------------------------------

SUBROUTINE model_mean_altitude ( xb, new_ob,  mean_h )  
      implicit none

! input : xb%h(:,:,:)
! output : mean_h(:) in kilometer

      type (xb_type), intent(in)      :: xb       ! first guess state.
      type (ob_in_mean_h), intent(in) :: new_ob
      real, dimension(kts:kte),intent(out) :: mean_h  ! mean altitude in each level of model domain

      integer :: i, j, k, is, ie, js, je, ks, ke
      real :: temp 
      real :: earth_radius_in_meter, horizontal_domain
      real, dimension(ids:ide,jds:jde,kds:kde) :: geo_h

      earth_radius_in_meter = new_ob%rfict*1000.
      is = ids
      ie = ide
      js = jds
      je = jde
      ks = kds
      ke = kde
! z is normalized potential height (m)= Phi/g0
      do k=ks,ke
      do j=js,je
      do i=is,ie
!     geo_h(i,j,k)=earth_radius_in_meter*xb%h(i,j,k)/  &   ! xb%h(i,j,k) is nomalized height in m
!                 (earth_radius_in_meter-xb%h(i,j,k))      ! geo_h is geometric height in m
      geo_h(i,j,k)=global_h(i,j,k)
      end do
      end do
      end do
!
! Calculation of the mean altitude grid
      horizontal_domain=float((ie-is+1)*(je-js+1))
      do k=ks,ke
         temp=0.0
         do j=js,je
         do i=is,ie
          temp=temp+geo_h(i,j,k)
         end do
         end do
         mean_h(k)=temp/horizontal_domain          ! in m
         mean_h(k)=mean_h(k)/1000.                 ! in km
      enddo

END SUBROUTINE model_mean_altitude

