subroutine da_residual_gpseph(iv, y, re, np_missing, np_bad_data, np_obs_used, np_available)

   !-----------------------------------------------------------------------
   ! Purpose: Calculate residiual for gpseph obs
   !-----------------------------------------------------------------------

   implicit none

   type (iv_type),intent(inout) :: iv     ! Innovation vector (O-B).
   type (y_type) ,intent(in)    :: y      ! y = H (xa)
   type (y_type) ,intent(inout) :: re     ! Residual vector (O-A).

   integer       ,intent(inout) :: np_available
   integer       ,intent(inout) :: np_obs_used
   integer       ,intent(inout) :: np_missing
   integer       ,intent(inout) :: np_bad_data

   type (bad_data_type) :: n_obs_bad
   integer              :: n, k
   integer              :: nlevels2
!
   real                 :: constant0, g_lat, g_height, weight1, weight2, &
                           gpseph_org

   if (trace_use_dull) call da_trace_entry("da_residual_gpseph")

! Assuming the highest latitude is 90.0 degree:
         constant0= sin(radian * 90.0)

   n_obs_bad % gpseph % num = number_type(0, 0, 0)

   do n=1, iv%info(gpseph)%nlocal
      nlevels2 = max(kde,iv%info(gpseph)%levels(n))
   !  do k=1, iv%info(gpseph)%levels(n)
      do k=1, nlevels2
         np_available = np_available + 1
!
! Weighted the GPSREF innovation with the latitude: 
         if (iv%gpseph(n)%eph(k)%qc >= obs_qc_pointer ) then

! depend on the height: above 7km, set to 1.0, below 7km, decrease to 0.0:
           g_height = iv%gpseph(n)% h(k)
           weight1 = 1.0 - (7000.0 - g_height) / 7000.0
           if ( weight1 > 1.0 ) weight1 = 1.0
! not depend on height:
           weight1 = 1.0

! depend on the latitude, at 90 degree, weight = 1.0, at 0 degree, weight = 0.0
           g_lat    = iv%info(gpseph)%lat(k,n)
           weight2  = abs(sin(radian * g_lat) / constant0)
! not depend on the latitude:
           weight2   = 1.0

           gpseph_org = iv%gpseph(n)%eph(k)%inv
           iv%gpseph(n)%eph(k)%inv = gpseph_org * weight1 * weight2
         endif
!.............................................................
         if (iv%gpseph(n)%eph(k)%qc >= obs_qc_pointer ) then
         re%gpseph(n)%eph(k) = &
            da_residual(n, k, y%gpseph(n)%eph(k), iv%gpseph(n)%eph(k), n_obs_bad%gpseph)
         endif
!
         if (iv%gpseph(n)%eph(k)%qc >= obs_qc_pointer ) &
           iv%gpseph(n)%eph(k)%inv = gpseph_org
         
      end do
   end do

   np_missing  = np_missing  + n_obs_bad % gpseph % num % miss
   np_bad_data = np_bad_data + n_obs_bad % gpseph % num % bad
   np_obs_used = np_obs_used + n_obs_bad % gpseph % num % use

   if (trace_use_dull) call da_trace_exit("da_residual_gpseph")

   
end subroutine da_residual_gpseph


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

!....................................................................................
SUBROUTINE lintp ( n1, x1, y1, n2, x2, y2 )

! **********************************************************
!  This subroutine provides linear interpolation of function
!  y(x) from one grid x1(i), onto another grid x2(i)
! 
!  y1(i), i=1,n1 are input grid function values
!  y2(i), i=1,n2 are output grid function values
! **********************************************************
!  both x1(i) and x2(i) should increase with 'i'
! **********************************************************
      implicit none
      integer :: i, j, k, l, m, n
      integer, intent(in)::n1
      integer, intent(in)::n2
      real,dimension(n1),intent(in)::x1
      real,dimension(n2),intent(in)::x2
      real,dimension(n1),intent(in)::y1
      real,dimension(n2),intent(out)::y2
      

      do i=1,n2
      if (x2(i).le.x1(1)) then
      y2(i)=y1(1)+(x2(i)-x1(1))*(y1(2)-y1(1))/(x1(2)-x1(1))
      else if (x2(i).ge.x1(n1)) then
      y2(i)=y1(n1)+(x2(i)-x1(n1))*(y1(n1)-y1(n1-1))/(x1(n1)-x1(n1-1))
      else

      do j=1,n1
      k=j
      if ((x1(j).le.x2(i)).and.(x1(j+1).ge.x2(i))) goto 1
      end do

1     y2(i)=y1(k)+(x2(i)-x1(k))*(y1(k+1)-y1(k))/(x1(k+1)-x1(k))
     

      end if
      end do

!      return
 
END SUBROUTINE lintp

!.........................................................................
SUBROUTINE obs_nonlocal ( xb, mean_h, lowest_level, new_ob, ilocal)

      implicit none

      type (xb_type), intent(in)     :: xb       ! first guess state.
      type (ob_in_mean_h), intent(inout)   :: new_ob
       
      real,dimension(kts:kte),intent(in) :: mean_h
      integer,intent(in) :: lowest_level
      integer :: i, j, k, l, m, n, i1, i2, i3, nbot, ntop
      integer :: is, ie, js, je, ks, ke
      integer :: i1l,i2l,i1r,i2r
      real :: dtr, rtd, step, dst, rox, glat, glon, h, rxabs, tmp_ref
      real :: xlatd, xlond
      real,dimension(ids:ide,jds:jde) :: x2, y2
      real,dimension(kts:kte) :: refm
      real,dimension(3):: rp, rn, rr, rx
      real,dimension(2) :: w1
      real,dimension(2,2) :: w2  
      integer,dimension(kts:kte),intent(out) :: ilocal
      integer :: ip1,ip2,ip3,je2
      real :: refp, ray_length

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

       ip1=0
       ip2=0
       ip3=0
       CALL find2(ie,je,rxabs,x2,y2,glat,glon,ip1,ip2,w2)

       CALL find1(ke,mean_h,h,ip3,w1)
       refp=0.
        do l=1,2
        do m=1,2
        do n=1,2
        refp=refp+new_ob%ref(ip3+n-1)*w2(l,m)*w1(n)
        end do
        end do
        end do
       
       ilocal(i)=1
       new_ob%eph(i)=0.

       i1l=0
       i2l=0
       i1r=0
       i2r=0
       i1= 0
       i2= 0
       
      
!       j=0
!1      j=j+1
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

! avoid boundry effect 
        if (i1==1 .or. i1>=ie-1 .or. i2==1 .or. i2>=je-1) then
        w2(:,:)=0.
        endif
! avoid mountain effect
       CALL find1(ke,mean_h,h,i3,w1)
       if (h < (global_terr(i1,i2)/1000.)) then
       !print*,'h<terr',i,j,k,i1,i2,xb%terr(i1,i2)/1000.
       w1(:)=0.
       endif

        if (k == -1) then
        i1l=i1
        i2l=i2
        endif
        if (k == 1) then
        i1r=i1
        i2r=i2
        endif


       tmp_ref=0.

        do l=1,2
        do m=1,2
        do n=1,2
        tmp_ref=tmp_ref+new_ob%ref(i3+n-1)*w2(l,m)*w1(n)
        end do
        end do
        end do
       new_ob%eph(i)=new_ob%eph(i)+step*tmp_ref

      END IF
       end do ! from left to right

!       goto 1

      END DO !from the begining to the end

      if ( .not. gpseph_local  ) then 
        if (new_ob%eph(i) == 0.) then
        new_ob%eph(i)=step*refp
        ilocal(i)=0
        endif
      else if ( gpseph_local ) then
       ilocal(i)=0
       new_ob%eph(i)=step*refp
      endif 

      end do
!
END SUBROUTINE obs_nonlocal

! ................................................................

SUBROUTINE splinx ( n, x, y, nn, xx, c, yy )

!*************************************************************
! This subroutine provides natural cubic spline interpolation
! of a given grid function onto another given grid.
! Inputs:  n - number of input grid points
!          x - array of input argument data
!          y - array of input function data
!          nn - number of output data
!          xx - array of output argument data
! Outputs: c - array of spline coefficients
!          yy - array of output function data
!*************************************************************
! If output grid is outside input grid
! then linear extrapolation is applied
!*************************************************************


      implicit none
      integer :: i, j, k, m
      integer,intent(in):: n 
      integer,intent(in):: nn
      real :: dxr, dyr, dxl, dyl, dx, dy, del, der, d, dyx
      real, dimension(n)  :: x,y
      real, dimension(nn),intent(in) ::xx
      real, dimension(nn),intent(out)::yy
      real, dimension(3,n),intent(out):: c

      c(1,1)=0.
      c(2,1)=0.
      c(2,n)=0.
      dxr=x(2)-x(1)
      dyr=(y(2)-y(1))/dxr
      do 1 i=2,n-1
      dxl=dxr
      dxr=x(i+1)-x(i)
      dyl=dyr
      dyr=(y(i+1)-y(i))/dxr
      dx=dxr+dxl
      dy=(dyr-dyl)/dx
      c(1,i)=-dxr/(2.*dx+dxl*c(1,i-1))
      c(2,i)=(6.*dx*dy-dxl*c(2,i-1))/(2.*dx+dxl*c(1,i-1))
1     continue
      do 2 i=n-1,2,-1
      c(2,i)=c(1,i)*c(2,i+1)+c(2,i)
2     continue
      do 3 i=1,n-1
      dx=x(i+1)-x(i)
      dy=y(i+1)-y(i)
      dyx=dy/dx
      c(1,i)=dyx-dx*(c(2,i)/3.+c(2,i+1)/6.)
      c(2,i)=c(2,i)/2.
      c(3,i)=(dy-c(1,i)*dx-c(2,i)*dx**2)/dx**3
3     continue

      do 5 j=1,nn
      if (xx(j).le.x(1)) then
      yy(j)=y(1)+c(1,1)*(xx(j)-x(1))
      else if (xx(j).ge.x(n)) then
      del=x(n)-x(n-1)
      der=c(1,n-1)+2.d0*c(2,n-1)*del+3.d0*c(3,n-1)*del**2
      yy(j)=y(n)+der*(xx(j)-x(n))
      else
      do 4 i=1,n
4     if (x(i).le.xx(j)) m=i
      d=xx(j)-x(m)
      yy(j)=y(m)+d*(c(1,m)+d*(c(2,m)+d*c(3,m)))
      end if
5     continue

END SUBROUTINE splinx

! ..............................................................
SUBROUTINE mod_nonlocal ( xb, new_ob, mean_h, lowest_level, ref_mean_h, model_eph)

      implicit none
! input : ref_mean_h, new_ob%lat, new_ob%lon, lowest_level, mean_h
! output : model_eph

      type (xb_type), intent(in)     :: xb       ! first guess state.
      type (ob_in_mean_h), intent(in)   :: new_ob
       
      real,dimension(kts:kte),intent(in) :: mean_h
      integer,intent(in) :: lowest_level
      real,dimension(ids:ide,jds:jde,kds:kde),intent(in)::ref_mean_h
      real,dimension(kts:kte),intent(out) :: model_eph

      integer :: i, j, k, l, m, n, i1, i2, i3, nbot, ntop
      integer :: is, ie, js, je, ks, ke
      integer :: i1l,i2l,i1r,i2r
      real :: dtr, rtd, step, dst, rox, glat, glon, h, rxabs, tmp_ref
      real :: xlatd, xlond
      real,dimension(ids:ide,jds:jde) :: x2, y2
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
       refp=0.
        do l=1,2
        do m=1,2
        do n=1,2
        refp=refp+ref_mean_h(ip1+l-1,ip2+m-1,ip3+n-1)*w2(l,m)*w1(n)
        end do
        end do
        end do
       ilocal(i)=1

       model_eph(i)=0.

       i1l=0
       i2l=0
       i1r=0
       i2r=0
       i1= 0
       i2= 0
      
!       j=0
!1      j=j+1
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

! avoid boundry effect 
        if (i1==1 .or. i1>=ie-1 .or. i2==1 .or. i2>=je-1) then
        w2(:,:)=0.
        endif
! avoid mountain effect
       CALL find1(ke,mean_h,h,i3,w1)
       if (h < (global_terr(i1,i2)/1000.)) then
       w1(:)=0.
       endif

        if (k == -1) then
        i1l=i1
        i2l=i2
        endif
        if (k == 1) then
        i1r=i1
        i2r=i2
        endif
       tmp_ref=0.
        do l=1,2
        do m=1,2
        do n=1,2
        tmp_ref=tmp_ref+ref_mean_h(i1+l-1,i2+m-1,i3+n-1)*w2(l,m)*w1(n)
        end do
        end do
        end do
       model_eph(i)=model_eph(i)+step*tmp_ref

      END IF
       end do

!       goto 1

      END DO
      if ( .not. gpseph_local ) then
        if (model_eph(i) == 0.) then
        model_eph(i)=step*refp
        ilocal(i)=0
        endif
      else if ( gpseph_local ) then
        ilocal(i)=0
        model_eph(i)=step*refp
      endif

      end do
!
END SUBROUTINE mod_nonlocal

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
      integer,intent(out)::i1, j1
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


