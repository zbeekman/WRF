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
SUBROUTINE obs_nonlocal ( n, xb, mean_h, lowest_level, new_ob, ilocal)

      implicit none

      integer, intent(in)             :: n
      type (xb_type), intent(in)     :: xb       ! first guess state.
      type (ob_in_mean_h), intent(inout)   :: new_ob
       
      real,dimension(kts:kte),intent(inout) :: mean_h
      integer,intent(in) :: lowest_level
      integer :: i, j, k, l, m, nn, i1, i2, i3, nbot, ntop
      integer :: is, ie, js, je, ks, ke
      integer :: i1l,i2l,i1r,i2r
      real :: dtr, rtd, step, dst, rox, glat, glon, h, rxabs, tmp_ref
      real :: xlatd, xlond
!zx   real,dimension(ids:ide,jds:jde) :: x2, y2
      real,dimension(kts:kte) :: refm
      real,dimension(3):: rp, rn, rr, rx
      real,dimension(2) :: w1
      real,dimension(2,2) :: w2  
      integer,dimension(kts:kte),intent(out) :: ilocal
      integer :: ip1,ip2,ip3,je2
      real :: refp, ray_length

!zx    dtr=pi/180.
!zx    rtd=180./pi
       step=5.
!zx    is = ids
!zx    ie = ide
!zx    js = jds
!zx    je = jde
       ks = kds
       ke = kde

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

!zx    nbot=lowest_level+1
       nbot = gps_rays(n)%nbot
!zx    ntop=ke-1
       ntop = gps_rays(n)%ntop

       mean_h = gps_rays(n)%mean_h

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
!zx    rp(1)=(new_ob%rfict+mean_h(i))*cos(dtr*new_ob%lat(i)) & 
!zx                                  *cos(dtr*new_ob%lon(i))
!zx    rp(2)=(new_ob%rfict+mean_h(i))*cos(dtr*new_ob%lat(i)) &
!zx                                  *sin(dtr*new_ob%lon(i))
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
       h = gps_rays(n)%ip123(i)%h(0)
!zx    ray_length=sqrt((mean_h(ke)-h)*(2*new_ob%rfict+mean_h(ke)+h))
!zx    je2=int(ray_length/step)+1
       je2 = gps_rays(n)%je2(i)

!zx    ip1=0
!zx    ip2=0
!zx    ip3=0
!zx    CALL find2(ie,je,rxabs,x2,y2,glat,glon,ip1,ip2,w2)
       ip1 = gps_rays(n)%ip123(i)%i1(0)
       ip2 = gps_rays(n)%ip123(i)%i2(0)
       w2 = gps_rays(n)%ip123(i)%w2(1:2,1:2,0)

!zx    CALL find1(ke,mean_h,h,ip3,w1)
       ip3 = gps_rays(n)%ip123(i)%i3(0)
       w1 = gps_rays(n)%ip123(i)%w1(1:2,0)
       refp=0.
        do l=1,2
        do m=1,2
        do nn=1,2
        refp=refp+new_ob%ref(ip3+nn-1)*w2(l,m)*w1(nn)
        end do
        end do
        end do
       
       ilocal(i)=1
       gps_rays(n)%ilocal(i)=1

       new_ob%eph(i)=0.

!zx    i1l=0
!zx    i2l=0
!zx    i1r=0
!zx    i2r=0
!zx    i1= 0
!zx    i2= 0
       
      
!       j=0
!1      j=j+1
      DO j=1,je2
!* calculate S from TP point and integrate to different direction (WRF)
!* S is asymmetric
       do k=-1,+1,2
!zx    dst=step*j*k
!zx     do l=1,3
!zx     rx(l)=rp(l)+dst*rr(l)
!zx     end do
!* transform coordiante from cartesian to sphere coordinate
!zx    rox=sqrt(rx(1)**2+rx(2)**2)
!zx    glat=rtd*atan(rx(3)/rox)
!zx    glon=rtd*atan2(rx(2),rx(1))
!zx    CALL absv(rx,rxabs)
!       h = rxabs - earth_radius
!zx    h = rxabs - new_ob%rfict
       h = gps_rays(n)%ip123(i)%h(k*j)
       
!        IF (h > mean_h(ke)) goto 2
       IF (h <= mean_h(ke-1))then
!       IF (h <= mean_h(27)) then
!       CALL find1(ke,mean_h,h,i3,w1)
!zx     if (k == -1) then
!zx     i1=i1l
!zx     i2=i2l
!zx     endif
!zx     if (k == 1) then
!zx     i1=i1r
!zx     i2=i2r
!zx     endif
!zx    CALL find2(ie,je,rxabs,x2,y2,glat,glon,i1,i2,w2)

! avoid boundry effect 
!zx     if (i1==1 .or. i1>=ie-1 .or. i2==1 .or. i2>=je-1) then
!zx     w2(:,:)=0.
!zx     endif
       i1 = gps_rays(n)%ip123(i)%i1(k*j)
       i2 = gps_rays(n)%ip123(i)%i2(k*j)
       w2 = gps_rays(n)%ip123(i)%w2(1:2,1:2,k*j)
! avoid mountain effect
!zx    CALL find1(ke,mean_h,h,i3,w1)
!zx    if (h < (global_terr(i1,i2)/1000.)) then
!zx    !print*,'h<terr',i,j,k,i1,i2,xb%terr(i1,i2)/1000.
!zx    w1(:)=0.
!zx    endif
       i3 = gps_rays(n)%ip123(i)%i3(k*j)
       w1 = gps_rays(n)%ip123(i)%w1(1:2,k*j)

!zx     if (k == -1) then
!zx     i1l=i1
!zx     i2l=i2
!zx     endif
!zx     if (k == 1) then
!zx     i1r=i1
!zx     i2r=i2
!zx     endif


       tmp_ref=0.

        do l=1,2
        do m=1,2
        do nn=1,2
        tmp_ref=tmp_ref+new_ob%ref(i3+nn-1)*w2(l,m)*w1(nn)
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
        gps_rays(n)%ilocal(i)=0
        endif
      else if ( gpseph_local ) then
        ilocal(i)=0
        gps_rays(n)%ilocal(i)=0
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
SUBROUTINE mod_nonlocal ( n, xb, new_ob, mean_h, lowest_level, ref_mean_h, model_eph)

      implicit none
! input : ref_mean_h, new_ob%lat, new_ob%lon, lowest_level, mean_h
! output : model_eph

      integer, intent(in)            :: n
      type (xb_type), intent(in)     :: xb       ! first guess state.
      type (ob_in_mean_h), intent(in)   :: new_ob
       
      real,dimension(kts:kte),intent(inout) :: mean_h
      integer,intent(in) :: lowest_level
      real,dimension(ids:ide,jds:jde,kds:kde),intent(in)::ref_mean_h
      real,dimension(kts:kte),intent(out) :: model_eph

      integer :: i, j, k, l, m, nn, i1, i2, i3, nbot, ntop
      integer :: is, ie, js, je, ks, ke
      integer :: i1l,i2l,i1r,i2r
      real :: dtr, rtd, step, dst, rox, glat, glon, h, rxabs, tmp_ref
      real :: xlatd, xlond
!zx   real,dimension(ids:ide,jds:jde) :: x2, y2
      real,dimension(kts:kte) :: refm
      real,dimension(3):: rp, rn, rr, rx
      real,dimension(2) :: w1
      real,dimension(2,2) :: w2  
      integer,dimension(kts:kte) :: ilocal
      integer :: ip1,ip2,ip3,je2
      real :: refp, ray_length

!zx    dtr=pi/180.
!zx    rtd=180./pi
       step=5.
!zx    is = ids
!zx    ie = ide
!zx    js = jds
!zx    je = jde
       ks = kds
       ke = kde

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
       nbot = gps_rays(n)%nbot
       ntop=ke-1
       ntop = gps_rays(n)%ntop

       mean_h = gps_rays(n)%mean_h

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
       h = gps_rays(n)%ip123(i)%h(0)
!zx    ray_length=sqrt((mean_h(ke)-h)*(2*new_ob%rfict+mean_h(ke)+h))
!zx    je2=int(ray_length/step)+1
       je2 = gps_rays(n)%je2(i)
!zx    ip1=0
!zx    ip2=0
!zx    ip3=0
!zx    CALL find2(ie,je,rxabs,x2,y2,glat,glon,ip1,ip2,w2)
       ip1 = gps_rays(n)%ip123(i)%i1(0)
       ip2 = gps_rays(n)%ip123(i)%i2(0)
       w2 = gps_rays(n)%ip123(i)%w2(1:2,1:2,0)

!zx    CALL find1(ke,mean_h,h,ip3,w1)
       ip3 = gps_rays(n)%ip123(i)%i3(0)
       w1 = gps_rays(n)%ip123(i)%w1(1:2,0)

       refp=0.
        do l=1,2
        do m=1,2
        do nn=1,2
        refp=refp+ref_mean_h(ip1+l-1,ip2+m-1,ip3+nn-1)*w2(l,m)*w1(nn)
        end do
        end do
        end do
       ilocal(i)=1

       model_eph(i)=0.

!zx    i1l=0
!zx    i2l=0
!zx    i1r=0
!zx    i2r=0
!zx    i1= 0
!zx    i2= 0
      
!       j=0
!1      j=j+1
      DO j=1,je2
!* calculate S from TP point and integrate to different direction (WRF)
!* S is asymmetric
       do k=-1,+1,2
!zx    dst=step*j*k
!zx     do l=1,3
!zx     rx(l)=rp(l)+dst*rr(l)
!zx     end do
!* transform coordiante from cartesian to sphere coordinate
!zx    rox=sqrt(rx(1)**2+rx(2)**2)
!zx    glat=rtd*atan(rx(3)/rox)
!zx    glon=rtd*atan2(rx(2),rx(1))
!zx    CALL absv(rx,rxabs)
!       h = rxabs - earth_radius
!zx    h = rxabs - new_ob%rfict
       h = gps_rays(n)%ip123(i)%h(k*j)

       IF (h <= mean_h(ke-1))then
!       IF (h <= mean_h(27)) then
!       CALL find1(ke,mean_h,h,i3,w1)
!zx     if (k == -1) then
!zx     i1=i1l
!zx     i2=i2l
!zx     endif
!zx     if (k == 1) then
!zx     i1=i1r
!zx     i2=i2r
!zx     endif
!zx    CALL find2(ie,je,rxabs,x2,y2,glat,glon,i1,i2,w2)

! avoid boundry effect 
!zx     if (i1==1 .or. i1>=ie-1 .or. i2==1 .or. i2>=je-1) then
!zx     w2(:,:)=0.
!zx     endif
       i1 = gps_rays(n)%ip123(i)%i1(k*j)
       i2 = gps_rays(n)%ip123(i)%i2(k*j)
       w2 = gps_rays(n)%ip123(i)%w2(1:2,1:2,k*j)

! avoid mountain effect
!zx    CALL find1(ke,mean_h,h,i3,w1)
!zx    if (h < (global_terr(i1,i2)/1000.)) then
!zx    w1(:)=0.
!zx    endif
       i3 = gps_rays(n)%ip123(i)%i3(k*j)
       w1 = gps_rays(n)%ip123(i)%w1(1:2,k*j)

!zx     if (k == -1) then
!zx     i1l=i1
!zx     i2l=i2
!zx     endif
!zx     if (k == 1) then
!zx     i1r=i1
!zx     i2r=i2
!zx     endif
       tmp_ref=0.
        do l=1,2
        do m=1,2
        do nn=1,2
        tmp_ref=tmp_ref+ref_mean_h(i1+l-1,i2+m-1,i3+nn-1)*w2(l,m)*w1(nn)
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

