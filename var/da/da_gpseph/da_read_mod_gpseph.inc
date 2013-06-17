! this program is for reading model N to calculate the excess phase
!                                     Shu-Ya Chen   Nov. 2005               
! 
SUBROUTINE da_read_mod_gpseph ( n, xb, lowest_level, new_ob, model_eph)
      IMPLICIT NONE

      TYPE (xb_type), INTENT(IN)       :: xb       ! first guess state.
      TYPE (ob_in_mean_h), INTENT(IN)  :: new_ob   ! Obs. data in the mean altitude of model
      integer, intent(in) :: n
      integer, intent(in) :: lowest_level  ! the first altitude of observation which related to model levels
      real, dimension(kts:kte),intent(out)   :: model_eph

      integer :: i, j, k, is, ie, js, je, ks, ke
      integer :: nlevel        ! the count of observational levels
      real :: horizontal_domain
      real :: earth_radius_in_meter 
      real :: temp
      real :: qq, ee
      real, dimension(kts:kte)   :: mean_h  ! mean altitude in each level of model domain
      real, dimension(kts:kte)   :: tmp_z
      real, dimension(kts:kte)   :: tmp_ref
      real, dimension(kts:kte)   :: refm
      real, dimension(3,kts:kte) :: cc
      real, dimension(its:ite,jts:jte,kts:kte) :: geo_h
      real, dimension(its:ite,jts:jte,kts:kte) :: ref_mean_h


!      earth_radius_in_meter = earth_radius *1000.
      earth_radius_in_meter = new_ob%rfict *1000.

      is = its
      ie = ite
      js = jts
      je = jte
      ks = kts
      ke = kte

      
!************************************************************
!variables from da_read_obs_gpseph: lowest_level,new_lat,new_lon,new_azim
!!! attention: if call this subroutine,
!!! one should change lowest_level into lowest_level+1 in the topper subroutine.
!************************************************************

! z is normalized potential height (m)= Phi/g0
! p is pressure (pa), t is temperature (k)
! q is specific humidity (kg/kg)
      do i=is,ie
      do j=js,je
      do k=ks,ke
!     geo_h(i,j,k)=earth_radius_in_meter*xb%h(i,j,k)/  &   ! xb%h(i,j,k) is nomalized height in m
!                 (earth_radius_in_meter-xb%h(i,j,k))      ! geo_h is geometric height in m
      geo_h(i,j,k)=xb%h(i,j,k)
! check if readed z is real geometric height...
      end do
      end do
      end do

!
! Calculation of the mean altitude grid
      horizontal_domain=float((ite-its+1)*(jte-jts+1))
      do k=ks,ke
         temp=0.0
         do j=js,je
         do i=is,ie
          temp=temp+geo_h(i,j,k)
         end do
         end do
         mean_h(k)=temp/horizontal_domain   ! in m
      enddo
! Interpolation of refractivity on the mean altitude grid

      do i=is,ie
      do j=js,je

       do k=ks,ke
       tmp_z(k)=geo_h(i,j,k)
       tmp_ref(k)=log(xb%ref(i,j,k))
       end do
       CALL splinx(ke,tmp_z,tmp_ref,ke,mean_h,cc,refm)
       do k=ks,ke
       ref_mean_h(i,j,k)=exp(refm(k))
       end do

      end do
      end do
     

      do k=ks,ke
       mean_h(k)=mean_h(k)/1000.   ! in km
      end do

      CALL mod_nonlocal(xb,new_ob,mean_h,lowest_level,ref_mean_h, model_eph)

       
END SUBROUTINE da_read_mod_gpseph 

