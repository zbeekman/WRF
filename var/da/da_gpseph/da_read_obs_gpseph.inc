SUBROUTINE da_read_obs_gpseph( n, xb, iv, lowest_level, new_ob, mean_h, ilocal) 

      IMPLICIT NONE
!-----------------------------------------------------------------------
! For calculating the Excess Phase from observation refertivity        |
!                                                                      |
! input : n, xb%h, xb%dimensions, iv%gpseph(n)%info%levels,            |
!         iv%gpseph(n)%ref(k)%inv, iv%gpseph(n)%lat(k)%inv,            |
!         iv%gpseph(n)%lon(k)%inv, iv%gpseph(n)%azim(k)%inv,           |
!         iv%gpseph(n)%h(k)                                            |
! output: lowest_level, new_ob%height, new_ob%ref,                     |
!         new_ob%lat, new_ob%lon, new_ob%azim, mean_h                  |
!                                                                      |
!         by Shu-Ya Chen         (Nov.2005)                            |
!-----------------------------------------------------------------------

      type (xb_type), intent(in)      :: xb       ! first guess state.
      type (iv_type), intent(in)      :: iv       ! O-B structure.
      type (ob_in_mean_h), intent(out):: new_ob   ! Obs. data in the 
                                                  ! mean altitude of model

      integer, intent(in)             :: n   
      integer, intent(out)            :: lowest_level  ! the first altitude of observation 
                                                       ! which related to model levels
      real, dimension(kts:kte),intent(out) :: mean_h  ! mean altitude in each level of model domain
      integer, dimension(kts:kte), intent(out):: ilocal

      integer                         :: i, j, k 
      integer                         :: nlevel        ! the count of observational levels
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
! --------------------------------------------------------------------------      
   deallocate(new_ob%height,new_ob%ref,new_ob%lat,new_ob%lon,new_ob%azim,new_ob%eph)
   allocate (new_ob%height(kds:kde))
   allocate (new_ob%ref(kds:kde))
   allocate (new_ob%lat(kds:kde))
   allocate (new_ob%lon(kds:kde))
   allocate (new_ob%azim(kds:kde))
   allocate (new_ob%eph(kds:kde))
   new_ob%height(:) = 0.0
   new_ob%ref(:) = 0.0
   new_ob%lat(:) = 0.0
   new_ob%lon(:) = 0.0
   new_ob%azim(:) = 0.0
   new_ob%eph(:) = 0.0

!      new_ob%rfict=iv%gpseph(n)%info%rfict
      temp_lat=(iv%info(gpseph)%lat(1,n)*pi)/180.
      new_ob%rfict=SQRT((requator**2*rpolar**2)/                    &
                        (requator**2*(sin(abs(temp_lat)))**2+  &
                         rpolar**2*(cos(abs(temp_lat)))**2)) 
      CALL model_mean_altitude(xb, new_ob, mean_h)          ! mean_h in km
     

      do k=kts,kte
       new_ob%height(k)=mean_h(k)
      enddo

      allocate( log_ref (iv%info(gpseph)%levels(n)) )
      allocate( lat_inv (iv%info(gpseph)%levels(n)) )
      allocate( lon_inv (iv%info(gpseph)%levels(n)) )
      allocate( azim_inv(iv%info(gpseph)%levels(n)) )
      allocate( tmp_h   (iv%info(gpseph)%levels(n)) )

! interpolate_level =2000, tmp_z(k) increase with constant interval(0.01 km) in the vertical direction
      do k=1,interpolate_level
      tmp_z(k) = (iv%gpseph(n)%h(1))/1000.0+0.01*(k-1)  ! iv%gpseph(n)%h(1) is in meter
      end do

      nlevel = iv%info(gpseph)%levels(n)
      do k=1,nlevel
      log_ref(k)  = log(iv%gpseph(n)%ref(k)%inv)
      lat_inv(k)  = iv%gpseph(n)%lat(k)%inv
      lon_inv(k)  = iv%gpseph(n)%lon(k)%inv
!     if ( use_nonlocal ) then
!     azim_inv(k) = iv%gpseph(n)%azim(k)%inv - 180.0
!     else
      azim_inv(k) = iv%gpseph(n)%azim(k)%inv
!     endif
      tmp_h(k)    = iv%gpseph(n)%h(k)/1000.    ! observational height in km
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


      do k=kts,kte
      if (mean_h(k) > tmp_h(1)) goto 1
      end do

1     lowest_level = k

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

 
      CALL obs_nonlocal(n, xb,mean_h,lowest_level,new_ob,ilocal)

      deallocate (log_ref, lat_inv, lon_inv, azim_inv, tmp_h ) 
     

END SUBROUTINE da_read_obs_gpseph

