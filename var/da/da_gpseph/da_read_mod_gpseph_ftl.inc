subroutine da_read_mod_gpseph_ftl(n, xa, xb, iv, lowest_level, new_ob, y)
      implicit none
      type (xb_type), intent(in)       :: xb       ! first guess state.
    ! type (x_type), intent(inout)     :: xa      ! gridded analysis increment. 
      type (x_type), intent(in)        :: xa      ! gridded analysis increment. 
      type (iv_type),intent(in)        :: iv       ! O-B structure.
      type (y_type),intent(inout)      ::  y
      type (ob_in_mean_h),intent(in) :: new_ob   ! Obs. data in the mean altitude of model
      integer,intent(in) :: lowest_level  ! the first altitude of observation which related to model levels
      integer, intent(in) :: n   
      integer :: i, j, k, is, ie, js, je, ks, ke
      integer :: nlevel        ! the count of observational levels
      real :: horizontal_domain
      real :: earth_radius_in_meter
      real :: temp
      real :: qq, ee, qq9, ee9
      real, dimension(kts:kte)   :: mean_h  ! mean altitude in each level of model domain
      real, dimension(kts:kte)   :: tmp_z
      real, dimension(kts:kte)   :: tmp_ref,tmp_ref9
      real, dimension(kts:kte)   :: refm,refm9
      real, dimension(3,kts:kte) :: cc,cc9
      real, dimension(its:ite,jts:jte,kts:kte) :: geo_h
      real, dimension(its:ite,jts:jte,kts:kte) :: ref_mean_h
      real, dimension(its:ite,jts:jte,kts:kte) :: ref_mean_h_ftl     
      real    :: partone, parttwo, dividnd
      real    :: partone9,parttwo9,dividnd9

!      earth_radius_in_meter = earth_radius *1000.0
      earth_radius_in_meter = new_ob%rfict *1000.0
      is = its
      ie = ite
      js = jts
      je = jte
      ks = kts
      ke = kte
      
      do i=is,ie
      do j=js,je
      do k=ks,ke
!      partone  = 0.776*xa%p(i,j,k)/xb%t(i,j,k) &
!               - 0.776*xb%p(i,j,k)*xa%t (i,j,k)/(xb%t(i,j,k)*xb%t(i,j,k))
!      partone9 = 0.776*xb%p(i,j,k)/xb%t(i,j,k)
!      dividnd  = xa%t (i,j,k)*(0.622+0.378*xb%q(i,j,k)) &
!               + xb%t(i,j,k)*       0.378*xa%q(i,j,k)
!      dividnd9 = xb%t(i,j,k)*(0.622+0.378*xb%q(i,j,k))
!      parttwo  =     coeff*xa%q(i,j,k)/dividnd9 &
!               -     coeff*xb%q(i,j,k)*dividnd /(dividnd9*dividnd9)
!      parttwo9 = 1.0+coeff*xb%q(i,j,k)/dividnd9
!      xa%ref(i,j,k)= partone9 * parttwo &
!                   + partone  * parttwo9
!      xb%ref(i,j,k)=partone9 * parttwo9
!     geo_h(i,j,k)=earth_radius_in_meter*xb%h(i,j,k)/  &   ! xb%h(i,j,k) is nomalized height in m
!                 (earth_radius_in_meter-xb%h(i,j,k))      ! geo_h is geometric height in m
      geo_h(i,j,k)=xb%h(i,j,k)
      end do
      end do
      end do

      
      horizontal_domain=float((ie-is+1)*(je-js+1))

      do k=ks,ke
         temp=0.0
         do j=js,je
         do i=is,ie
          temp=temp+geo_h(i,j,k)
         end do
         end do
         mean_h(k)=temp/horizontal_domain   ! in m
      enddo
      do i=is,ie
      do j=js,je

       do k=ks,ke
       tmp_z(k)=geo_h(i,j,k)
       tmp_ref(k)=xa%ref(i,j,k)*(1./xb%ref(i,j,k))
       tmp_ref9(k)=log(xb%ref(i,j,k))
       enddo

      call splinx_ftl(ke,tmp_z,tmp_ref9,tmp_ref,cc9,cc,ke,mean_h,refm9,refm)
       do k=ks,ke
       ref_mean_h_ftl(i,j,k)=refm(k)*exp(refm9(k))
       ref_mean_h(i,j,k)=exp(refm9(k))
       end do

      enddo
      enddo
    
      do k=ks,ke
       mean_h(k)=mean_h(k)/1000.   ! in km
      end do

      CALL mod_nonlocal_ftl(n, xb, new_ob, mean_h, lowest_level,  &
                            ref_mean_h_ftl, y)
      
END SUBROUTINE da_read_mod_gpseph_ftl
