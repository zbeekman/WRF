subroutine da_read_mod_gpseph_adj(num, jo_grad_y, xb, new_ob, &
                                  lowest_level, jo_grad_x, ilocal)

      implicit none

      type (y_type),  intent(in)       :: jo_grad_y   
      type (x_type),  intent(inout)    :: jo_grad_x   
      type (xb_type), intent(in)       :: xb       ! first guess state.
      type (ob_in_mean_h), intent(in)  :: new_ob   ! Obs. data in the mean altitude of model
      integer, intent(in)              :: lowest_level  ! the first altitude of observation which related to model levels
      integer , intent(in) :: num
      real, dimension(kts:kte)  :: mean_h  ! mean altitude in each level of model domain
      integer, dimension(kts:kte), intent(in) :: ilocal

      INTEGER :: i, j, k, is, ie, js, je, ks, ke
      INTEGER :: nlevel        ! the count of observational levels
      REAL    :: horizontal_domain
      REAL    :: earth_radius_in_meter
      REAL    :: temp
      REAL    :: dtr
      REAL    :: qq9, ee9
      REAL, DIMENSION(kts:kte)   :: tmp_refm_adj
      REAL, DIMENSION(kts:kte)   :: tmp_ref_adj
      REAL, DIMENSION(kts:kte)   :: tmp_z
      REAL, DIMENSION(kts:kte)   :: tmp_ref_fw
      REAL, DIMENSION(kts:kte)   :: tmp_z_fw
      REAL, DIMENSION(kts:kte)   :: tmp_ref
      REAL, DIMENSION(kts:kte)   :: refm
      REAL, DIMENSION(kts:kte)   :: tmp_refm
      REAL, DIMENSION(kts:kte)   :: tmp_refm_fw
      REAL, DIMENSION(3,kts:kte) :: cc, a_cc
      REAL, DIMENSION(ids:ide,jds:jde,kds:kde) :: geo_h9
      REAL, DIMENSION(ids:ide,jds:jde,kds:kde) :: ref_keep
      REAL, DIMENSION(ids:ide,jds:jde,kds:kde) ::ref_mean_h
      REAL, DIMENSION(ids:ide,jds:jde,kds:kde) ::ref_mean_h_adj
      REAL, DIMENSION(ids:ide,jds:jde,kds:kde) ::refm_keep
      REAL, DIMENSION(kts:kte)   :: mean_h_meter 
      real    :: partone, parttwo, dividnd
      real    :: partone9,parttwo9,dividnd9
! ----------------------------------------------------------------------------      
      is=ids
      ie=ide
      js=jds
      je=jde
      ks=kds
      ke=kde
   

      dtr = pi/180.
!      earth_radius_in_meter = earth_radius *1000.0
      earth_radius_in_meter = new_ob%rfict *1000.0
!      jo_grad_x%ref(:,:,:) = 0.0
      a_cc(:,:) = 0.0
      tmp_ref_adj(:) = 0.0
      tmp_refm_adj(:) = 0.0
     
      do k=ks,ke
      do j=js,je
      do i=is,ie
!     geo_h9(i,j,k)=earth_radius_in_meter*xb%h(i,j,k)/  & 
!                 (earth_radius_in_meter-xb%h(i,j,k)) 
      geo_h9(i,j,k)=global_h(i,j,k)
      enddo
      enddo
      enddo

      horizontal_domain = float((ie-is+1)*(je-js+1))
      do k=ks,ke
         temp=0.0
         do j=js,je
         do i=is,ie
          temp=temp+geo_h9(i,j,k)
         end do
         end do
         mean_h_meter(k)=temp/horizontal_domain   ! in m
      enddo

! THE INPUT mean_h IS IN KILOMETER
!      ref_keep(:,:,:)= xb%ref(:,:,:) 
!      mean_h_keep(:) = mean_h(:)*1000.     ! in meter
!      mean_h_keep(:) = mean_h(:)     ! in meter
      
! Interpolation of refractivity on the mean altitude grid (basic)

      do j=js,je
      do i=is,ie

       do k=ks,ke
       tmp_z_fw(k)=geo_h9(i,j,k)
       tmp_ref_fw(k)=log(global_ref(i,j,k))
       end do

       CALL splinx(ke,tmp_z_fw,tmp_ref_fw,ke,mean_h_meter,cc,tmp_refm_fw)

       do k=ks,ke
       refm_keep(i,j,k) = tmp_refm_fw(k)
       ref_mean_h(i,j,k)=exp(tmp_refm_fw(k))
       end do

      end do
      end do
 
      do k=ks,ke
       mean_h(k)=mean_h_meter(k)/1000.   ! in km
      end do

      call mod_nonlocal_adj(num,jo_grad_y,xb,new_ob,mean_h, &
                            lowest_level,ref_mean_h_adj,ilocal)
     
!      xb%ref(:,:,:) = ref_keep(:,:,:)
!      mean_h(:) = mean_h_keep(:)


    do j =  je,js,-1
    do i =  ie,is,-1
      do k = ke,ks,-1
       refm(k)=refm_keep(i,j,k)
        tmp_refm_adj(k) = tmp_refm_adj(k)+ref_mean_h_adj(i,j,k)*exp(refm(k))
        ref_mean_h_adj(i,j,k) = 0.0
      end do

      do k=ks,ke
      tmp_z(k)=geo_h9(i,j,k)
      tmp_ref(k)=log(global_ref(i,j,k))
      end do

      call splinx_adj( ke,tmp_z,tmp_ref,tmp_ref_adj,ke,mean_h_meter,cc,a_cc,tmp_refm_adj )
      
      do k =ke,ks,-1
!      jo_grad_x%ref(i,j,k) = jo_grad_x%ref(i,j,k)+tmp_ref_adj(k)*(1./global_ref(i,j,k))
       global_adj_ref(i,j,k) = global_adj_ref(i,j,k)+tmp_ref_adj(k)*(1./global_ref(i,j,k))
       tmp_ref_adj(k) = 0.
      end do
    end do
    end do
    
end SUBROUTINE da_read_mod_gpseph_adj
