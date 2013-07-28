subroutine da_get_innov_vector_gpseph(it, num_qcstat_conv, grid, ob, iv)
!==========================================================================
!  Purpose: Calculate invovation (O-B) for gpseph
!       Updated for Analysis on Arakawa-C grid
!       Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
!               Shu-Ya Chen/Shu-Hua Chen        Date: 01/27/2011
!=========================================================================
   implicit none

   integer,          intent(in)    :: it       ! External iteration.
   type(domain),     intent(in)    :: grid     ! first guess state.
   type(y_type),     intent(inout) :: ob       ! Observation structure.
   type(iv_type),    intent(inout) :: iv       ! O-B structure.
   TYPE (ob_in_mean_h)             :: new_ob
   integer,          intent(inout) :: num_qcstat_conv(:,:,:,:)       

   integer :: n        ! Loop counter.
   integer :: i, j, k  ! Index dimension.
   INTEGER :: num_levs ! Number of obs levels.
   real    :: dx, dxm  ! Interpolation weights.
   real    :: dy, dym  ! Interpolation weights.
   real, dimension(:,:), allocatable :: model_eph !Model gpseph at ob loc
   real, dimension(:,:), allocatable :: model_ref !Model gpsref at ob loc
   real, dimension(:,:), allocatable :: iv_ref    !obs gpsref at ob loc
   real, dimension(:,:), allocatable :: cc
   REAL, DIMENSION(kds:kde)          :: mean_h
   integer, DIMENSION(kds:kde)       :: ilocal
   real, dimension(:), allocatable   ::err_height,err_pencent,intp_err_pencent
   integer :: kk, count_lev, err_lev
   integer :: lowest_level, nlevels, nlevels2
   integer :: Iu_eph, l, ntotal
   character(len=2) :: c_it

   if (trace_use_dull) call da_trace_entry("da_get_innov_vector_gpseph")
!=============================================================================
   ! GPS EPH Pseudo OBS test:

   if ( (pseudo_var(1:3) == 'eph' .and. num_pseudo > 0) .and.  &
        it == 1 ) then

      ! Deallocate:
      if (iv%info(gpseph)%nlocal > 0) then
         do n = 1, iv%info(gpseph)%nlocal
            deallocate (iv % gpseph(n) %  h )
            deallocate (iv % gpseph(n) % eph)
            deallocate (iv % gpseph(n) % ref)
            deallocate (iv % gpseph(n) %   p)
            deallocate (iv % gpseph(n) %   t)
            deallocate (iv % gpseph(n) %   q)
            deallocate (iv % gpseph(n) % azim)
            deallocate (iv % gpseph(n) % lat)
            deallocate (iv % gpseph(n) % lon)
            deallocate (ob % gpseph(n) % eph)
         end do
         deallocate(iv % gpseph)
         deallocate(ob % gpseph)
      end if
      use_gpsephobs = .true.
 
      iv%info(gpseph)%nlocal = num_pseudo
      iv%info(gpseph)%plocal(1) = num_pseudo
      iv%info(gpseph)%ntotal   = num_pseudo
      iv%info(gpseph)%max_lev = 1
      iv%info(pseudo)%nlocal = 0

      call da_allocate_observations (iv)
      iv%info(gpseph)%n1 = 1
      iv%info(gpseph)%n2 = 1

      allocate(iv%gpseph(num_pseudo)%  h(1:1))
      allocate(iv%gpseph(num_pseudo)%eph(1:1))
      allocate(iv%gpseph(num_pseudo)%ref(1:1))
      allocate(iv%gpseph(num_pseudo)%  h(1:1))
      allocate(iv%gpseph(num_pseudo)%  p(1:1))
      allocate(iv%gpseph(num_pseudo)%  t(1:1))
      allocate(iv%gpseph(num_pseudo)%  q(1:1))
      allocate(iv%gpseph(num_pseudo)%azim(1:1))
      allocate(iv%gpseph(num_pseudo)%lat(1:1))
      allocate(iv%gpseph(num_pseudo)%lon(1:1))
      allocate(ob%gpseph(1:num_pseudo))
      allocate(ob%gpseph(num_pseudo)%eph(1:1))

      write(stdout,'(a,i2)') '==> GPS REF pseudo OBS test: num_pseudo=',num_pseudo

      iv%info(gpseph)%levels(1) = 1

      iv%info(gpseph)%x(:,1) = pseudo_x
      iv%info(gpseph)%y(:,1) = pseudo_y

      iv%info(gpseph)%i(:,1) = int(pseudo_x)
      iv%info(gpseph)%j(:,1) = int(pseudo_y)
      iv % gpseph(1) %  h(1) = pseudo_z

      iv%info(gpseph)%dx(:,1) = pseudo_x-real(iv%info(gpseph)%i(1,1))
      iv%info(gpseph)%dy(:,1) = pseudo_y-real(iv%info(gpseph)%j(1,1))
      iv%info(gpseph)%dxm(:,1)=1.0-iv%info(gpseph)%dx(1,1)
      iv%info(gpseph)%dym(:,1)=1.0-iv%info(gpseph)%dy(1,1)

      iv % gpseph(1) %eph(1) % inv = pseudo_val
      iv % gpseph(1) %eph(1) % error = pseudo_err
      iv % gpseph(1) %eph(1) % qc = 0
      ! Set halo:
      if ((iv%info(gpseph)%i(1,1) < its-1) .or.(iv%info(gpseph)%i(1,1) > ite) .or. &
          (iv%info(gpseph)%j(1,1) < jts-1) .or.(iv%info(gpseph)%j(1,1) > jte)) then
         iv%info(gpseph)%proc_domain(:,1) = .false.
      else
         iv%info(gpseph)%proc_domain(:,1) = .false.

         if (iv%info(gpseph)%i(1,1) >= its .and. iv%info(gpseph)%i(1,1) <= ite .and. &
             iv%info(gpseph)%j(1,1) >= jts .and. iv%info(gpseph)%j(1,1) <= jte) then
            iv%info(gpseph)%proc_domain(:,1) = .true.
         end if
      end if

      write(stdout,'(a4,2f15.5)') pseudo_var, pseudo_val, pseudo_err
      write(stdout,'(3f15.5)')    pseudo_x, pseudo_y, pseudo_z
   end if
!=============================================================================

   if (iv%info(gpseph)%nlocal < 1) return

   ntotal = 0
   nlevels=max(kde,iv%info(gpseph)%max_lev)
   allocate (model_eph(nlevels,iv%info(gpseph)%n1:iv%info(gpseph)%n2))
   model_eph(:,:) = 0.0

   if (.not.(pseudo_var(1:3) == 'eph' .and. num_pseudo > 0)) then
      do n=iv%info(gpseph)%n1,iv%info(gpseph)%n2
         do k=1, iv%info(gpseph)%levels(n)
            if( iv%gpseph(n)%eph(k)%qc == fails_error_max .and. it > 1 ) &
                iv%gpseph(n)%eph(k)%qc = 0
         end do
         ! Get cross pt. horizontal interpolation weights:
         i   = iv%info(gpseph)%i(1,n)
         j   = iv%info(gpseph)%j(1,n)
         dx  = iv%info(gpseph)%dx(1,n)
         dy  = iv%info(gpseph)%dy(1,n)
         dxm = iv%info(gpseph)%dxm(1,n)
         dym = iv%info(gpseph)%dym(1,n)
      end do
   endif

! Before gross error check, fill in eph data into iv and ob

   if ( (.not.(pseudo_var(1:3) == 'eph' .and. num_pseudo > 0)) .or. &
        it > 1 ) then
   if ( .not. allocated ( gps_rays ) ) &
      allocate ( gps_rays(iv%info(gpseph)%n1:iv%info(gpseph)%n2) )
   do n=iv%info(gpseph)%n1,iv%info(gpseph)%n2
          if ( iv%info(gpseph)%levels(n) < 1 ) cycle
     
          allocate (new_ob%height(kds:kde))
          allocate (new_ob%ref(kds:kde))
          allocate (new_ob%lat(kds:kde))
          allocate (new_ob%lon(kds:kde))
          allocate (new_ob%azim(kds:kde))
          allocate (new_ob%eph(kds:kde))
          new_ob%height(:) = 0.0
          new_ob%   ref(:) = 0.0
          new_ob%   lat(:) = 0.0
          new_ob%   lon(:) = 0.0
          new_ob%  azim(:) = 0.0
          new_ob%   eph(:) = 0.0
   CALL da_gpseph_rays(n, grid%xb, iv, mean_h, lowest_level, new_ob)
   CALL da_read_obs_gpseph(n, grid%xb, iv, lowest_level, new_ob, mean_h, ilocal)
   CALL da_read_mod_gpseph(n, grid%xb, lowest_level, new_ob, model_eph(:,n))
       do k=kds,kde
       ob%gpseph(n)%eph(k)=new_ob%eph(k)
       enddo
       if ( iv%info(gpseph)%levels(n) > kde ) then
       do k=kde+1,iv%info(gpseph)%levels(n)
       ob%gpseph(n)%eph(k)= missing_r
       enddo
       endif
          deallocate(new_ob%height,new_ob%ref,new_ob%lat,new_ob%lon,new_ob%azim,new_ob%eph)

     iv%gpseph(n)%level1=lowest_level+1
     iv%gpseph(n)%level2=kde-1
     if (iv%gpseph(n)%level1 .gt. iv%gpseph(n)%level2) goto 100

! set quality control (QC) ---------------------------------------------------
   nlevels2=max(kde,iv%info(gpseph)%levels(n))
        do k = 1, nlevels2
           if ((k >= iv%gpseph(n)%level1) .and. (k <= iv%gpseph(n)%level2)) then
              iv%gpseph(n)%eph(k)%qc = 0
           else
              iv%gpseph(n)%eph(k)%qc = missing_data
           endif
        enddo

        do k = 1, nlevels2
           iv%gpseph(n)%eph(k)%inv = 0.0
           if (ob%gpseph(n)%eph(k) > missing_r .AND. &
               iv%gpseph(n)%eph(k)%qc >= obs_qc_pointer) then
           iv%gpseph(n)%eph(k)%inv = ob%gpseph(n)%eph(k) - model_eph(k,n)
           end if
        end do
 100    continue
   end do
   else
     do n=iv%info(gpseph)%n1,iv%info(gpseph)%n2
        ob % gpseph(1)%eph(1) = model_eph(1,n) + iv %gpseph(1)%eph(1)%inv 
     end do
   end if

   if (.not.(pseudo_var(1:3) == 'eph' .and. num_pseudo > 0) ) then

! -------------------------------------------------------------------
! [3.0] read obs. error
!  one month statistics (2003/08/15~2003/09/15) and
!  use Hollingsworth-Lonngberg to calculate the observation error
!   OBSERR = percentage(%) * Excess Phase
!                                           Shu-Ya Chen     Nov. 2006
!-------------------------------------------------------------------
!     open(1001,file='OBSERROR.TBL',status='old')
!     allocate(err_pencent(num_levs))
!        do kk=1,num_levs
!           err_pencent(kk)=missing_r
!        enddo
!     count_lev=0
!     read(1001,*)
!     do kk=1,10000
!        read(1001,*,END=1005) err_lev,err_height(kk), err_pencent(kk)
!        err_pencent(kk)=err_pencent(kk)*0.01
!        count_lev=count_lev+1
!     enddo
!1005 close(1001)
      do n=iv%info(gpseph)%n1,iv%info(gpseph)%n2
         num_levs = iv%info(gpseph)%levels(n)
         allocate(err_height(num_levs))
         allocate(err_pencent(num_levs))
         do kk=1,num_levs
            err_pencent(kk)=missing_r
         enddo
         count_lev=0
         open(1001,file='OBSERROR.TBL',status='old')
         read(1001,*)
         do kk=1,num_levs
            read(1001,*,END=1005) err_lev,err_height(kk), err_pencent(kk)
            err_pencent(kk)=err_pencent(kk)*0.01
            count_lev=count_lev+1
         enddo
 1005    close(1001)
         allocate(intp_err_pencent(kde))
         allocate(cc(3,count_lev))
!        write(1006,'(2(A7),A12,3(A20))') 'number','level','mean_h','error(%)','Excess Phase','OBS ERROR VALSE'
!        call lintp(count_lev,err_height,err_pencent,kde,mean_h,intp_err_pencent)
         call splinx(count_lev,err_height,err_pencent,kde,mean_h,cc,intp_err_pencent)
         do kk=kds,kde
         iv%gpseph(n)%eph(kk)%error = intp_err_pencent(kk)*ob%gpseph(n)%eph(kk)
         write(1006,'(2(I7),F12.5,3(F20.7))') n,kk,mean_h(kk),intp_err_pencent(kk),ob%gpseph(n)%eph(kk),iv%gpseph(n)%eph(kk)%error
         enddo
         if (iv%info(gpseph)%levels(n) > kde) then
         do k=kde+1,num_levs
            iv%gpseph(n)%eph(k)%error = missing_r
         enddo
         endif

         deallocate(intp_err_pencent,err_height,err_pencent,cc)
      enddo
   endif

!   ! Quality check 1: Gross error(departure from the background) check 

   if (.not.(pseudo_var(1:3) == 'eph' .and. num_pseudo > 0) .and. check_max_iv ) &
      call da_check_max_iv_gpseph(iv, it, num_qcstat_conv)

! ------------------------------------------------------------------------------
   if (.not. anal_type_verify) then
   if (.not.(pseudo_var(1:3) == 'eph' .and. num_pseudo > 0)) then

  goto 12345

! Write out GPS Eph data:
     DO n=iv%info(gpseph)%n1,iv%info(gpseph)%n2
     Iu_eph = 337
     write(c_it,'(I2.2)') it
     open (unit=Iu_eph, file='RO_Innov_'//iv%info(gpseph)%date_char(n)//'_'//c_it, &
           form='formatted')
           write(unit=Iu_eph, fmt='(/i5,2x,a,2x,a,2x,4f10.3,i5)') n, &
              iv%info(gpseph)%date_char(n), iv%info(gpseph)%id(n),  &
              iv%info(gpseph)%lat(1,n)  , iv%info(gpseph)%lon(1,n), &
              iv%info(gpseph)%x(1,n)  , iv%info(gpseph)%y(1,n), &
              iv%info(gpseph)%levels(n)
           write(unit=Iu_eph, fmt='(a5,3x,6a14)') 'level','     height   ', &
                       '    Obs_eph   ','  model_eph   ','  Innov_eph   ', &
                       '  error_eph   ',' qc_eph       '
        !  do k = 1, iv%info(gpseph)%levels(n)
        !  do k = 1, nlevels
           do k = kds,kde
               write(unit=Iu_eph, fmt='(i5,1x,5f14.3,i10)')  k, &
        !         iv%gpseph(n)%h(k), ob%gpseph(n)%eph(k),   &
                  mean_h(k)*1000., ob%gpseph(n)%eph(k),   &
                  model_eph(k,n), iv%gpseph(n)%eph(k)%inv,   &
                  iv%gpseph(n)%eph(k)%error, iv%gpseph(n)%eph(k)%qc
           end do
     close(Iu_eph)
     END DO
12345 continue
   ! .........................................................................
   end if  ! end of pseudo test
   end if  ! end of verify check

   deallocate (model_eph)

   if (trace_use_dull) call da_trace_exit("da_get_innov_vector_gpseph")

end subroutine da_get_innov_vector_gpseph
