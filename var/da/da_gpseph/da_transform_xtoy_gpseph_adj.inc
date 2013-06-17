subroutine da_transform_xtoy_gpseph_adj(grid, iv, jo_grad_y, jo_grad_x)

   !-------------------------------------------------------------------------
   ! Purpose: TBD
   !-------------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(in)    :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)
   type (ob_in_mean_h)           :: new_ob
   integer                       :: n  ! Loop counter.
   INTEGER                       :: lowest_level
   REAL, DIMENSION(kds:kde)      :: mean_h
   integer, DIMENSION(kds:kde)   :: ilocal

   if (trace_use_dull) call da_trace_entry("da_transform_xtoy_gpseph_adj")

! Have to allocate new_ob (different part with da_get_innov_vector_Gpseph.inc)
!
   do n=iv%info(gpseph)%n1,iv%info(gpseph)%n2
     allocate (new_ob%height(kds:kde))
     allocate (new_ob%ref(kds:kde)   )
     allocate (new_ob%lat(kds:kde)   )
     allocate (new_ob%lon(kds:kde)   )
     allocate (new_ob%azim(kds:kde)  )
     allocate (new_ob%eph(kds:kde)   )
     new_ob%height(:) = 0.0
     new_ob%ref(:)    = 0.0
     new_ob%lat(:)    = 0.0
     new_ob%lon(:)    = 0.0
     new_ob%azim(:)   = 0.0
     new_ob%eph(:)    = 0.0
     call da_read_obs_gpseph( n, grid%xb, iv, lowest_level, new_ob, mean_h, ilocal)
     call da_read_mod_gpseph_adj(n, jo_grad_y, grid%xb, new_ob, lowest_level, jo_grad_x, ilocal)

     deallocate (new_ob%height)
     deallocate (new_ob%ref   )
     deallocate (new_ob%lat   )
     deallocate (new_ob%lon   )
     deallocate (new_ob%azim  )
     deallocate (new_ob%eph   )
   end do

   if (trace_use_dull) call da_trace_exit("da_transform_xtoy_gpseph_adj")

end subroutine da_transform_xtoy_gpseph_adj


