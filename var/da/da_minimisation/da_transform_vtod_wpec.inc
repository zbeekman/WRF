subroutine da_transform_vtod_wpec(cv_size, be, ep, cv, vp, vv, vp6, vv6, xbx, grid, config_flags)

   !----------------------------------------------------------------------
   ! Purpose:  Transform control variable (V) to observation-space (Y)
   !----------------------------------------------------------------------

   implicit none

   integer,                    intent(in)    :: cv_size ! Size of cv array.
   type(be_type),              intent(in)    :: be     ! background error structure.
   type(ep_type),              intent(in)    :: ep     ! Ensemble perturbation structure.
   real,                       intent(in)    :: cv(1:cv_size)     ! control variables.
   type(vp_type),              intent(inout) :: vp     ! Grdipt/level CV.
   type(vp_type),              intent(inout) :: vp6     ! Grdipt/level CV for 6h.
   type(vp_type),              intent(inout) :: vv     ! Grdipt/EOF CV.
   type(vp_type),              intent(inout) :: vv6     ! Grdipt/EOF CV for 6h.
   type(xbx_type),             intent(inout) :: xbx    ! For header & non-grid arrays.
   type(domain),               intent(inout) :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags

   type(x_type) :: shuffle
   integer :: nobwin, jl_start, jl_end, time_step_seconds

   character(len=4) :: filnam
   character(len=256) :: timestr, timestr1
   real, dimension(:,:,:), allocatable :: hr_rainc, hr_rainnc
   real, dimension(:,:,:), allocatable :: g_rainc, g_rainnc

   if (trace_use) call da_trace_entry("da_transform_vtod_wpec")

      call da_transform_vtox(grid, cv_size, xbx, be, ep, cv, vv, vp)
      call da_transform_xtoxa_all(grid)

      call da_wpec_constraint_lin(grid, xbx)


   if (trace_use) call da_trace_exit("da_transform_vtod_wpec")

end subroutine da_transform_vtod_wpec


