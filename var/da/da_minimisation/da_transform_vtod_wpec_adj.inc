subroutine da_transform_vtod_wpec_adj(cv_size, be, ep, cv, vp, vv, vp6, vv6, xbx, &
   grid, config_flags, jcdf_flag)

   !-------------------------------------------------------------------------
   ! Purpose:  Does Adjoint of control variable (V) transform to Obs-space(Y)
   !-------------------------------------------------------------------------

   implicit none

   integer,                    intent(in)    :: cv_size ! Size of cv array.
   type(be_type),              intent(in)    :: be     ! background error structure.
   type(ep_type),              intent(in)    :: ep     ! ensemble perturbation structure.
   real,                       intent(out)   :: cv(1:cv_size) ! control variables.
   type(vp_type),              intent(inout) :: vp     ! Grdipt/level CV.
   type(vp_type),              intent(inout) :: vp6    ! Grdipt/level CV for 6h.
   type(vp_type),              intent(inout) :: vv     ! Grdipt/EOF CV.
   type(vp_type),              intent(inout) :: vv6    ! Grdipt/EOF CV for 6h.
   type(xbx_type),             intent(inout) :: xbx    ! For header & non-grid arrays.
   type(domain),               intent(inout) :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags
   logical,                    intent(in)    :: jcdf_flag       ! additional flag to switch off JcDF, used to switch off JcDF in calculation of J.

   type(x_type) :: shuffle
   integer :: nobwin,ndynopt, jl_start, jl_end, time_step_seconds
   character*256 :: timestr, timestr1
   integer       :: i, j, k
   real          :: subarea, whole_area
   real          :: f_flag
   real, dimension(:,:,:), allocatable :: a_hr_rainc, a_hr_rainnc
   real, dimension(:,:,:), allocatable :: a_rainc, a_rainnc

   character(len=4) :: filnam

   call da_trace_entry("da_transform_vtod_wpec_adj")

   cv = 0.0

      call da_zero_x(grid%xa)

      call da_wpec_constraint_adj(grid, xbx)

      call da_transform_xtoxa_adj_all(grid)

      call da_transform_vtox_adj(grid, cv_size, xbx, be, ep, vp, vv, cv)

   if (trace_use) call da_trace_exit("da_transform_vtod_wpec_adj")

end subroutine da_transform_vtod_wpec_adj


