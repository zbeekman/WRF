module da_gpseph

   use module_domain, only : domain, x_type, xb_type
   use module_dm, only : wrf_dm_sum_real

   use da_control, only : obs_qc_pointer,max_ob_levels,missing_r, &
      v_interp_p, v_interp_h, check_max_iv_print, pi, radian, interpolate_level, &
      missing, max_error_uv, max_error_t, rootproc,fails_error_max, &
      max_error_p,max_error_q, check_max_iv_unit,check_max_iv, qcstat_conv_unit, &
      max_stheight_diff,missing_data,max_error_bq,max_error_slp, ob_vars, &
      max_error_bt, max_error_buv, gpseph,max_error_thickness,gpseph_local, &
      pseudo_var, num_pseudo, kms,kme,kts,kte, trace_use_dull, &
      anal_type_verify,fails_error_max,pseudo_err,pseudo_x, pseudo_y, stdout, &
      use_gpsephobs, pseudo_z,pseudo_val,max_error_eph, &
      pseudo, jts, jte,its,ite,ids,ide,jds,jde,kds,kde, gpseph_nonlocal
   use da_define_structures, only : maxmin_type, iv_type, y_type, jo_type, &
      bad_data_type, number_type, bad_data_type, &
      maxmin_type, da_allocate_observations
   use da_interpolation, only : da_interp_lin_3d,da_interp_lin_3d_adj, &
      da_to_zk
   use da_par_util, only : da_proc_stats_combine
   use da_par_util1, only : da_proc_sum_int
   use da_statistics, only : da_stats_calculate
   use da_tools, only : da_max_error_qc, da_residual, da_convert_zk,da_get_print_lvl
   use da_tracing, only : da_trace_entry, da_trace_exit

   ! The "stats_gpseph_type" is ONLY used locally in da_gpseph:

   type residual_gpseph1_type
      real :: eph                   ! GPS linear excess phase
      real :: ref                   ! GPS Refractivity
      real :: lat
      real :: lon
      real :: azim
      real ::   p                   ! Retrieved from GPS Refractivity
      real ::   t                   ! Retrieved from GPS Refractivity
      real ::   q                   ! Used in GPS Refra. retrieval.
   end type residual_gpseph1_type

   type maxmin_gpseph_stats_type
      type (maxmin_type)         :: eph          ! GPS Refractivity
   end type maxmin_gpseph_stats_type

   type stats_gpseph_type
      type (maxmin_gpseph_stats_type)  :: maximum, minimum
      type (residual_gpseph1_type)     :: average, rms_err
   end type stats_gpseph_type

   TYPE ob_in_mean_h    ! calculate Obs. data in the mean altitude of model
      real                       :: rfict
      REAL             , pointer :: height(:)
      REAL             , pointer :: lat(:)   ! Latitude
      REAL             , pointer :: lon(:)   ! Longitude
      REAL             , pointer :: ref(:)   ! Refractivity
      REAL             , pointer :: azim(:)  ! Azimuth
      REAL             , pointer :: eph(:)   ! Excess Phase
   END TYPE ob_in_mean_h

contains

#include "da_ao_stats_gpseph.inc"
#include "da_calculate_grady_gpseph.inc"
#include "da_jo_and_grady_gpseph.inc"
#include "da_oi_stats_gpseph.inc"
#include "da_print_stats_gpseph.inc"
#include "da_transform_xtoy_gpseph.inc"
#include "da_transform_xtoy_gpseph_adj.inc"
#include "da_check_max_iv_gpseph.inc"
#include "da_get_innov_vector_gpseph.inc"
#include "da_read_mod_gpseph.inc"
#include "da_read_obs_gpseph.inc"
#include "da_residual_gpseph.inc"
#include "da_read_mod_gpseph_ftl.inc"
#include "da_read_mod_gpseph_adj.inc"
#include "da_residual_gpseph_ftl.inc"
#include "da_residual_gpseph_adj.inc"

end module da_gpseph

