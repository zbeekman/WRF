subroutine da_print_stats_gpseph(stats_unit, ngpseph, GPSEph)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   integer,                  intent(in)    :: stats_unit
   integer,                  intent(inout) :: ngpseph
   type (stats_gpseph_type), intent(in)    :: gpseph

   if (trace_use_dull) call da_trace_entry("da_print_stats_gpseph")
   
   write (unit=stats_unit, fmt='(a/)') '   var           eph(m)  n    k'  

   write (unit=stats_unit, fmt='(a,i16)') '  Number: ', ngpseph

   if (ngpseph < 1) ngpseph = 1
   
   write(unit=stats_unit, fmt='((a,f12.4,2i5))') &
      ' Minimum(n,k): ', GPSEph%minimum%eph,    &
      ' Maximum(n,k): ', GPSEph%maximum%eph
   write(unit=stats_unit, fmt='((a,f12.4,10x))') &
      ' Average     : ', GPSEph%average%eph/real(ngpseph),    &
      '    RMSE     : ', sqrt(GPSEph%rms_err%eph/real(ngpseph))

   if (trace_use_dull) call da_trace_exit("da_print_stats_gpseph")

end subroutine da_print_stats_gpseph


