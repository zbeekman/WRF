subroutine da_gpseph_final ()

   ! -------------------------------------------------------------------
   ! Purpose: Release global  lat, lon, terr, h
   ! -------------------------------------------------------------------

   implicit none

   if (trace_use_dull) call da_trace_entry("da_gpseph_final")

   if (gpseph_nonlocal)   then

      deallocate ( global_lat )
      deallocate ( global_lon )
      deallocate ( global_terr )
      deallocate ( global_h )
      deallocate ( global_ref )
      deallocate ( global_xa_ref )
      deallocate ( global_adj_ref )

   endif

   if (trace_use_dull) call da_trace_exit("da_gpseph_final")

end subroutine da_gpseph_final


