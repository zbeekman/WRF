subroutine da_gpseph_final ()

   ! -------------------------------------------------------------------
   ! Purpose: Release global  lat, lon, terr, h
   ! -------------------------------------------------------------------

   implicit none

   integer n1, n2, n, k

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

   if ( allocated (gps_rays ) ) then
      n1 = Lbound ( gps_rays, 1 )
      n2 = ubound ( gps_rays, 1 )
      do n = n1, n2
         if ( allocated (gps_rays(n)%je2) ) deallocate (gps_rays(n)%je2)
         if ( allocated (gps_rays(n)%mean_h) ) deallocate (gps_rays(n)%mean_h)
         do k = gps_rays(n)%nbot, gps_rays(n)%ntop
            deallocate ( gps_rays(n)%ip123(k)%i1 )
            deallocate ( gps_rays(n)%ip123(k)%i2 )
            deallocate ( gps_rays(n)%ip123(k)%i3 )
            deallocate ( gps_rays(n)%ip123(k)%w1 )
            deallocate ( gps_rays(n)%ip123(k)%w2 )
         enddo
         deallocate (gps_rays(n)%ip123)
      enddo
      deallocate (gps_rays)
   endif

   if (trace_use_dull) call da_trace_exit("da_gpseph_final")

end subroutine da_gpseph_final


