subroutine da_gpseph_init (grid)

   ! -------------------------------------------------------------------
   ! Purpose: Gather the lat, lon, terr, h
   ! -------------------------------------------------------------------

   implicit none

   type (domain), intent(in)      :: grid       ! first guess state

   integer :: i, j, k, ij, ijk

   if (trace_use_dull) call da_trace_entry("da_gpseph_init")

   ! Gather lat, lon, h and terr for gpseph data
   if (gpseph_nonlocal)   then

      allocate ( global_lat (ids:ide,jds:jde) )
      allocate ( global_lon (ids:ide,jds:jde) )
      allocate ( global_terr (ids:ide,jds:jde) )
      allocate ( global_h (ids:ide,jds:jde,kds:kde) )
      allocate ( global_ref (ids:ide,jds:jde,kds:kde) )
      allocate ( global_xa_ref (ids:ide,jds:jde,kds:kde) )
      allocate ( global_adj_ref (ids:ide,jds:jde,kds:kde) )
      global_lat = 0.0
      global_lon = 0.0
      global_terr = 0.0
      global_h = 0.0
      global_ref = 0.0
      global_xa_ref = 0.0
      global_adj_ref = 0.0

#ifdef DM_PARALLEL
      ij = ( ide- ids+1)*( jde- jds+1)
      ijk = ( ide- ids+1)*( jde- jds+1)*(kde- kds+1)
   !  Collect xb component of h into global buffer.
      call da_patch_to_global( grid, grid%xb%h, global_h )
      call wrf_dm_bcast_real( global_h, ijk )
   !  Collect xb component of ref into global buffer.
      call da_patch_to_global( grid, grid%xb%ref, global_ref )
      call wrf_dm_bcast_real( global_ref, ijk )
   !  Collect xb component of lat into global buffer.
      call da_patch_to_global( grid, grid%xb%lat, global_lat )
      call wrf_dm_bcast_real( global_lat, ij )
   !  Collect xb component of lon into global buffer.
      call da_patch_to_global( grid, grid%xb%lon, global_lon )
      call wrf_dm_bcast_real( global_lon, ij )
   !  Collect xb component of terr into global buffer.
      call da_patch_to_global( grid, grid%xb%terr, global_terr )
      call wrf_dm_bcast_real( global_terr, ij )
#else
      do k = kds, kde
         do j = jds, jde
            do i = ids, ide
               global_h(i,j,k) = grid%xb%h(i,j,k)
               global_ref(i,j,k) = grid%xb%ref(i,j,k)
            enddo
         enddo
      enddo
      do j = jds, jde
         do i = ids, ide
            global_lat(i,j) = grid%xb%lat(i,j)
            global_lon(i,j) = grid%xb%lon(i,j)
            global_terr(i,j) = grid%xb%terr(i,j)
         enddo
      enddo
#endif

   endif

   if (trace_use_dull) call da_trace_exit("da_gpseph_init")

end subroutine da_gpseph_init


