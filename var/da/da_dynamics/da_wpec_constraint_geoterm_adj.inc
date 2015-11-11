subroutine da_wpec_constraint_geoterm_adj( cori, rho, term_x, term_y, u, v)
 
   !---------------------------------------------------------------------------
   ! Purpose: Adjoint of da_wpec_constraint_geoterm.
   !---------------------------------------------------------------------------

   implicit none
   
   real, intent(in)    :: cori(ims:ime,jms:jme)   ! Coriolis factor.
   real, intent(in)    :: rho(ims:ime,jms:jme)    ! Density
   real, intent(in)    :: term_x(ims:ime,jms:jme) ! x component of term.
   real, intent(in)    :: term_y(ims:ime,jms:jme) ! y component of term.
   real, intent(inout) :: u(ims:ime,jms:jme)      ! u wind increment
   real, intent(inout) :: v(ims:ime,jms:jme)      ! v wind increment
   integer :: is, ie                       ! 1st dim. end points.
   integer :: js, je                       ! 2nd dim. end points.
   if (trace_use) call da_trace_entry("da_wpec_constraint_geoterm_adj")

   !---------------------------------------------------------------------------
   ! [1.0] Initialise:
   !---------------------------------------------------------------------------

   ! Computation to check for edge of domain:
!   is = its; ie = ite; js = jts; je = jte
   is = ims; ie = ime; js = jms; je = jme
   if (.not. global .and. its == ids) is = ids+1
   if (.not. global .and. ite == ide) ie = ide-1
   if (jts == jds) js = jds+1; if (jte == jde) je = jde-1
   !---------------------------------------------------------------------------
   ! [2.0] Calculate term_y = f rho u~:
   !---------------------------------------------------------------------------

   u(is:ie,js:je) = u(is:ie,js:je) + rho(is:ie,js:je) * cori(is:ie,js:je) &
      * term_y(is:ie,js:je)

   !---------------------------------------------------------------------------
   ! [1.0] Calculate term_x = -f rho v~:
   !---------------------------------------------------------------------------

   v(is:ie,js:je) = v(is:ie,js:je) - rho(is:ie,js:je) * cori(is:ie,js:je) &
      * term_x(is:ie,js:je)

   if (trace_use) call da_trace_exit("da_wpec_constraint_geoterm_adj")

end subroutine da_wpec_constraint_geoterm_adj


