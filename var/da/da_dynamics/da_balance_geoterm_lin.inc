subroutine da_balance_geoterm_lin( cori, rho, u, v, term_x, term_y)

   !---------------------------------------------------------------------------
   !  Purpose: calculates linearised geostrophic term in balance equation.
   !
   !  method:  term is k x rho f u on a single level.
   !
   !  assumptions: Various (see documentation).
   !---------------------------------------------------------------------------

   implicit none
   
   real, intent(in)    :: cori(ims:ime,jms:jme)   ! Coriolis factor.
   real, intent(in)    :: rho(ims:ime,jms:jme)    ! Density
   real, intent(in)    :: u(ims:ime,jms:jme)      ! u wind increment
   real, intent(in)    :: v(ims:ime,jms:jme)      ! v wind increment
   real, intent(inout) :: term_x(ims:ime,jms:jme) ! x component of term.
   real, intent(inout) :: term_y(ims:ime,jms:jme) ! y component of term.

   integer :: is, js              ! i,j lower loop limits
   integer :: ie, je              ! i,j upper loop limits

   if (trace_use) call da_trace_entry("da_balance_geoterm_lin")

   !---------------------------------------------------------------------------
   ! [1.0] Initialise:
   !---------------------------------------------------------------------------

   ! Computation to check for edge of domain:
  ! is = its; ie = ite; js = jts; je = jte
   is = ims; ie = ime; js = jms; je = jme
   if (.not. global .and. its == ids) is = ids+1
   if (.not. global .and. ite == ide) ie = ide-1
   if (jts == jds) js = jds+1; if (jte == jde) je = jde-1

   !---------------------------------------------------------------------------
   ! [2.0] Calculate term_x = -f rho v~:
   !---------------------------------------------------------------------------

   term_x(is:ie,js:je) = term_x(is:ie,js:je) - rho(is:ie,js:je) * cori(is:ie,js:je) * v(is:ie,js:je)

   !---------------------------------------------------------------------------
   ! [3.0] Calculate term_y = f rho u~:
   !---------------------------------------------------------------------------

   term_y(is:ie,js:je) = term_y(is:ie,js:je) + rho(is:ie,js:je) * cori(is:ie,js:je) * u(is:ie,js:je)

   if (trace_use) call da_trace_exit("da_balance_geoterm_lin")

end subroutine da_balance_geoterm_lin


