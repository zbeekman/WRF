subroutine da_wpec_constraint_cycloterm_lin( rho, ub, vb, u, v, coefx, coefy, term_x, term_y)
                                     
   !---------------------------------------------------------------------------
   ! Purpose: Calculates linearised cyclostrophic term in wpec constraint equation.
   !
   ! Method:  Term is rho (u.grad) u on a single level.
   !---------------------------------------------------------------------------

   implicit none
   
   real, intent(in)    :: rho(ims:ime,jms:jme)    ! Density.
   real, intent(in)    :: ub(ims:ime,jms:jme)     ! Background u wind
   real, intent(in)    :: vb(ims:ime,jms:jme)     ! Background u wind
   real, intent(in)    :: u(ims:ime,jms:jme)      ! u wind increment
   real, intent(in)    :: v(ims:ime,jms:jme)      ! v wind increment
   real, intent(in)    :: coefx(ims:ime,jms:jme)  ! Multiplicative const.
   real, intent(in)    :: coefy(ims:ime,jms:jme)
   real, intent(inout) :: term_x(ims:ime,jms:jme) ! x component of term.
   real, intent(inout) :: term_y(ims:ime,jms:jme) ! y component of term.

   integer :: i, j                         ! Loop counters.
   integer :: is, ie                       ! 1st dim. end points.
   integer :: js, je                       ! 2nd dim. end points.
   
   real    :: data(ims:ime,jms:jme)        ! Temporary storage.

   real    :: varb, var

   if (trace_use) call da_trace_entry("da_wpec_constraint_cycloterm_lin")

   !---------------------------------------------------------------------------
   ! [1.0] Initialise:
   !---------------------------------------------------------------------------

   is = ims; ie = ime; js = jms; je = jme
   if (.not. global .and. its == ids) is = ids+1
   if (.not. global .and. ite == ide) ie = ide-1
   if (jts == jds) js = jds+1; if (jte == jde) je = jde-1

   !---------------------------------------------------------------------------
   ! [2.0] Calculate term_x = rho M ( u'du/dx + v'du/dy + udu'/dx + vdu'/dy):
   !---------------------------------------------------------------------------

   ! [2.1] Interior points:

   do j = js, je
      do i = is, ie
         data(i,j) = coefx(i,j)*u(i,j) * ( ub(i+1,j) - ub(i-1,j)) + &
                     coefy(i,j)*v(i,j) * ( ub(i,j+1) - ub(i,j-1)) + &
                     coefx(i,j)*ub(i,j) * ( u(i+1,j) - u(i-1,j)) + &
                     coefy(i,j)*vb(i,j) * ( u(i,j+1) - u(i,j-1))
      end do
   end do

      
   ! [2.7] Multiply by rho and add to term_x:

   term_x(is:ie,js:je) = rho(is:ie,js:je) * data(is:ie,js:je) + &
                             term_x(is:ie,js:je)

   !---------------------------------------------------------------------------
   ! [3.0] Calculate term_y = rho M ( u'dv/dx + v'dv/dy + udv'/dx + vdv'/dy):
   !---------------------------------------------------------------------------

   ! [3.1] Interior points:

   do j = js, je
      do i = is, ie
         data(i,j) = coefx(i,j)*u(i,j) * ( vb(i+1,j) - vb(i-1,j)) + &
                     coefy(i,j)*v(i,j) * ( vb(i,j+1) - vb(i,j-1)) + &
                     coefx(i,j)*ub(i,j) * ( v(i+1,j) - v(i-1,j)) + &
                     coefy(i,j)*vb(i,j) * ( v(i,j+1) - v(i,j-1))
      end do
   end do
   

   ! [3.7] Multiply by  rho and add to term_y

   term_y(is:ie,js:je) = rho(is:ie,js:je) * data(is:ie,js:je) + term_y(is:ie,js:je)

   if (trace_use) call da_trace_exit("da_wpec_constraint_cycloterm_lin")

end subroutine da_wpec_constraint_cycloterm_lin


