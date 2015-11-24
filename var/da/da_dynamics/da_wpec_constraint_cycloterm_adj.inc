subroutine da_wpec_constraint_cycloterm_adj (rho, ub, vb, u, v, coefx, coefy, term_x, term_y)

   !---------------------------------------------------------------------------
   ! Purpose: Adjoint of da_wpec_constraint_cycloterm
   !---------------------------------------------------------------------------

   implicit none
   
   real, intent(in)    :: rho(ims:ime,jms:jme)    ! Density.
   real, intent(in)    :: ub(ims:ime,jms:jme)     ! Background u wind
   real, intent(in)    :: vb(ims:ime,jms:jme)     ! Background u wind
   real, intent(in)    :: term_x(ims:ime,jms:jme) ! x component of term
   real, intent(in)    :: term_y(ims:ime,jms:jme) ! y component of term
   real, intent(in)    :: coefx(ims:ime,jms:jme)
   real, intent(in)    :: coefy(ims:ime,jms:jme)  ! Mulplicative coeff. 

   real, intent(inout) :: u(ims:ime,jms:jme)      ! u wind increment
   real, intent(inout) :: v(ims:ime,jms:jme)      ! v wind increment

   integer              :: i, j                         ! Loop counters.
   integer              :: is, ie                       ! 1st dim. end points.
   integer              :: js, je                       ! 2nd dim. end points.
   real                 :: data(ims:ime,jms:jme)        ! Work array.

   real                 :: var, uar

   if (trace_use) call da_trace_entry("da_wpec_constraint_cycloterm_adj")

   !---------------------------------------------------------------------------
   ! [1.0] Initialise:
   !---------------------------------------------------------------------------

   is = ims; ie = ime; js = jms; je = jme
   if (.not. global .and. its == ids) is = ids+1
   if (.not. global .and. ite == ide) ie = ide-1
   if (jts == jds) js = jds+1
   if (jte == jde) je = jde-1

   !---------------------------------------------------------------------------
   ! [3.0] Calculate term_y = rho M ( u'dv/dx + v'dv/dy + udv'/dx + vdv'/dy ):
   !---------------------------------------------------------------------------

   ! [3.7] Multiply by rho and add to term_y

   data(is:ie,js:je) = rho(is:ie,js:je) * term_y(is:ie,js:je)
   
   !  [3.1] Interior points:

   do j = je, js, -1
      do i = ie, is, -1
         uar = coefx(i,j) * ub(i,j) * data(i,j)  
         var = coefy(i,j) * vb(i,j) * data(i,j)  

         u(i,j) = u(i,j) + coefx(i,j)*data(i,j)*( vb(i+1,j) - vb(i-1,j) ) 
         v(i,j) = v(i,j) + coefy(i,j)*data(i,j)*( vb(i,j+1) - vb(i,j-1) ) 
         v(i+1,j) = v(i+1,j) + uar
         v(i-1,j) = v(i-1,j) - uar
         v(i,j+1) = v(i,j+1) + var
         v(i,j-1) = v(i,j-1) - var
      end do
   end do
   
   !---------------------------------------------------------------------------
   ! [2.0] Calculate term_x = rho M ( u'du/dx + v'du/dy + udu'/dx + vdu'/dy ):
   !---------------------------------------------------------------------------

   ! [2.7] Multiply by rho and add to term_x:

   data(is:ie,js:je) = rho(is:ie,js:je) * term_x(is:ie,js:je)

   ! [2.1] Interior points:

   do j = je, js, -1
      do i = ie, is, -1
         uar = coefx(i,j) * ub(i,j) * data(i,j)
         var = coefy(i,j) * vb(i,j) * data(i,j)

         u(i,j) = u(i,j) + coefx(i,j)*( ub(i+1,j) - ub(i-1,j) ) * data(i,j)
         v(i,j) = v(i,j) + coefy(i,j)*( ub(i,j+1) - ub(i,j-1) ) * data(i,j)
         u(i+1,j) = u(i+1,j) + uar
         u(i-1,j) = u(i-1,j) - uar
         u(i,j+1) = u(i,j+1) + var
         u(i,j-1) = u(i,j-1) - var  
      end do
   end do

   if (trace_use) call da_trace_exit("da_wpec_constraint_cycloterm_adj")

end subroutine da_wpec_constraint_cycloterm_adj

