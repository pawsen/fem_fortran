module link1

   ! This module contains subroutines specific to the link1 element,
   ! and the definition of the link1-datatype.

   implicit none
 
   private
   public :: link1_ke, link1_ss

contains

   subroutine link1_ke(xe, young, area, ke)

   ! This subroutine constructs the stiffness matrix for
   ! a truss element.

   implicit none

   integer :: i, j
   real(8), intent(in) :: young, area
   real(8), dimension(:), intent(in) :: xe
   real(8), dimension(:,:), intent(out) :: ke

   real(8) :: delx, dely, l0, delx2, dely2, delxy, coef

   delx = xe(3) - xe(1)
   dely = xe(4) - xe(2)
   delx2 = delx * delx
   dely2 = dely * dely
   delxy = delx * dely
   l0 = sqrt(delx2 + dely2)

   coef = young*area/l0**3

   ke(1, 1) = coef * delx2
   ke(1, 2) = coef * delxy
   ke(1, 3) = -coef * delx2
   ke(1, 4) = -coef * delxy
   ke(2, 1) = coef * delxy
   ke(2, 2) = coef * dely2
   ke(2, 3) = -coef * delxy
   ke(2, 4) = -coef * dely2
   ke(3, 1) = -coef * delx2
   ke(3, 2) = -coef * delxy
   ke(3, 3) = coef * delx2
   ke(3, 4) = coef * delxy
   ke(4, 1) = -coef * delxy
   ke(4, 2) = -coef * dely2
   ke(4, 3) = coef * delxy
   ke(4, 4) = coef * dely2

end subroutine link1_ke

subroutine link1_ss(xe, de, young, estress, estrain)

   ! This subrotuine constructs the element stress and strain

   implicit none

   integer :: i, j
   real(8), intent(in) :: young
   real(8), dimension(:), intent(in) :: xe, de
   real(8), dimension(:), intent(out) ::  estress, estrain

   real(8) :: delx, dely, l0, delx2, dely2, delu, delv

   delx = xe(3) - xe(1)
   dely = xe(4) - xe(2)
   delx2 = delx * delx
   dely2 = dely * dely
   l0 = sqrt(delx2 + dely2)

   delu = de(3) - de(1)
   delv = de(4) - de(2)

   estrain(1) = (delu*delx + delv*dely)/l0**2
   estrain(2:3) = 0.
   estress(1) = young * estrain(1)
   estress(2:3) = 0.

end subroutine link1_ss

end module link1
