MODULE plate41rect 

 ! This module contains subroutines specific to the plate41 element.

 IMPLICIT NONE

 PRIVATE
 PUBLIC :: plate41rect_ke,plate41rect_me, plate41rect_re, plate41rect_ss

CONTAINS

 SUBROUTINE plate41rect_ke(xe, young, nu,thk, ke)

  ! This subroutine constructs the stiffness matrix for
  ! a rectangular 4-noded quad element.

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: young, nu, thk
  REAL(8), DIMENSION(:), INTENT(IN) :: xe
  REAL(8), DIMENSION(:,:), INTENT(OUT) :: ke

  REAL(8) :: aa, bb
  
  aa = (xe(4)-xe(1))/2.
  bb = (xe(11)-xe(2))/2.
  
  ! Build element stiffness matrix:
      ke(1,1) = -young*thk**3/aa**3/bb**3*(7*aa**2*bb**2-2*bb**2*aa**2*nu+10*bb**4+10*aa**4)/(-1+nu**2)/120
      ke(1,2) = -young*thk**3/aa**2/bb*(10*bb**2+4*aa**2*nu+aa**2)/(-1+nu**2)/120
      ke(1,3) = -young*thk**3/aa/bb**2*(bb**2+4*nu*bb**2+10*aa**2)/(-1+nu**2)/120
      ke(1,4) = -young*thk**3/aa**3/bb**3*(-7*aa**2*bb**2+2*bb**2*aa**2*nu-10*bb**4+5*aa**4)/(-1+nu**2)/120
      ke(1,5) = young*thk**3/aa**2/bb*(-10*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/120
      ke(1,6) = -young*thk**3/aa/bb**2*(-bb**2-4*nu*bb**2+5*aa**2)/(-1+nu**2)/120
      ke(1,7) = young*thk**3/aa**3/bb**3*(-7*aa**2*bb**2+2*bb**2*aa**2*nu+5*bb**4+5*aa**4)/(-1+nu**2)/120
      ke(1,8) = -young*thk**3/aa**2/bb*(5*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/120
      ke(1,9) = -young*thk**3/aa/bb**2*(-bb**2+nu*bb**2+5*aa**2)/(-1+nu**2)/120
      ke(1,10) = young*thk**3/aa**3/bb**3*(7*aa**2*bb**2-2*bb**2*aa**2*nu-5*bb**4+10*aa**4)/(-1+nu**2)/120
      ke(1,11) = young*thk**3/aa**2/bb*(-5*bb**2+4*aa**2*nu+aa**2)/(-1+nu**2)/120
      ke(1,12) = -young*thk**3/aa/bb**2*(bb**2-nu*bb**2+10*aa**2)/(-1+nu**2)/120
      ke(2,1) = -young*thk**3/aa**2/bb*(10*bb**2+4*aa**2*nu+aa**2)/(-1+nu**2)/120
      ke(2,2) = young*thk**3/aa/bb*(-5*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/45
      ke(2,3) = -nu*young*thk**3/(-1+nu**2)/12
      ke(2,4) = -young*thk**3/aa**2/bb*(-10*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/120
      ke(2,5) = -young*thk**3/aa/bb*(10*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/180
      ke(2,6) = 0
      ke(2,7) = young*thk**3/aa**2/bb*(5*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/120
      ke(2,8) = young*thk**3/aa/bb*(-5*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/180
      ke(2,9) = 0
      ke(2,10) = young*thk**3/aa**2/bb*(-5*bb**2+4*aa**2*nu+aa**2)/(-1+nu**2)/120
      ke(2,11) = -young*thk**3/aa/bb*(5*bb**2-2*aa**2+2*aa**2*nu)/(-1+nu**2)/90
      ke(2,12) = 0
      ke(3,1) = -young*thk**3/aa/bb**2*(bb**2+4*nu*bb**2+10*aa**2)/(-1+nu**2)/120
      ke(3,2) = -nu*young*thk**3/(-1+nu**2)/12
      ke(3,3) = -1/bb*young*thk**3*(bb**2-nu*bb**2+5*aa**2)/(-1+nu**2)/aa/45
      ke(3,4) = -young*thk**3/aa/bb**2*(-bb**2-4*nu*bb**2+5*aa**2)/(-1+nu**2)/120
      ke(3,5) = 0
      ke(3,6) = -1/bb*young*thk**3*(-2*bb**2+2*nu*bb**2+5*aa**2)/(-1+nu**2)/aa/90
      ke(3,7) = young*thk**3/aa/bb**2*(-bb**2+nu*bb**2+5*aa**2)/(-1+nu**2)/120
      ke(3,8) = 0
      ke(3,9) = -1/bb*young*thk**3*(bb**2-nu*bb**2+5*aa**2)/(-1+nu**2)/aa/180
      ke(3,10) = young*thk**3/aa/bb**2*(bb**2-nu*bb**2+10*aa**2)/(-1+nu**2)/120
      ke(3,11) = 0
      ke(3,12) = -1/bb*young*thk**3*(-bb**2+nu*bb**2+10*aa**2)/(-1+nu**2)/aa/180
      ke(4,1) = -young*thk**3/aa**3/bb**3*(-7*aa**2*bb**2+2*bb**2*aa**2*nu-10*bb**4+5*aa**4)/(-1+nu**2)/120
      ke(4,2) = -young*thk**3/aa**2/bb*(-10*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/120
      ke(4,3) = -young*thk**3/aa/bb**2*(-bb**2-4*nu*bb**2+5*aa**2)/(-1+nu**2)/120
      ke(4,4) = -young*thk**3/aa**3/bb**3*(7*aa**2*bb**2-2*bb**2*aa**2*nu+10*bb**4+10*aa**4)/(-1+nu**2)/120
      ke(4,5) = young*thk**3/aa**2/bb*(10*bb**2+4*aa**2*nu+aa**2)/(-1+nu**2)/120
      ke(4,6) = -young*thk**3/aa/bb**2*(bb**2+4*nu*bb**2+10*aa**2)/(-1+nu**2)/120
      ke(4,7) = young*thk**3/aa**3/bb**3*(7*aa**2*bb**2-2*bb**2*aa**2*nu-5*bb**4+10*aa**4)/(-1+nu**2)/120
      ke(4,8) = -young*thk**3/aa**2/bb*(-5*bb**2+4*aa**2*nu+aa**2)/(-1+nu**2)/120
      ke(4,9) = -young*thk**3/aa/bb**2*(bb**2-nu*bb**2+10*aa**2)/(-1+nu**2)/120
      ke(4,10) = young*thk**3/aa**3/bb**3*(-7*aa**2*bb**2+2*bb**2*aa**2*nu+5*bb**4+5*aa**4)/(-1+nu**2)/120
      ke(4,11) = young*thk**3/aa**2/bb*(5*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/120
      ke(4,12) = -young*thk**3/aa/bb**2*(-bb**2+nu*bb**2+5*aa**2)/(-1+nu**2)/120
      ke(5,1) = young*thk**3/aa**2/bb*(-10*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/120
      ke(5,2) = -young*thk**3/aa/bb*(10*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/180
      ke(5,3) = 0
      ke(5,4) = young*thk**3/aa**2/bb*(10*bb**2+4*aa**2*nu+aa**2)/(-1+nu**2)/120
      ke(5,5) = young*thk**3/aa/bb*(-5*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/45
      ke(5,6) = nu*young*thk**3/(-1+nu**2)/12
      ke(5,7) = -young*thk**3/aa**2/bb*(-5*bb**2+4*aa**2*nu+aa**2)/(-1+nu**2)/120
      ke(5,8) = -young*thk**3/aa/bb*(5*bb**2-2*aa**2+2*aa**2*nu)/(-1+nu**2)/90
      ke(5,9) = 0
      ke(5,10) = -young*thk**3/aa**2/bb*(5*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/120
      ke(5,11) = young*thk**3/aa/bb*(-5*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/180
      ke(5,12) = 0
      ke(6,1) = -young*thk**3/aa/bb**2*(-bb**2-4*nu*bb**2+5*aa**2)/(-1+nu**2)/120
      ke(6,2) = 0
      ke(6,3) = -1/bb*young*thk**3*(-2*bb**2+2*nu*bb**2+5*aa**2)/(-1+nu**2)/aa/90
      ke(6,4) = -young*thk**3/aa/bb**2*(bb**2+4*nu*bb**2+10*aa**2)/(-1+nu**2)/120
      ke(6,5) = nu*young*thk**3/(-1+nu**2)/12
      ke(6,6) = -1/bb*young*thk**3*(bb**2-nu*bb**2+5*aa**2)/(-1+nu**2)/aa/45
      ke(6,7) = young*thk**3/aa/bb**2*(bb**2-nu*bb**2+10*aa**2)/(-1+nu**2)/120
      ke(6,8) = 0
      ke(6,9) = -1/bb*young*thk**3*(-bb**2+nu*bb**2+10*aa**2)/(-1+nu**2)/aa/180
      ke(6,10) = young*thk**3/aa/bb**2*(-bb**2+nu*bb**2+5*aa**2)/(-1+nu**2)/120
      ke(6,11) = 0
      ke(6,12) = -1/bb*young*thk**3*(bb**2-nu*bb**2+5*aa**2)/(-1+nu**2)/aa/180
      ke(7,1) = young*thk**3/aa**3/bb**3*(-7*aa**2*bb**2+2*bb**2*aa**2*nu+5*bb**4+5*aa**4)/(-1+nu**2)/120
      ke(7,2) = young*thk**3/aa**2/bb*(5*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/120
      ke(7,3) = young*thk**3/aa/bb**2*(-bb**2+nu*bb**2+5*aa**2)/(-1+nu**2)/120
      ke(7,4) = young*thk**3/aa**3/bb**3*(7*aa**2*bb**2-2*bb**2*aa**2*nu-5*bb**4+10*aa**4)/(-1+nu**2)/120
      ke(7,5) = -young*thk**3/aa**2/bb*(-5*bb**2+4*aa**2*nu+aa**2)/(-1+nu**2)/120
      ke(7,6) = young*thk**3/aa/bb**2*(bb**2-nu*bb**2+10*aa**2)/(-1+nu**2)/120
      ke(7,7) = -young*thk**3/aa**3/bb**3*(7*aa**2*bb**2-2*bb**2*aa**2*nu+10*bb**4+10*aa**4)/(-1+nu**2)/120
      ke(7,8) = young*thk**3/aa**2/bb*(10*bb**2+4*aa**2*nu+aa**2)/(-1+nu**2)/120
      ke(7,9) = young*thk**3/aa/bb**2*(bb**2+4*nu*bb**2+10*aa**2)/(-1+nu**2)/120
      ke(7,10) = -young*thk**3/aa**3/bb**3*(-7*aa**2*bb**2+2*bb**2*aa**2*nu-10*bb**4+5*aa**4)/(-1+nu**2)/120
      ke(7,11) = -young*thk**3/aa**2/bb*(-10*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/120
      ke(7,12) = young*thk**3/aa/bb**2*(-bb**2-4*nu*bb**2+5*aa**2)/(-1+nu**2)/120
      ke(8,1) = -young*thk**3/aa**2/bb*(5*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/120
      ke(8,2) = young*thk**3/aa/bb*(-5*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/180
      ke(8,3) = 0
      ke(8,4) = -young*thk**3/aa**2/bb*(-5*bb**2+4*aa**2*nu+aa**2)/(-1+nu**2)/120
      ke(8,5) = -young*thk**3/aa/bb*(5*bb**2-2*aa**2+2*aa**2*nu)/(-1+nu**2)/90
      ke(8,6) = 0
      ke(8,7) = young*thk**3/aa**2/bb*(10*bb**2+4*aa**2*nu+aa**2)/(-1+nu**2)/120
      ke(8,8) = young*thk**3/aa/bb*(-5*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/45
      ke(8,9) = -nu*young*thk**3/(-1+nu**2)/12
      ke(8,10) = young*thk**3/aa**2/bb*(-10*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/120
      ke(8,11) = -young*thk**3/aa/bb*(10*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/180
      ke(8,12) = 0
      ke(9,1) = -young*thk**3/aa/bb**2*(-bb**2+nu*bb**2+5*aa**2)/(-1+nu**2)/120
      ke(9,2) = 0
      ke(9,3) = -1/bb*young*thk**3*(bb**2-nu*bb**2+5*aa**2)/(-1+nu**2)/aa/180
      ke(9,4) = -young*thk**3/aa/bb**2*(bb**2-nu*bb**2+10*aa**2)/(-1+nu**2)/120
      ke(9,5) = 0
      ke(9,6) = -1/bb*young*thk**3*(-bb**2+nu*bb**2+10*aa**2)/(-1+nu**2)/aa/180
      ke(9,7) = young*thk**3/aa/bb**2*(bb**2+4*nu*bb**2+10*aa**2)/(-1+nu**2)/120
      ke(9,8) = -nu*young*thk**3/(-1+nu**2)/12
      ke(9,9) = -1/bb*young*thk**3*(bb**2-nu*bb**2+5*aa**2)/(-1+nu**2)/aa/45
      ke(9,10) = young*thk**3/aa/bb**2*(-bb**2-4*nu*bb**2+5*aa**2)/(-1+nu**2)/120
      ke(9,11) = 0
      ke(9,12) = -1/bb*young*thk**3*(-2*bb**2+2*nu*bb**2+5*aa**2)/(-1+nu**2)/aa/90
      ke(10,1) = young*thk**3/aa**3/bb**3*(7*aa**2*bb**2-2*bb**2*aa**2*nu-5*bb**4+10*aa**4)/(-1+nu**2)/120
      ke(10,2) = young*thk**3/aa**2/bb*(-5*bb**2+4*aa**2*nu+aa**2)/(-1+nu**2)/120
      ke(10,3) = young*thk**3/aa/bb**2*(bb**2-nu*bb**2+10*aa**2)/(-1+nu**2)/120
      ke(10,4) = young*thk**3/aa**3/bb**3*(-7*aa**2*bb**2+2*bb**2*aa**2*nu+5*bb**4+5*aa**4)/(-1+nu**2)/120
      ke(10,5) = -young*thk**3/aa**2/bb*(5*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/120
      ke(10,6) = young*thk**3/aa/bb**2*(-bb**2+nu*bb**2+5*aa**2)/(-1+nu**2)/120
      ke(10,7) = -young*thk**3/aa**3/bb**3*(-7*aa**2*bb**2+2*bb**2*aa**2*nu-10*bb**4+5*aa**4)/(-1+nu**2)/120
      ke(10,8) = young*thk**3/aa**2/bb*(-10*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/120
      ke(10,9) = young*thk**3/aa/bb**2*(-bb**2-4*nu*bb**2+5*aa**2)/(-1+nu**2)/120
      ke(10,10) = -young*thk**3/aa**3/bb**3*(7*aa**2*bb**2-2*bb**2*aa**2*nu+10*bb**4+10*aa**4)/(-1+nu**2)/120
      ke(10,11) = -young*thk**3/aa**2/bb*(10*bb**2+4*aa**2*nu+aa**2)/(-1+nu**2)/120
      ke(10,12) = young*thk**3/aa/bb**2*(bb**2+4*nu*bb**2+10*aa**2)/(-1+nu**2)/120
      ke(11,1) = young*thk**3/aa**2/bb*(-5*bb**2+4*aa**2*nu+aa**2)/(-1+nu**2)/120
      ke(11,2) = -young*thk**3/aa/bb*(5*bb**2-2*aa**2+2*aa**2*nu)/(-1+nu**2)/90
      ke(11,3) = 0
      ke(11,4) = young*thk**3/aa**2/bb*(5*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/120
      ke(11,5) = young*thk**3/aa/bb*(-5*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/180
      ke(11,6) = 0
      ke(11,7) = -young*thk**3/aa**2/bb*(-10*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/120
      ke(11,8) = -young*thk**3/aa/bb*(10*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/180
      ke(11,9) = 0
      ke(11,10) = -young*thk**3/aa**2/bb*(10*bb**2+4*aa**2*nu+aa**2)/(-1+nu**2)/120
      ke(11,11) = young*thk**3/aa/bb*(-5*bb**2-aa**2+aa**2*nu)/(-1+nu**2)/45
      ke(11,12) = nu*young*thk**3/(-1+nu**2)/12
      ke(12,1) = -young*thk**3/aa/bb**2*(bb**2-nu*bb**2+10*aa**2)/(-1+nu**2)/120
      ke(12,2) = 0
      ke(12,3) = -1/bb*young*thk**3*(-bb**2+nu*bb**2+10*aa**2)/(-1+nu**2)/aa/180
      ke(12,4) = -young*thk**3/aa/bb**2*(-bb**2+nu*bb**2+5*aa**2)/(-1+nu**2)/120
      ke(12,5) = 0
      ke(12,6) = -1/bb*young*thk**3*(bb**2-nu*bb**2+5*aa**2)/(-1+nu**2)/aa/180
      ke(12,7) = young*thk**3/aa/bb**2*(-bb**2-4*nu*bb**2+5*aa**2)/(-1+nu**2)/120
      ke(12,8) = 0
      ke(12,9) = -1/bb*young*thk**3*(-2*bb**2+2*nu*bb**2+5*aa**2)/(-1+nu**2)/aa/90
      ke(12,10) = young*thk**3/aa/bb**2*(bb**2+4*nu*bb**2+10*aa**2)/(-1+nu**2)/120
      ke(12,11) = nu*young*thk**3/(-1+nu**2)/12
      ke(12,12) = -1/bb*young*thk**3*(bb**2-nu*bb**2+5*aa**2)/(-1+nu**2)/aa/45
      
 end SUBROUTINE plate41rect_ke

SUBROUTINE plate41rect_me(xe, young, nu, dens, thk, me)
  ! This subroutine constructs the stiffness matrix for
  ! a rectangular 4-noded quad element.

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: young, nu, thk, dens
  REAL(8), DIMENSION(:), INTENT(IN) :: xe
  REAL(8), DIMENSION(:,:), INTENT(OUT) :: me

  REAL(8) :: aa, bb
  
  aa = (xe(4)-xe(1))/2.
  bb = (xe(11)-xe(2))/2.
   ! Computing the element mass matrix
     me = 0.0d0

     ! Mass matrix - Calculated
      Me(1,1) = 1727.D0/3150.D0*aa*bb
      Me(1,2) = 461.D0/3150.D0*bb*aa**2
      Me(1,3) = 461.D0/3150.D0*bb**2*aa
      Me(1,4) = 613.D0/3150.D0*aa*bb
      Me(1,5) = -137.D0/1575.D0*bb*aa**2
      Me(1,6) = 199.D0/3150.D0*bb**2*aa
      Me(1,7) = 197.D0/3150.D0*aa*bb
      Me(1,8) = -58.D0/1575.D0*bb*aa**2
      Me(1,9) = -58.D0/1575.D0*bb**2*aa
      Me(1,10) = 613.D0/3150.D0*aa*bb
      Me(1,11) = 199.D0/3150.D0*bb*aa**2
      Me(1,12) = -137.D0/1575.D0*bb**2*aa
      Me(2,1) = 461.D0/3150.D0*bb*aa**2
      Me(2,2) = 16.D0/315.D0*bb*aa**3
      Me(2,3) = bb**2*aa**2/25
      Me(2,4) = 137.D0/1575.D0*bb*aa**2
      Me(2,5) = -4.D0/105.D0*bb*aa**3
      Me(2,6) = 2.D0/75.D0*bb**2*aa**2
      Me(2,7) = 58.D0/1575.D0*bb*aa**2
      Me(2,8) = -2.D0/105.D0*bb*aa**3
      Me(2,9) = -4.D0/225.D0*bb**2*aa**2
      Me(2,10) = 199.D0/3150.D0*bb*aa**2
      Me(2,11) = 8.D0/315.D0*bb*aa**3
      Me(2,12) = -2.D0/75.D0*bb**2*aa**2
      Me(3,1) = 461.D0/3150.D0*bb**2*aa
      Me(3,2) = bb**2*aa**2/25
      Me(3,3) = 16.D0/315.D0*bb**3*aa
      Me(3,4) = 199.D0/3150.D0*bb**2*aa
      Me(3,5) = -2.D0/75.D0*bb**2*aa**2
      Me(3,6) = 8.D0/315.D0*bb**3*aa
      Me(3,7) = 58.D0/1575.D0*bb**2*aa
      Me(3,8) = -4.D0/225.D0*bb**2*aa**2
      Me(3,9) = -2.D0/105.D0*bb**3*aa
      Me(3,10) = 137.D0/1575.D0*bb**2*aa
      Me(3,11) = 2.D0/75.D0*bb**2*aa**2
      Me(3,12) = -4.D0/105.D0*bb**3*aa
      Me(4,1) = 613.D0/3150.D0*aa*bb
      Me(4,2) = 137.D0/1575.D0*bb*aa**2
      Me(4,3) = 199.D0/3150.D0*bb**2*aa
      Me(4,4) = 1727.D0/3150.D0*aa*bb
      Me(4,5) = -461.D0/3150.D0*bb*aa**2
      Me(4,6) = 461.D0/3150.D0*bb**2*aa
      Me(4,7) = 613.D0/3150.D0*aa*bb
      Me(4,8) = -199.D0/3150.D0*bb*aa**2
      Me(4,9) = -137.D0/1575.D0*bb**2*aa
      Me(4,10) = 197.D0/3150.D0*aa*bb
      Me(4,11) = 58.D0/1575.D0*bb*aa**2
      Me(4,12) = -58.D0/1575.D0*bb**2*aa
      Me(5,1) = -137.D0/1575.D0*bb*aa**2
      Me(5,2) = -4.D0/105.D0*bb*aa**3
      Me(5,3) = -2.D0/75.D0*bb**2*aa**2
      Me(5,4) = -461.D0/3150.D0*bb*aa**2
      Me(5,5) = 16.D0/315.D0*bb*aa**3
      Me(5,6) = -bb**2*aa**2/25
      Me(5,7) = -199.D0/3150.D0*bb*aa**2
      Me(5,8) = 8.D0/315.D0*bb*aa**3
      Me(5,9) = 2.D0/75.D0*bb**2*aa**2
      Me(5,10) = -58.D0/1575.D0*bb*aa**2
      Me(5,11) = -2.D0/105.D0*bb*aa**3
      Me(5,12) = 4.D0/225.D0*bb**2*aa**2
      Me(6,1) = 199.D0/3150.D0*bb**2*aa
      Me(6,2) = 2.D0/75.D0*bb**2*aa**2
      Me(6,3) = 8.D0/315.D0*bb**3*aa
      Me(6,4) = 461.D0/3150.D0*bb**2*aa
      Me(6,5) = -bb**2*aa**2/25
      Me(6,6) = 16.D0/315.D0*bb**3*aa
      Me(6,7) = 137.D0/1575.D0*bb**2*aa
      Me(6,8) = -2.D0/75.D0*bb**2*aa**2
      Me(6,9) = -4.D0/105.D0*bb**3*aa
      Me(6,10) = 58.D0/1575.D0*bb**2*aa
      Me(6,11) = 4.D0/225.D0*bb**2*aa**2
      Me(6,12) = -2.D0/105.D0*bb**3*aa
      Me(7,1) = 197.D0/3150.D0*aa*bb
      Me(7,2) = 58.D0/1575.D0*bb*aa**2
      Me(7,3) = 58.D0/1575.D0*bb**2*aa
      Me(7,4) = 613.D0/3150.D0*aa*bb
      Me(7,5) = -199.D0/3150.D0*bb*aa**2
      Me(7,6) = 137.D0/1575.D0*bb**2*aa
      Me(7,7) = 1727.D0/3150.D0*aa*bb
      Me(7,8) = -461.D0/3150.D0*bb*aa**2
      Me(7,9) = -461.D0/3150.D0*bb**2*aa
      Me(7,10) = 613.D0/3150.D0*aa*bb
      Me(7,11) = 137.D0/1575.D0*bb*aa**2
      Me(7,12) = -199.D0/3150.D0*bb**2*aa
      Me(8,1) = -58.D0/1575.D0*bb*aa**2
      Me(8,2) = -2.D0/105.D0*bb*aa**3
      Me(8,3) = -4.D0/225.D0*bb**2*aa**2
      Me(8,4) = -199.D0/3150.D0*bb*aa**2
      Me(8,5) = 8.D0/315.D0*bb*aa**3
      Me(8,6) = -2.D0/75.D0*bb**2*aa**2
      Me(8,7) = -461.D0/3150.D0*bb*aa**2
      Me(8,8) = 16.D0/315.D0*bb*aa**3
      Me(8,9) = bb**2*aa**2/25
      Me(8,10) = -137.D0/1575.D0*bb*aa**2
      Me(8,11) = -4.D0/105.D0*bb*aa**3
      Me(8,12) = 2.D0/75.D0*bb**2*aa**2
      Me(9,1) = -58.D0/1575.D0*bb**2*aa
      Me(9,2) = -4.D0/225.D0*bb**2*aa**2
      Me(9,3) = -2.D0/105.D0*bb**3*aa
      Me(9,4) = -137.D0/1575.D0*bb**2*aa
      Me(9,5) = 2.D0/75.D0*bb**2*aa**2
      Me(9,6) = -4.D0/105.D0*bb**3*aa
      Me(9,7) = -461.D0/3150.D0*bb**2*aa
      Me(9,8) = bb**2*aa**2/25
      Me(9,9) = 16.D0/315.D0*bb**3*aa
      Me(9,10) = -199.D0/3150.D0*bb**2*aa
      Me(9,11) = -2.D0/75.D0*bb**2*aa**2
      Me(9,12) = 8.D0/315.D0*bb**3*aa
      Me(10,1) = 613.D0/3150.D0*aa*bb
      Me(10,2) = 199.D0/3150.D0*bb*aa**2
      Me(10,3) = 137.D0/1575.D0*bb**2*aa
      Me(10,4) = 197.D0/3150.D0*aa*bb
      Me(10,5) = -58.D0/1575.D0*bb*aa**2
      Me(10,6) = 58.D0/1575.D0*bb**2*aa
      Me(10,7) = 613.D0/3150.D0*aa*bb
      Me(10,8) = -137.D0/1575.D0*bb*aa**2
      Me(10,9) = -199.D0/3150.D0*bb**2*aa
      Me(10,10) = 1727.D0/3150.D0*aa*bb
      Me(10,11) = 461.D0/3150.D0*bb*aa**2
      Me(10,12) = -461.D0/3150.D0*bb**2*aa
      Me(11,1) = 199.D0/3150.D0*bb*aa**2
      Me(11,2) = 8.D0/315.D0*bb*aa**3
      Me(11,3) = 2.D0/75.D0*bb**2*aa**2
      Me(11,4) = 58.D0/1575.D0*bb*aa**2
      Me(11,5) = -2.D0/105.D0*bb*aa**3
      Me(11,6) = 4.D0/225.D0*bb**2*aa**2
      Me(11,7) = 137.D0/1575.D0*bb*aa**2
      Me(11,8) = -4.D0/105.D0*bb*aa**3
      Me(11,9) = -2.D0/75.D0*bb**2*aa**2
      Me(11,10) = 461.D0/3150.D0*bb*aa**2
      Me(11,11) = 16.D0/315.D0*bb*aa**3
      Me(11,12) = -bb**2*aa**2/25
      Me(12,1) = -137.D0/1575.D0*bb**2*aa
      Me(12,2) = -2.D0/75.D0*bb**2*aa**2
      Me(12,3) = -4.D0/105.D0*bb**3*aa
      Me(12,4) = -58.D0/1575.D0*bb**2*aa
      Me(12,5) = 4.D0/225.D0*bb**2*aa**2
      Me(12,6) = -2.D0/105.D0*bb**3*aa
      Me(12,7) = -199.D0/3150.D0*bb**2*aa
      Me(12,8) = 2.D0/75.D0*bb**2*aa**2
      Me(12,9) = 8.D0/315.D0*bb**3*aa
      Me(12,10) = -461.D0/3150.D0*bb**2*aa
      Me(12,11) = -bb**2*aa**2/25
      Me(12,12) = 16.D0/315.D0*bb**3*aa

      Me = Me*thk*dens

  
 END SUBROUTINE plate41rect_me

 SUBROUTINE plate41rect_re(xe, eface, fe, thk, re)

  ! This subroutine assembles the element surface loads.

  INTEGER, INTENT(IN) :: eface
  REAL(8), INTENT(IN) :: fe, thk
  REAL(8), DIMENSION(:), INTENT(IN) :: xe
  REAL(8), INTENT(OUT) :: re(12)

  REAL(8) :: aa, bb
  
  !      
  !  l ______k
  !   |     |
  !   |  1  |  
  !   |_____|
  !  i       j
  !      
  !
  ! node numbers: i, j, k, l
  ! face numbers: 1(j->i), 2(k->j), 3(l->k), 4(i->l)

  aa = (xe(4)-xe(1))/2.
  bb = (xe(11)-xe(2))/2.
  
  !Påført moment i x-retning
  IF (eface == 1) THEN
   re(2) = fe
   re(5) = fe
 
  !Påført moment i y-retning
  ELSEIF (eface == 2) THEN
   re(6) = fe
   re(9) = fe

  ELSEIF (eface == 3) THEN
   re(8) = fe
   re(11)= fe

  ELSEIF (eface == 4) THEN
   re(12)= fe
   re(3) = fe

  !fordelt last på latteral side
  ELSEIF (eface == 5) THEN
   !Forces in z-dir.
   re(1) = aa*bb*fe !Node 1
   re(4) = aa*bb*fe !Node 2
   re(7) = aa*bb*fe !Node 3
   re(10)= aa*bb*fe !Node 4
      
   !Moments x-dir.
   re(2) = aa*bb*fe*bb !Node 1
   re(5) = aa*bb*fe*bb !Node 2
   re(8) = -aa*bb*fe*bb !Node 3
   re(11)= -aa*bb*fe*bb !Node 4

   !Moments y-dir.
   re(3) = -aa*bb*fe*aa !Node 1
   re(6) = aa*bb*fe*aa !Node 2
   re(9) = aa*bb*fe*aa !Node 3
   re(12)= -aa*bb*fe*aa !Node 4

   ELSE
     write (*, *) 'Error: Wrong face number!!'
     stop

  ENDIF

  !!!!
  !!!!! GANG MED TYKKELSE!!!
  !!!!!
  re = re*thk

 END SUBROUTINE plate41rect_re

 SUBROUTINE plate41rect_ss(xe, de, young, nu, thk, estress, estrain)

  ! This subrotuine constructs the element stress and strain
 
  REAL(8), INTENT(IN) :: young, nu, thk 
  REAL(8), DIMENSION(:), INTENT(IN)  :: xe, de
!  REAL(8), INTENT(OUT) :: sigmavm, psi, sigma1, sigma2
  REAL(8), DIMENSION(:), INTENT(OUT) :: estress, estrain

  REAL(8) :: Bmat(3, 12), Dmat(3,3), Emat(3,3)
  REAL(8) :: aa, bb, x, y, z, c, s
  
  aa = (xe(4)-xe(1))/2.
  bb = (xe(11)-xe(2))/2.
  x = 0
  y = 0
  z = thk/2. ! distance from mid-plane.

  Bmat=0.

  ! Build strain-displacement matrix
  bmat(1,1) = 3.D0/4.D0*x/aa**3-3.D0/4.D0*x*y/bb/aa**3
  bmat(1,2) = -1/aa/4+3.D0/4.D0*x/aa**2+y/bb/aa/4-3.D0/4.D0*x*y/bb/aa**2
  bmat(1,3) = 0.
  bmat(1,4) = -3.D0/4.D0*x/aa**3+3.D0/4.D0*x*y/bb/aa**3
  bmat(1,5) = 1/aa/4+3.D0/4.D0*x/aa**2-y/bb/aa/4-3.D0/4.D0*x*y/bb/aa**2
  bmat(1,6) = 0.
  bmat(1,7) = -3.D0/4.D0*x/aa**3-3.D0/4.D0*x*y/bb/aa**3
  bmat(1,8) = 1/aa/4+3.D0/4.D0*x/aa**2+y/bb/aa/4+3.D0/4.D0*x*y/bb/aa**2
  bmat(1,9) = 0.
  bmat(1,10) = 3.D0/4.D0*x/aa**3+3.D0/4.D0*x*y/bb/aa**3
  bmat(1,11) = -1/aa/4+3.D0/4.D0*x/aa**2-y/bb/aa/4+3.D0/4.D0*x*y/bb/aa**2
  bmat(1,12) = 0.
  bmat(2,1) = 3.D0/4.D0*y/bb**3-3.D0/4.D0*y*x/bb**3/aa
  bmat(2,2) = 0.
  bmat(2,3) = -1/bb/4+x/bb/aa/4+3.D0/4.D0*y/bb**2-3.D0/4.D0*y*x/bb**2/aa
  bmat(2,4) = 3.D0/4.D0*y/bb**3+3.D0/4.D0*y*x/bb**3/aa
  bmat(2,5) = 0.
  bmat(2,6) = -1/bb/4-x/bb/aa/4+3.D0/4.D0*y/bb**2+3.D0/4.D0*y*x/bb**2/aa
  bmat(2,7) = -3.D0/4.D0*y/bb**3-3.D0/4.D0*y*x/bb**3/aa
  bmat(2,8) = 0.
  bmat(2,9) = 1/bb/4+x/bb/aa/4+3.D0/4.D0*y/bb**2+3.D0/4.D0*y*x/bb**2/aa
  bmat(2,10) = -3.D0/4.D0*y/bb**3+3.D0/4.D0*y*x/bb**3/aa
  bmat(2,11) = 0.
  bmat(2,12) = 1/bb/4-x/bb/aa/4+3.D0/4.D0*y/bb**2-3.D0/4.D0*y*x/bb**2/aa
  bmat(3,1) = 1/bb/aa-3.D0/4.D0*x**2/bb/aa**3-3.D0/4.D0*y**2/bb**3/aa
  bmat(3,2) = 1/bb/4+x/bb/aa/2-3.D0/4.D0*x**2/bb/aa**2
  bmat(3,3) = 1/aa/4+y/bb/aa/2-3.D0/4.D0*y**2/bb**2/aa
  bmat(3,4) = -1/bb/aa+3.D0/4.D0*x**2/bb/aa**3+3.D0/4.D0*y**2/bb**3/aa
  bmat(3,5) = 1/bb/4-x/bb/aa/2-3.D0/4.D0*x**2/bb/aa**2
  bmat(3,6) = -1/aa/4-y/bb/aa/2+3.D0/4.D0*y**2/bb**2/aa
  bmat(3,7) = 1/bb/aa-3.D0/4.D0*x**2/bb/aa**3-3.D0/4.D0*y**2/bb**3/aa
  bmat(3,8) = -1/bb/4+x/bb/aa/2+3.D0/4.D0*x**2/bb/aa**2
  bmat(3,9) = -1/aa/4+y/bb/aa/2+3.D0/4.D0*y**2/bb**2/aa
  bmat(3,10) = -1/bb/aa+3.D0/4.D0*x**2/bb/aa**3+3.D0/4.D0*y**2/bb**3/aa
  bmat(3,11) = -1/bb/4-x/bb/aa/2+3.D0/4.D0*x**2/bb/aa**2
  bmat(3,12) = 1/aa/4-y/bb/aa/2-3.D0/4.D0*y**2/bb**2/aa

  ! Compute element strain
  estrain = -z*MATMUL(Bmat, de)

  ! Build constitutive matrix
  Dmat = 0.
  Dmat(1, 1) = young*thk**3./(12.-12.*nu**2.)
  Dmat(1, 2) = nu*young*thk**3./(12.-12.*nu**2.)
  Dmat(1, 3) = 0.
  Dmat(2, 1) = nu*young*thk**3./(12.-12.*nu**2.)
  Dmat(2, 2) = young*thk**3./(12.-12.*nu**2.)
  Dmat(2, 3) = 0.
  Dmat(3, 1) = 0.
  Dmat(3, 2) = 0.
  Dmat(3, 3) = (1./2.)*(1.-nu)*young*thk**3./(12.-12.*nu**2.)
  
  Emat=Dmat*12./thk**3
  ! Compute element stress
  estress = MATMUL(Emat, estrain)

!$$$$$$   ! Compute principal stress and direction
!$$$$$$   sigma1=1/2*(estress(1)+estress(2))+sqrt(((estress(1)-estress(2))/2)*((estress(1)-estress(2))/2)+estress(3)*estress(3))
!$$$$$$   sigma2=1/2*(estress(1)+estress(2))-sqrt(((estress(1)-estress(2))/2)*((estress(1)-estress(2))/2)+estress(3)*estress(3))
!$$$$$$   ! Direction:
!$$$$$$   c = (estress(1)-estress(2))/(sigma1-sigma2)
!$$$$$$   s = (-2*estress(3))/(sigma1-sigma2)
!$$$$$$   psi = datan2(s,c)/2
!$$$$$$   
!$$$$$$   ! Von Mises stress
  !$$$$$$   sigmavm=sqrt(estress(1)*estress(1)+estress(2)*estress(2)-estress(1)*estress(2)+3*estress(3)*estress(3))
  
END SUBROUTINE plate41rect_ss

!!$subroutine shape(xe,bmat)
!!$  ! bmat found analytical
!!$  
!!$  real(8), intent(out) :: bmat(3,12), xe(:)
!!$  REAL(8) :: aa, bb
!!$  
!!$  aa = (xe(4)-xe(1))/2d0
!!$  bb = (xe(11)-xe(2))/2d0
!!$
!!$  bmat(1,1) = 3.D0/4.D0*x/aa**3-3.D0/4.D0*x*y/bb/aa**3
!!$  bmat(1,2) = -1/aa/4+3.D0/4.D0*x/aa**2+y/bb/aa/4-3.D0/4.D0*x*y/bb/aa**2
!!$  bmat(1,3) = 0.
!!$  bmat(1,4) = -3.D0/4.D0*x/aa**3+3.D0/4.D0*x*y/bb/aa**3
!!$  bmat(1,5) = 1/aa/4+3.D0/4.D0*x/aa**2-y/bb/aa/4-3.D0/4.D0*x*y/bb/aa**2
!!$  bmat(1,6) = 0.
!!$  bmat(1,7) = -3.D0/4.D0*x/aa**3-3.D0/4.D0*x*y/bb/aa**3
!!$  bmat(1,8) = 1/aa/4+3.D0/4.D0*x/aa**2+y/bb/aa/4+3.D0/4.D0*x*y/bb/aa**2
!!$  bmat(1,9) = 0.
!!$  bmat(1,10) = 3.D0/4.D0*x/aa**3+3.D0/4.D0*x*y/bb/aa**3
!!$  bmat(1,11) = -1/aa/4+3.D0/4.D0*x/aa**2-y/bb/aa/4+3.D0/4.D0*x*y/bb/aa**2
!!$  bmat(1,12) = 0.
!!$  bmat(2,1) = 3.D0/4.D0*y/bb**3-3.D0/4.D0*y*x/bb**3/aa
!!$  bmat(2,2) = 0.
!!$  bmat(2,3) = -1/bb/4+x/bb/aa/4+3.D0/4.D0*y/bb**2-3.D0/4.D0*y*x/bb**2/aa
!!$  bmat(2,4) = 3.D0/4.D0*y/bb**3+3.D0/4.D0*y*x/bb**3/aa
!!$  bmat(2,5) = 0.
!!$  bmat(2,6) = -1/bb/4-x/bb/aa/4+3.D0/4.D0*y/bb**2+3.D0/4.D0*y*x/bb**2/aa
!!$  bmat(2,7) = -3.D0/4.D0*y/bb**3-3.D0/4.D0*y*x/bb**3/aa
!!$  bmat(2,8) = 0.
!!$  bmat(2,9) = 1/bb/4+x/bb/aa/4+3.D0/4.D0*y/bb**2+3.D0/4.D0*y*x/bb**2/aa
!!$  bmat(2,10) = -3.D0/4.D0*y/bb**3+3.D0/4.D0*y*x/bb**3/aa
!!$  bmat(2,11) = 0.
!!$  bmat(2,12) = 1/bb/4-x/bb/aa/4+3.D0/4.D0*y/bb**2-3.D0/4.D0*y*x/bb**2/aa
!!$  bmat(3,1) = 1/bb/aa-3.D0/4.D0*x**2/bb/aa**3-3.D0/4.D0*y**2/bb**3/aa
!!$  bmat(3,2) = 1/bb/4+x/bb/aa/2-3.D0/4.D0*x**2/bb/aa**2
!!$  bmat(3,3) = 1/aa/4+y/bb/aa/2-3.D0/4.D0*y**2/bb**2/aa
!!$  bmat(3,4) = -1/bb/aa+3.D0/4.D0*x**2/bb/aa**3+3.D0/4.D0*y**2/bb**3/aa
!!$  bmat(3,5) = 1/bb/4-x/bb/aa/2-3.D0/4.D0*x**2/bb/aa**2
!!$  bmat(3,6) = -1/aa/4-y/bb/aa/2+3.D0/4.D0*y**2/bb**2/aa
!!$  bmat(3,7) = 1/bb/aa-3.D0/4.D0*x**2/bb/aa**3-3.D0/4.D0*y**2/bb**3/aa
!!$  bmat(3,8) = -1/bb/4+x/bb/aa/2+3.D0/4.D0*x**2/bb/aa**2
!!$  bmat(3,9) = -1/aa/4+y/bb/aa/2+3.D0/4.D0*y**2/bb**2/aa
!!$  bmat(3,10) = -1/bb/aa+3.D0/4.D0*x**2/bb/aa**3+3.D0/4.D0*y**2/bb**3/aa
!!$  bmat(3,11) = -1/bb/4-x/bb/aa/2+3.D0/4.D0*x**2/bb/aa**2
!!$  bmat(3,12) = 1/aa/4-y/bb/aa/2-3.D0/4.D0*y**2/bb**2/aa
!!$
!!$  
!!$
!!$
!!$end subroutine shape

END MODULE plate41rect
