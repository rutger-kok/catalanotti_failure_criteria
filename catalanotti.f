! This is a subroutine implementing a set of failure criteria for
! composite materials.
! Copyright (C) 2021 Rutger Wouter Kok

! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.

! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.

! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
! 02110-1301 USA

! These subroutines are an implementation of the three-dimensional
! failure criteria for composite materials developed by Catalanotti et
! al. in [1]. The subroutines can be used as part of an Abaqus
! subroutine, called from a standalone Fortran program, or called from
! a Python script using F2PY.

! [1] G. Catalanotti, P.P. Camanho, A.T. Marques
! Three-dimensional failure criteria for fiber-reinforced
! laminates
! Composite Structures 95 (2013) 63–79
! http://dx.doi.org/10.1016/j.compstruct.2012.07.016   
      
      subroutine fail_initiation(trialStress,ST,SL,etaL,etaT,lambda,
     1                           kappa,FI,a_fail)
        ! Purpose: Iterates over failure plane angles from 0 to 180
        !   degrees to determine the failure angle and max failure
        !   index (using Catalanotti failure criteria)
        ! Variable Dictionary:
        ! FI = failure index
        ! a = angle of failure plane in degrees
        ! aR = angle of failure plane in radians
        implicit none
        ! input variables
        real*8, dimension(6), intent(in) :: trialStress
        real*8, intent(in) :: ST,SL,etaL,etaT,lambda,kappa
        ! local variables
        real*8 :: pi,aR,trialFI
        integer :: a
        ! output variables
        real*8, intent(out) :: FI,a_fail
    
        pi = 4.0d0*atan(1.0d0) ! determine value of pi
    
        FI = 0.0d0 ! initialize failure index
        a = 0 ! initial failure plane angle
        ! iterate over angles to determine angle which maximizes FIs
        do while (a <= 180) ! iterate over angles from 0 to 180 degrees
          aR = a*(pi/180.0d0)  ! angle in radians
          call catalanotti(trialStress,ST,SL,etaL,etaT,lambda,kappa,
     1                     aR,trialFI)
          ! Update max value if current value > max value
          if (trialFI > FI) then
            FI = trialFI
            a_fail = aR ! record failure plane 
          end if
          ! Update angle
          a = a + 1
        end do
      end subroutine fail_initiation

      subroutine catalanotti(trialStress,ST,SL,etaL,etaT,lambda,kappa,
     1                       aR,trialFI)
           ! Purpose: Failure criteria from Catalanotti et al. (2013)
           ! Variable Dictionary:
           ! trialStress = stress calculated using the undamaged
           !   stiffness tensor.
           ! trialStress_crack = trialStress in crack coordinate system
           ! etaT = friction coefficient in the transverse direction
           ! etaL = friction coefficient in the longitudinal direction
           ! kappa = parameter used to calculate failure indices
           ! lambda = parameter used to calculate failure indices
           ! ST = in-situ transverse shear strength
           ! SL = in-situ longitudinal shear strength
           ! trialFI_MT = stores trial values of tensile failure index
           ! trialFI_MC = stores trial values of compression failure index
           ! trialFI = max of tensile and compressive failure indices
           ! tN, tT, tL = tractions on failure plane
           ! aR = angle in radians
           ! R_GC = transformation matrix from global to crack (failure
           !   plane) coordinate system
           implicit none
           ! input variables
           real*8, dimension(6), intent(in) :: trialStress
           real*8, intent(in) :: ST,SL,etaL,etaT,lambda,kappa,aR
           ! local variables
           real*8 :: tN,tT,tL,trialFI_MT,trialFI_MC
           real*8, dimension(6,6) :: R_GC
           real*8, dimension(6) :: trialStress_crack
           ! output variables
           real*8, intent(out) :: trialFI
           call rot_matrix(aR, R_GC)
           ! rotate stresses to crack coordinate frame
           trialStress_crack = matmul(R_GC(1:6,1:6), trialStress(1:6))
           ! define tractions
           tN = trialStress_crack(2)
           tT = trialStress_crack(5)
           tL = trialStress_crack(4)
           ! Calculate value of failure indices at current angle
           if (tN >= 0.0d0) then
             trialFI_MT = (tN/ST)**2.0d0 + (tL/SL)**2.0d0 +
     1                 (tT/ST)**2.0d0 + lambda*(tN/ST)*(tL/SL)**2.0d0 +
     2                 kappa*(tN/ST) ! Eq. 42 [1]
             trialFI_MC = 0.0d0
           else
             trialFI_MC = (tL/(SL-etaL*tN))**2.0d0 +
     1                 (tT/(ST-etaT*tN))**2.0d0 ! Eq. 5 [1]
             trialFI_MT = 0.0d0
           end if
           ! Return the maximum trial failure index
           trialFI = max(trialFI_MT,trialFI_MC)
        end subroutine catalanotti

      subroutine rot_matrix(angle, R)
        ! Purpose: defines transformation matrix to rotate from failure
        !   plane coordinates system and back.
        ! Variable Dictionary:
        ! angle = angle through which to rotate
        ! R = transformation matrix
        ! m = cosine of angle
        ! n = sine of angle
        implicit none
        ! input variables
        real*8, intent(in) :: angle
        ! local variables
        real*8 m,n
        ! output variables
        real*8, dimension(6,6), intent(out) :: R
        m = cos(angle)
        n = sin(angle)
        R = 0.0d0
        R(1,1) = 1.0d0
        R(2,2) = m**2.0d0
        R(2,3) = n**2.0d0
        R(2,5) = 2.0d0*m*n
        R(3,2) = n**2.0d0
        R(3,3) = m**2.0d0
        R(3,5) = -2.0d0*m*n
        R(4,4) = m
        R(4,6) = n
        R(5,2) = -m*n
        R(5,3) = m*n
        R(5,5) = (m**2.0d0) - (n**2.0d0)
        R(6,4) = n
        R(6,6) = m
      end subroutine rot_matrix

      subroutine rotate_stress(trialStress,phiC,XC,g12,trialStressP)
        ! Purpose: Rotate stresses to misaligned coordinate for
        !   calculation of longitudinal compressive failure index.
        ! Variable dictionary:
        ! trialStress = stress calculated using the undamaged stiffness
        !   tensor.
        ! trialStressP = trialStress rotated into the misalignment
        !   coordinate frame to calculate the failure index for
        !   longitudinal compression.
        ! trialStressT = trialStress rotated by in-plane angle theta.
        ! phiC = initial misalignment angle for compressive failure
        !   calculation
        ! g12 = longitudinal compressive failure stress
        ! phi = angle of kink band
        ! theta = kinking plane angle
        ! phi0 = initial misalignment angle
        ! m = cosine of theta
        ! n = sine of theta
        ! u = cosine of phi
        ! v = sine of phi
        ! tS4EQ0 = boolean - checks if trialStress(4) = 0
        ! tS6EQ0 = boolean - checks if trialStress(6) = 0
        ! gammaM = shear stress induced misalignment angle
        ! gammaMC = gammaM under axial compression loading
        ! eps = tolerance for conditional statements
        implicit none
        ! input variables
        real*8, intent(in) :: phiC,XC,g12
        real*8, dimension(6), intent(in) :: trialStress
        ! local variables
        real*8, dimension(6) :: trialStressT
        real*8  m,n,u,v,phi,theta
        real*8  gammaMC,gammaM,phi0
        real*8, parameter :: eps=1.d-8
        logical tS4EQ0, tS6EQ0
        ! output variables
        real*8, dimension(6), intent(out) :: trialStressP

        ! first determine kink plane angle theta
        tS4EQ0 = (abs(trialStress(4)-0.d0) < eps)
        tS6EQ0 = (abs(trialStress(6)-0.d0) < eps)

        if (tS4EQ0.and.tS6EQ0) then
          if (abs(trialStress(2)-trialStress(3)) < eps) then
            theta = atan(1.0d0) ! pi/4
          else
            theta = 0.5d0*atan((2.0d0*trialStress(5))/(trialStress(2)-
     1              trialStress(3))) !Eq.55 [1]
          end if 
        else
          if (tS4EQ0) then
            theta = 2.0d0*atan(1.0d0) ! pi/2
          else 
            theta = atan(trialStress(6)/trialStress(4)) !Eq. 56 [1]
          end if
        end if

        ! Rotate stresses by angle theta
        m = cos(theta)
        n = sin(theta)
        trialStressT(1) = trialStress(1)
        trialStressT(2) = trialStress(2)*m**2 +
     1                    2.0d0*trialStress(5)*m*n +
     2                    trialStress(3)*n**2
        trialStressT(3) = trialStress(3)*m**2 -
     1                    2.0d0*trialStress(5)*m*n +
     2                    trialStress(2)*n**2
        trialStressT(4) = trialStress(4)*m + trialStress(6)*n
        trialStressT(5) = trialStress(5)*(m**2 - n**2) -
     1                    trialStress(2)*n*m + trialStress(3)*n*m
        trialStressT(6) = trialStress(6)*m - trialStress(4)*n

        ! Determine kink band angle phi
        gammaMC = (sin(2.0d0*phiC)*XC)/(2.0d0*g12)  ! Eq. 74 [1]
        phi0 = phiC - gammaMC  ! Eq. 75
        gammaM = ((phi0*g12 + abs(trialStressT(4)))/
     1           (g12+trialStressT(1)-trialStressT(2)) - phi0)  ! Eq. 81 [1]
        phi = sign(1.0d0, trialStress(4))*(phi0 + gammaM)  ! Eq. 77 [1]

        ! Rotate stresses by angle phi
        u = cos(phi)
        v = sin(phi)
        trialStressP(1) = trialStressT(1)*u**2 +
     1                    2.0d0*trialStressT(4)*u*v +
     2                    trialStressT(2)*v**2
        trialStressP(2) = trialStressT(2)*u**2 -
     1                    2.0d0*trialStressT(4)*v*u +
     2                    trialStressT(1)*v**2     
        trialStressP(3) = trialStressT(3)      
        trialStressP(4) = trialStressT(4)*(u**2 -v**2) +
     1                    trialStressT(2)*v*u - trialStressT(1)*v*u
        trialStressP(5) = trialStressT(5)*u - trialStressT(6)*v     
        trialStressP(6) = trialStressT(6)*u + trialStressT(5)*v

        end subroutine rotate_stress