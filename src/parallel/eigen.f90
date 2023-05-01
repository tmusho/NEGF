!  eigen.f90 
!
!  Authors: G. Walker - Vanderbilt University
!           T. Musho - Vanderbilt University         
!
!  Created: 12/09/08
!  Last Modified: 12/18/08
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, version 3 of the License.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!
!****************************************************************************
!
!  PROGRAM: Eigen(Module)
!
!  PURPOSE: Contains subroutine to determine eigenvalues of Hamiltonian
!
!  SUBROUTINES: Eigen_Main - Main Startup, subroutine calls
!               Eigen_Plot - Plot Eigenvalues to Eig.dat
!               Eigen_Calc - Calculate Eigenvalues - Lapack Call
!               Eigen_Allocate - Allocate arrays, define LAPACK variables 
!
!****************************************************************************
!
module Eigen

!Use global variables from module
use Startup
use Device
use Hamil

implicit none

integer(kind=8) :: LWORK, LDA
integer(kind=8), private, allocatable, dimension(:) :: IPIV
real(kind=8), allocatable, dimension(:) :: DH, DHO, DH_p
complex(kind=8), private, allocatable, dimension(:) :: WORK

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Eigen_Main - Main Startup, subroutine calls
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Eigen_Main
  
   call Eigen_Allocate ! Allocate device matrix
   call Eigen_Fill
   if(rank .eq. 0) then
      call Eigen_Calc ! Input values into device matrix
      !call Eigen_Shift
      if(debug_l>5) call Eigen_Plot ! Plot hamiltonian diag and raw
   end if

   deallocate(DHO,WORK,IPIV)

end subroutine Eigen_Main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Eigen_Allocate - Allocate arrays, define LAPACK variables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Eigen_Allocate

   integer(kind=4) :: err
   character*120 :: var

   allocate(DH(Np),DHO(Np), stat=err) ! Input - Diagonal Values, Output - Allocation Eigenvalue Vector
   write(var,*) 'DH(Np), DHO(Np)'
   call Errors_Allocate(err,var)

  !! LAPACK Variables !!
   allocate(WORK(4*Np), IPIV(Np), stat=err) !Allocate LAPACK Variables (Work will be reallocated in Main Loop for Optimal Size)
   write(var,*) 'WORK(4*Np), IPIV(Np)'
   call Errors_Allocate(err,var)!Allocate LAPACK Variables - Leading Dimension of Array

   LWORK = 4*Np !Allocate LAPACK Variables - LWORK >= max(1,N) Optimal Work Calculated below in Main Loop
   LDA = 1
   DH(:) = 0.0
   DHO(:) = 0.0

end subroutine Eigen_Allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Eigen_Fill - Fill LAPACK Arrays
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Eigen_Fill

   integer(kind=4) :: r

   !! Put Diagonal and Offdiagonal Tridiagonal Part of Hamiltonian Matrix into Arrays !!
   do r = 1, Np
      DH(r) = H(r,r) ! Diagonal Terms from Hamiltonian Matrix
      if (r .ne. Np) then
         DHO(r) = H(r,r+1) ! Offdiagonal Terms from Hamiltonian Matrix
      end if
   end do

end subroutine Eigen_Fill

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Eigen_Shift - Fill LAPACK Arrays
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Eigen_Shift

    integer(kind=4) :: r, err
    character*120 :: var

    DHO(1) = Ecs

    if(Np>Nb) then
      do r = 2, Nb+1
        DHO(r) = DH(r-1)
      end do
    else
      do r = 2, Np
        DHO(r) = DH(r-1)
      end do
      Nb = Np
    end if

    deallocate(DH)
    allocate(DH(Nb), stat=err)
    write(var,*) 'DH(Nb)'
    call Errors_Allocate(err,var)

    DH(:) = 0.0

    do r = 1, Nb
      DH(r) = DHO(r)
    end do


end subroutine Eigen_Shift


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Eigen_Calc - Calculate Eigenvalues
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Eigen_Calc

    integer(kind=4)::INFO

    !! Call LAPACK Routine to Determine Eigenvalues of the Tridiagonal Hamiltonian Matrix !!
    call DSTEV('N',Np,DH,DHO,IPIV,LDA,WORK,INFO) ! Calculates energy eigenvalues of Hamiltonian (Eigenvectors are not calculated) 
    if (INFO <  0) then
       write (*,*) 'the i-th argument had an illegal value'
       stop 'abort'
    else if (INFO > 0) then
       if (INFO <= Np) then
          write (*,*) 'the i-th argument had an illegal value'
       else
          write (*,*) 'the algorithm failed to converge i off-diagonal elements of E did not converge to zero.'
       end if 
       stop 'abort'
    end if

end subroutine Eigen_Calc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Eigen_Plot - Plot Eigenvalues
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Eigen_Plot

   integer(kind=4) :: r, err
   character*120 :: var

   open(unit=13, file='Eig.dat', status='new', action='write', iostat=err)
   write(var,*) 'Eig.dat'
   call Errors_Fileopen(err,var)

   do r = 1, Nb
      write(13,'(I2,E15.6)') r, DH(r)
   end do

   close(unit=13)

end subroutine Eigen_Plot

end module
