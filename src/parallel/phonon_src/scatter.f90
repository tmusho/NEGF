!  scatter.f90 
!
!  Authors: G. Walker - Vanderbilt University
!           T. Musho - Vanderbilt University         
!
!  Created: 8/15/10
!  Last Modified: 8/15/10
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
!  PROGRAM: Scatter(Module)
!
!  PURPOSE: Contains subroutine to calculate scattering mixing
!
!  SUBROUTINES: Potent_lin_smix - Hartree Approx, Simple Mixing
!               Potent_initchg - Force convergence of initial electron distribution
!               
!                
!
!****************************************************************************
!
module Scatter

!Use global variables from module
use Startup
use Matxopt
use Anderson_mixing
use Errors
use Debug

implicit none

   type :: ptr_to_array
     type(Anderson_Context), private, pointer :: ST
   end type ptr_to_array
   type(ptr_to_array), private, allocatable, dimension(:) :: SinpT, SoutT

   real(kind=8), private, allocatable, dimension(:) :: Unew, dU

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Scatter_init - Initalize scattering
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Scatter_init

   !Allocation
   select case (coc_s)
     case (1) ! Simple mixing
       !do nothing 
     case (2) ! Anderson mixing
       call Scatter_anderson_allocate
   end select

end subroutine Scatter_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Potent_reset - Reset potential variables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Scatter_reset

   integer(kind=4) :: n

   !Reset anderson mixing table
   if(coc_s .eq. 2) then
     do n = 1, Np
       call Anderson_ResetMix(SinpT(n)%ST)
       call Anderson_ResetMix(SoutT(n)%ST)
     end do
   end if

end subroutine Scatter_reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Scatter_anderson_allocate - Allocate temp vectors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Scatter_anderson_allocate

   integer(kind=4) :: err, n
   character*100 :: var

   allocate(SinpT(Np), SoutT(Np), stat=err)
   write(var,*) 'SinpT(Nd), SoutT(Np)'
   call Errors_Allocate(err,var)

   do n = 1, Np
     call initAnderson(SinpT(n)%ST, Nank_s, NE, smix_s, wmix_s)
     call initAnderson(SoutT(n)%ST, Nank_s, NE, smix_s, wmix_s)
   end do

end subroutine Scatter_anderson_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Scatter_smix - Calculate scattering matricies
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Scatter_lin_smix(SigO,SigI)

   real(kind=8), dimension(:,:), intent(in) :: SigO
   real(kind=8), dimension(:,:), intent(inout) :: SigI

    ! Simple mixing of 
    SigI = ((1-smix_s)*SigI) + (smix_s*SigO)

end subroutine Scatter_lin_smix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Scatter_anderson - Calculate scattering matricies
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Scatter_anderson(SigInpNew_t,SigOutpNew_t,SigInp_t,SigOutp_t)

   integer(kind=4) :: n

   real(kind=8), dimension(:,:), intent(in) :: SigInpNew_t, SigOutpNew_t
   real(kind=8), dimension(:,:), intent(inout) :: SigInp_t, SigOutp_t

    ! Anderson mix
    do n = 1, Np
      call Anderson_Mix(SinpT(n)%ST, SigInp_t(n,:), (SigInpNew_t(n,:)-SigInp_t(n,:)))
      call Anderson_Mix(SoutT(n)%ST, SigOutp_t(n,:), (SigOutpNew_t(n,:)-SigOutp_t(n,:)))
    end do

end subroutine Scatter_anderson

end module
