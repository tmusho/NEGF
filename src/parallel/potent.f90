!  potent.f90 
!
!  Authors: G. Walker - Vanderbilt University
!           T. Musho - Vanderbilt University         
!
!  Created: 5/15/09
!  Last Modified: 5/15/09
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
!  PROGRAM: Potent(Module)
!
!  PURPOSE: Contains subroutine to calculate potential
!
!  SUBROUTINES: Potent_lin_smix - Hartree Approx, Simple Mixing
!               Potent_initchg - Force convergence of initial electron distribution
!               
!                
!
!****************************************************************************
!
module Potent

!Use global variables from module
use Startup
use Device
use Matxopt
use Anderson_mixing
use Errors
use Debug

implicit none

   type(Anderson_Context), private, pointer :: PT
   real(kind=8), private, allocatable, dimension(:) :: Unew, dU

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Potent_init - Initialize potential
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Potent_init

   !Allocation
   select case (coc)
     case (1) ! Simple mixing
       call Potent_smix_allocate 
     case (2) ! Anderson mixing
       call Potent_anderson_allocate
   end select

end subroutine Potent_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Potent_reset - Reset potential variables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Potent_reset

   dU(:) = 0.0; Unew(:) = 0.0;
   !Reset anderson mixing table
   if(coc .eq. 2) call Anderson_ResetMix(PT)

end subroutine Potent_reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Potent_smix_allocate - Allocate temp vectors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Potent_smix_allocate

   integer(kind=4) :: err
   character*120 :: var

   allocate(Unew(Np), dU(Np), stat=err)
   write(var,*) 'Unew(Nd), dU(Np)'
   call Errors_Allocate(err,var)

   Unew(:) = 0.0; dU(:) = 0.0

end subroutine Potent_smix_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Potent_anderson_allocate - Allocate temp vectors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Potent_anderson_allocate

   integer(kind=4) :: err
   character*120 :: var

   allocate(Unew(Np), dU(Np), stat=err)
   write(var,*) 'Unew(Nd), dU(Np)'
   call Errors_Allocate(err,var)

   Unew(:) = 0.0; dU(:) = 0.0

  call initAnderson(PT, Nank, Np, smix, wmix)

end subroutine Potent_anderson_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Potent_smix_deallocate - Deallocate temp vectors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Potent_deallocate

   integer(kind=4) :: err
   character*120 :: var

   deallocate(dU, Unew, stat=err)
   write(var,*) 'Unew(Nd), dU(Np)'
   call Errors_Allocate(err,var)

end subroutine Potent_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Potent_lin_smix - Hartree Approx, Simple Mixing
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Potent_lin_smix(NN,chg,Ud,No)

   real(kind=8), intent(out) :: chg
   real(kind=8), dimension(:), intent(in) :: NN
   real(kind=8), dimension(:), intent(inout) :: Ud
   real(kind=8), dimension(:), intent(inout) :: No

   select case (cop)
     case (0) ! Hartree-Fock approx
       Unew = -Uo*(NN-No) ! Channel potential (eV)
     case (1) ! Full poisson potential, more exspensive
       write(*,*) !! Full potential function not implemented !!'
       !Unew = -matmul(Uo_f,(NN-No))!+matmul(Uo_bdy,Ubdy)
   end select

   dU = Unew-Ud ! Change in channel potential due to new electron density 
   chg = maxval(dabs(dU)) ! Find Max Absolute Value
   Ud = Ud + smix*dU ! Updating channel potential with respect to old potential               

   if(debug_l>2) write(6,'(/,A,ES11.4)') 'Potential chg = ',chg ! Write change of charge to screen

   call Potent_initchg(V,NN,chg,No)

end subroutine Potent_lin_smix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Potent_anderson - Anderson Mixing potential calc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Potent_anderson(NN,chg,Ud,No)

   real(kind=8), intent(out) :: chg
   real(kind=8), dimension(:), intent(in) :: NN
   real(kind=8), dimension(:), intent(inout) :: Ud
   real(kind=8), dimension(:), intent(inout) :: No

   select case (cop)
     case (0) ! Hartree-Fock approx
       Unew = -Uo*(NN-No) ! Channel potential (eV)
     case (1) ! Full poisson potential, more exspensive
       write(*,*) !! Full potential function not implemented !!'
       !Unew = -matmul(Uo_f,(NN-No))!+matmul(Uo_bdy,Ubdy)
   end select

   dU = Unew-Ud ! Change in channel potential due to new electron density 
   call Anderson_Mix(PT, Unew, dU)
   chg = maxval(dabs(dU)) ! Find Max Absolute Value               
   Ud = Unew

   if(debug_l>2) write(6,'(/,A,ES11.4)') 'Potential chg = ',chg ! Write change of charge to screen

   call Potent_initchg(V,NN,chg,No)

end subroutine Potent_anderson

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Potent_initchg - Force convergence of initial electron distribution
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Potent_initchg(V,NN,chg,No)

   real(kind=8), intent(inout) :: chg
   real(kind=8), intent(in) :: V
   real(kind=8), dimension(:), intent(in) :: NN
   real(kind=8), dimension(:), intent(out) :: No

   if (V .eq. 0.0) then !First Iteration Enter (It was V==0 which is not definite because V is Floating Point)
     if(chg .le. sc_ip) then
        No = NN
        chg = -1.0 ! Exit loop
     else
        chg = 1.0 ! Stay in loop
     end if
   end if

end subroutine Potent_initchg

end module
