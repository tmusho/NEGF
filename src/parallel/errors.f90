!  errors.f90 
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
!  PROGRAM: Errors(Module)
!
!  PURPOSE: Contains subroutine to display errors
!
!  SUBROUTINES: Errors_Allocate - Allocation Error
!               Errors_Deallocate - Deallocation Error
!               Errors_Fileopen - File open error
!
!****************************************************************************
!
module Errors

implicit none

integer(kind=4) :: num_w = 0

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Errors_Allocate - Allocation Error
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Errors_Allocate(inerr,var)

  integer(kind=4), intent(in) :: inerr
  character*100, intent(in) :: var

  if(inerr .ne. 0) then
       write(6,*) 'Error Allocation Failed - ',var
       stop 'abort'
  end if

end subroutine Errors_Allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Errors_Deallocate - Deallocation Error
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Errors_Deallocate(inerr,var)

  integer(kind=4), intent(in) :: inerr
  character*100, intent(in) :: var

  if(inerr .ne. 0) then
       write(6,*) 'Error Deallocation Failed - ',var
       stop !call MPI_Abort(MPI_COMM_WORLD,3)
  end if

end subroutine Errors_Deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Errors_Fileopen - File open error
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Errors_Fileopen(inerr,var)

  integer(kind=4), intent(in) :: inerr
  character*100, intent(in) :: var

  if(inerr .ne. 0) then
       write(6,*) 'Error Opening File - ',var
       stop 'abort'
  end if

end subroutine Errors_Fileopen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Errors_HamilE - Check hamiltonian energy 
! ref pg. 159 Datta
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Errors_HamilE(t,Ec,Ef)

  real(kind=8), intent(in) :: t, Ec, Ef

  if(abs(Ef-Ec) .gt. t .and. num_w .eq. 0) then
    write(6,*) ' !! Warning grid spacing too large - Ef-Ec > to !!'
    write(6,*) ' !!Ef-Ec = ',(Ef-Ec),' t = ',t,' !!'
    num_w = num_w + 1
  end if

end subroutine Errors_HamilE

end module
