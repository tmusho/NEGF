!  energy_amr.f90 
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
!  PURPOSE: Contains subroutine to calculate energy amr
!
!  SUBROUTINES: Energy_amr - 
!           
!              
!
!****************************************************************************
!
module Energy_amr

use Startup

implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Energy_amr_refine_log - AMR refine based on log ratio of max value
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Energy_amr_refine_log(IT,Er)

  integer(kind=4) :: k, n
  real(kind=8) :: Imv
  real(kind=8), dimension(:), intent(in) :: IT
  integer(kind=4), dimension(:), intent(inout) :: Er

  Imv = maxval(abs(IT)) ! Determine max current transmission
  Er(:) = 1 !Set refinement to one
  do k = 1, NE
    do n = 1, amr_d
      if(abs(log(abs(IT(k))/Imv)) .le. real(n)) then
        Er(k) = Er(k) + 1
      end if
    end do
  end do

end subroutine Energy_amr_refine_log
end module
