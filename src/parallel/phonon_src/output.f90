!  output.f90
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
!  PROGRAM: Output(Module)
!
!  PURPOSE: Contains subroutine to output data to text file
!
!  SUBROUTINES: 
!               
!               
!                
!
!****************************************************************************
!
module Output

!Use global variables from module
use Startup
use Errors
use Debug

implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Output_iv - Plot ASCII output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Output_iv(IV)

    integer(kind=4) :: n, err
    character*100 :: var

    real(kind=8), dimension(:,:), intent(in) :: IV

    open (unit=28, file='IV.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'IV.dat'
    call Errors_Fileopen(err,var)

    !IV Plot
    do n = 1, NV
      write(28,*) IV(n,1), -IV(n,2)
    end do 

    close(28)
end subroutine Output_iv
end module
