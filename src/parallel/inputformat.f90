!  InputFormat.f90 
!
!  Authors: G. Walker - Vanderbilt University
!           T. Musho - Vanderbilt University
!
!  Created: 1/20/11
!  Last Modified: 1/20/11
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
!  PROGRAM: InputFormat(Module)
!
!  PURPOSE: Checks Input Format
!
!  SUBROUTINES: 
!              
!            
!
!****************************************************************************
!
module InputFormat
use Utility, only : LINE_LENGTH,Utility_InitErrorMsg
use Input, only: Input_GetLine
implicit none

!  Format Versions
character(12), private :: slpfkey='format 00003'

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! InputFormat_Check - Check Format
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine InputFormat_Check(Code)
!  Input_Load - check format code and read input

   character(len=*), intent(in) :: Code

   character(LINE_LENGTH) :: NewLine 

   NewLine=Input_GetLine(1)

   if(index(Code,'slpf')/=0)then
      if(index(NewLine,slpfkey)==0)then
         call Utility_InitErrorMsg('Bad input file Format Key - correct key is: '//slpfkey,0,1)
      end if
   end if

   return
end subroutine InputFormat_Check

end module InputFormat
