!  utility.f90 
!
!  Authors: G. Walker - Vanderbilt University
!           T. Musho - Vanderbilt University
!
!  Created: 3/30/09
!  Last Modified: 3/30/09
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
!  PROGRAM: Utility(Module)
!
!  PURPOSE: Contains debug routines
!
!  SUBROUTINES: 
!               
!                
!
!****************************************************************************
!
module Utility

implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Utility_SetCycle - Add step number to filename
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Utility_SetCycle(filename,start,istepno,idigit)
!convert the integer istepno to a string and replace last 
!idigit characters filename with it.
implicit none

!  ----- arguments -----
   character(len=*), intent(inout) :: filename  ! target filename
   integer(4), intent(in) :: istepno             ! step number
   integer(4), intent(in) :: start               ! starting character number
   integer(4), intent(in) :: idigit              ! number of digits in step number

!  ----- local variables -----
   integer(4) :: icount,ipos,inum

   icount=istepno
   do ipos=start+idigit-1,start,-1
      inum=mod(icount,10)
      filename(ipos:ipos)=char(ichar('0')+inum)
      icount=icount/10
   enddo

   return
end subroutine Utility_SetCycle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Utility_MakeFileName - Makes filename
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character(len=80) function Utility_MakeFileName(Prefix,NumFluidProcs,step,ID,PostFix)
   character(len=*), intent(in) :: PreFix
   integer(4), intent(in) :: NumFluidProcs
   integer(4), intent(in) :: step
   integer(4), intent(in) :: ID
   character(len=*), intent(in) :: PostFix

   character(len=80) :: FileName   
   integer(4) :: Len1

   Len1=len(Prefix)
   
   if(NumFluidProcs==1)then
      FileName=Prefix//'_xxxxxx.'//PostFix
      call Utility_SetCycle(FileName,Len1+2,step,6)  
   else
      FileName=Prefix//'_xxxxxx_xxx.'//PostFix
      call Utility_SetCycle(FileName,Len1+2,step,6) 
      call Utility_SetCycle(FileName,Len1+9,ID,3)
   endif
   
   Utility_MakeFileName=FileName        

   return
end function Utility_MakeFileName

end module


