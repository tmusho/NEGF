!  Input.f90 
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
!  PROGRAM: Input(Module)
!
!  PURPOSE: Contains input routines
!
!  SUBROUTINES: 
!              
!            
!
!****************************************************************************
!
module Input

use Utility, only : UnitIn, &
                    UnitOut, &
                    StdOut, &
                    UnitLog, &
                    LINE_LENGTH
use Utility, only : Utility_ErrorMsgChar,Utility_ErrorMsg, &
                    Utility_InitErrorMsg,Utility_Inputerror

implicit none

character(LINE_LENGTH), dimension(:), pointer :: InputLine
integer(4) :: NumLines

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Input_Load - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Input_Load(FileName,FluidID,NumFluidProcs)
!  Input_Load - check format code and read input

   character(len=*), intent(in) :: FileName
   integer(4), intent(in) :: FluidID
   integer(4), intent(in) :: NumFluidProcs
         
   !  --- local data ---
   character(LINE_LENGTH) :: NewLine 
   integer(4) :: ios,error,LineNo
   logical(4) :: here

!  Open Input File   
   inquire(file=FileName,exist=here)
   if(.not.here)then
      call Utility_ErrorMsg('Input_Load::input file '//trim(Filename)//' missing.')
   endif
      
   open(UnitIn,file=Filename,action='READ')  

   ! allocate space   
   NumLines=-1
   ios=0
   do while (ios==0)
      NumLines=NumLines+1     
      NewLine=Input_NextLine(ios)      
   enddo
   rewind(UnitIn)  

   allocate(InputLine(NumLines),STAT=error)
   if(error/=0)then
      call Utility_InitErrorMsg('Input_Load:: allocation failed',FluidID,NumFluidProcs)
   endif

   !  read and store input
   LineNo=0
   NewLine=Input_NextLine(ios)
   do while (ios==0)
      LineNo=LineNo+1
      InputLine(LineNo)=NewLine
      NewLine=Input_NextLine(ios)
   enddo

   close(UnitIn)

   return
end subroutine Input_Load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Input_Echo - echo input 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Input_Echo(FileName,FluidID)

   character(len=*), intent(in) :: FileName
   integer(4), intent(in) :: FluidID

!   --- local vars ---
   logical(4):: here
   integer(4) :: ios
   character(LINE_LENGTH) :: NewLine 

!  Open Input File   
   inquire(file=FileName,exist=here)
   if(.not.here)then
      call Utility_ErrorMsg('Input_Echo:input file '//trim(Filename)//' missing.')
   endif
      
   open(UnitIn,file=Filename,action='READ')  

   if(FluidID==0)then

   !  echo input file
      write(UnitOut,*)
      write(UnitOut,*)'        +++++ INPUT FILE +++++ '
      read(UnitIn,'(1a120)',IOSTAT=ios)NewLine
      do while (ios==0)
         write(UnitOut,'(1a120)')NewLine
         read(UnitIn,'(1a120)',IOSTAT=ios)NewLine
      enddo
      write(UnitOut,*)'      +++++ END INPUT FILE +++++ '
      write(UnitOut,*)      
      rewind(UnitIn)

   endif

   flush(UnitOut)
   close(UnitIn)
   
   return
   
end subroutine Input_Echo 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Input_FindSection - Locates the start and end line number of section 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Input_FindSection(echo,string1,string2,MinLines,MaxLines,NumLines,StartLineNo,fatal)

   logical(4), intent(in) :: echo
   character(*), intent(in) :: string1,string2   ! target string
   integer(4), intent(in) :: MinLines,MaxLines
   integer(4), intent(out) :: NumLines,StartLineNo
   logical(4), intent(in) :: fatal       ! if .true., abort if word not found

   integer(4) :: EndLineNo

   StartLineNo=Input_FindKeyWord(string1,echo,fatal)
   if(StartLineNo>0)then
   !  Key word found 
      EndLineNo=Input_FindKeyWord(string2,.false.,.true.)
      NumLines=EndLineNo-StartLineNo-1
      if(NumLines>0)then
         if((NumLines<MinLines).or.(NumLines>MaxLines))then
            call Utility_ErrorMsgChar('Input_FindSection:: Wrong number of lines in section ',string1)
         endif
      endif
   else
   !  Key word missing
      NumLines=0
   endif

   return
end subroutine Input_FindSection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Input_GetLine - returns specified line
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character(LINE_LENGTH) function Input_GetLine(LineNo)

   integer(4), intent(in) :: LineNo

   if((LineNo>=1).and.(LineNo<=NumLines))then
      Input_GetLine=InputLine(LineNo)
   else
      call Utility_ErrorMsg('Input_GetLine:: Input Line Number out of bounds')
   endif

   return
end function Input_GetLine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Input_CheckLine - check next input line for character string
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer(4) function Input_FindKeyWord(string,echo,fatal)

implicit none

!  ----- arguments -----
   character(*), intent(in) :: string   ! target string
   logical(4), intent(in) :: echo        ! echo results flag
   logical(4), intent(in) :: fatal       ! if .true., abort if word not found

!  ----- local variables -----
   integer(4) :: LineNo

   LineNo=1
   do while((LineNo<=NumLines).and.(index(InputLine(min(LineNo,NumLines)),string)==0))
      LineNo=LineNo+1
   enddo
   if(LineNo<=NumLines)then
   ! line found 
      if(echo)then
         write(UnitOut,*)
         write(UnitOut,*)string
      endif
      Input_FindKeyWord=LineNo
   else if(.not.fatal)then
   !  not found & not fatal
      Input_FindKeyWord=0
   else
   ! not found & fatal
      call Utility_ErrorMsgChar('Input_KeyWordLine :: Required Key Phrase Not Found: ',string)
   endif    

   return

end function Input_FindKeyWord

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Input_NextLine - returns the next valid line; skips blank lines, comment 
!                  lines and removes comments
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
character(LINE_LENGTH) function Input_NextLine(ios)

implicit none
   integer(4), intent(out) :: ios ! input unit number

   character(LINE_LENGTH) :: NewLine
   logical :: done
   integer(4) :: pound,i,length,j

   length=len(NewLine)
   done=.false.
   ios=0
   do while ((ios==0).and.( .not. done))
   !   +++++ change 1a120 to match line length ++++++
       read(UnitIn,'(1a120)',IOSTAT=ios)NewLine
       if(ios==0)then
       !  not EOF
          call RemoveCR(NewLine)          
          call Input_Upper2Lower(NewLine)
          pound=index(NewLine,'#')
          if(pound==0)then
          !  is not a comment card; check to see if its blank
             i=0
             do j=1,Length
                if(NewLine(j:j)/=' ') i=i+1
             enddo
             if(i>0) then
             ! not a blank line
               done=.true.
             endif
          else
          !  contains comments
             i=1
             ! check if characters other than blanks appear before #
             do while ((i < pound).and.(NewLine(i:i)==' '))
                i=i+1
             enddo
             if(i<pound)then
             !  characters before #, not a comment card
                done=.true.
             endif
          endif
       else
          pound=0
       endif   
    enddo
    if(pound/=0)then
       do i=pound,LINE_LENGTH
          NewLine(i:i)=' '
       enddo
    endif
    Input_NextLine=NewLine

   return

end function Input_NextLine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Upper2Lower - converts all characters on line from upper case to lower
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Input_Upper2Lower(Line)
   character(len=*), intent(inout) :: Line

   integer(4), parameter :: Cap2Small=iachar("a")-iachar("A")
   integer(4) :: i,length

   length=len(Line)
   do i=1,Length
      if((line(i:i)>='A').and.(line(i:i)<='Z'))then
      !  upper case letter
         line(i:i)=achar(iachar(line(i:i))+Cap2Small)
      endif
   enddo

   return

end subroutine Input_Upper2Lower

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! RemoveCR - Removes Carriage Returns from string
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RemoveCR(String)    
!  DOS  EOR: LF CR
!  unix EOR: LF    
   character(*),intent(inout):: String  !The string to be transformed
   integer(4):: i,Length
  
   Length=len(String)
   do i=Length,1,-1
     if(String(i:i)==achar(13))then
        String(i:)=String(i+1:)//' '
     endif 
   enddo

   return
end subroutine RemoveCR 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Input_Fread - return real8 input from line and checks that number is within bounds
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(8) function Input_Fread(LineNo,word,xlow,xhigh)

implicit none      

!  ----- arguments -----
   integer(4), intent(in) :: LineNo
   character(len=*), intent(in) :: word  ! input variable title
   real(8), intent(in) :: xlow,xhigh      ! upper and lower bounds for fnum

!  -----local variable-----
   integer :: ios
   real(8) :: fnum

   read (InputLine(LineNo),*,IOSTAT=ios)fnum
   if(ios/=0)then
      call Utility_InputError(ios,'Input_Fread','at variable '//word,InputLine(LineNo))
   endif
   write(UnitOut,*)word,' = ',fnum

   if(xhigh.gt.xlow)THEN
     if(fnum.lt.xlow.or.fnum.gt.xhigh)then
        call Utility_ErrorMsgChar('Input_Fread:: input out of bounds, variable = ',word)
     endif
   endif

   Input_Fread=fnum

   return
end function Input_Fread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Input_Iread - return integer input from line and checks that number is within bounds
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer(4) function Input_Iread(LineNo,word,ilow,ihigh)

!  -----arguments -----
   integer(4), intent(in) :: LineNo
   character(len=*), intent(in) :: word    ! input variable title
   integer(4), intent(in) :: ilow,ihigh     ! upper and lower bounds for inum

!  -----local variable----- 
   integer(4) :: inum        ! input variable 
   integer(4) :: ios         ! read status
   real(4) :: rnum

   read(InputLine(LineNo),*,IOSTAT=ios)inum                         
   if(ios/=0)then
   !  integer read failed, try real read
      read(InputLine(LineNo),*,IOSTAT=ios)rnum 
      if(ios/=0)then                             
         call Utility_InputError(ios,'Input_Iread','at variable '//word,InputLine(LineNo))
      endif
      inum=NINT(rnum)
   endif

   write(UnitOut,*)word,' = ',inum

   if(ihigh.gt.ilow)then
      if(inum.lt.ilow.or.inum.gt.ihigh)then
         call Utility_ErrorMsgChar('Input_Iread:: input out of bounds, variable = ',word)
      endif
   endif

   Input_Iread=inum
   return
end function Input_Iread

end module
