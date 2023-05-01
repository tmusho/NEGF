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
integer(4), parameter :: UnitStop=1         ! UnitStop
integer(4), parameter :: UnitDebug=2        ! Debug Output 
integer(4), parameter :: UnitLog=3          ! Logfile OutPut
integer(4), parameter :: UnitIn=4           ! user input file
integer(4), parameter :: StdIn=5            ! Restart File
integer(4), parameter :: StdOut=6           ! StdOut 
integer(4), parameter :: UnitOut=7          ! slpf_xxxxxx.out (text output)
integer(4), parameter :: UnitRestart=8      ! restart_xxxxxx.bin (binary output)
integer(4), parameter :: LINE_LENGTH=120    ! input file line length
integer(4), parameter :: LEN_TITLE=53 ! title length
integer(4), parameter :: FILENAME_LENGTH=120    ! length of char str for file names
logical(4):: PauseOn=.true.     !#unix uncomment  !#linux uncomment  !#win comment   
character*40 :: Input_fname = 'slpf.inp' !default
character*40 :: Restart_fname = 'restart_000001.bin' !default

character(LEN_TITLE) :: version
integer(4) :: FluidID = 0
integer(4) :: NumFluidProcs = 1
logical(4), private :: ErrorFileFlag=.false.

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Utility_GetFname - Grap input filename from restart file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Utility_GetFname(Input_Fname)

   character*40, intent(out) :: Input_Fname
   logical(4) :: UnitInOpen 

   open(UnitRestart,file=trim(Restart_Fname),form='unformatted')
   read(UnitRestart) Input_Fname
   close(UnitRestart)

   inquire(Unit=UnitIn,opened=UnitInOpen)
   if(UnitInOpen)then
     write(StdOut,*)
     write(StdOut,*)'+++++ ERROR MESSAGE +++++'   
     write(StdOut,*)'Cannot find restart input file ',Input_Fname
     call Utility_PauseCode
   endif

end subroutine utility_GetFname

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Utility_Load - Load Input File
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Utility_Load(VersionIn,NumFluidProcsIn,FluidIDIn)
!  Utility_Load - Load private variables

!  --- arguements ---
   character(LEN_TITLE),intent(in):: VersionIn
   integer(4), intent(in) :: FluidIDin
   integer(4), intent(in) :: NumFluidProcsIn

   Version=VersionIn
   FluidID=FluidIdin
   NumFluidProcs=NumFluidProcsIn
   ErrorFileFlag=.true. 

   return   
end subroutine Utility_Load

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
   end do

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

!++++++++++++++++++++++++++++++++++ Utility_ErrorMsg +++++++++++++++++++++++++++
subroutine Utility_InitErrorMsg(Msg,ID,NumFluidProcs)
!  write error message and terminate

!  ----- arguments -----
   character(*), intent(in):: Msg  ! error message
   integer(4), intent(in) :: ID
   integer(4), intent(in) :: NumFluidProcs

!  ----- local data -----
   character(20) :: filename

   if(NumFluidProcs==1)then
      FileName='StartError.out'
   else
      FileName='StartError'//'_xxx.out'
      call Utility_SetCycle(FileName,12,ID,3)
   endif
   open(UnitOut,file=filename)

   write(UnitOut,*)
   write(UnitOut,*)'+++++ ERROR MESSAGE +++++'     
   write(UnitOut,*)Msg

   write(StdOut,*)
   write(StdOut,*)'+++++ ERROR MESSAGE +++++'  
   write(StdOut,*)Msg

   call WriteErrorFile(Msg)

   call Utility_PauseCode

   stop 
end subroutine Utility_InitErrorMsg

!++++++++++++++++++++++++++++++++++ WriteFileFlag +++++++++++++++++++++++++++
subroutine WriteErrorFile(Msg)
!  write error message and terminate

!  ----- arguments -----
   character(*), intent(in):: Msg  ! error message

   character(13) :: filename

   if(ErrorFileFlag)then
      if(NumFluidProcs>1)then
         filename='error_xxx.log'
         call Utility_SetCycle(filename,7,FluidID,3)
      else
         filename='error.log'
      endif
         
      open(UnitStop,file=trim(filename))
      write(UnitStop,'(a)') trim(Version) 
      write(UnitStop,'(a,i8)')'Processor = ',FluidID
      write(UnitStop,*)
      write(UnitStop,*) '+++++ ERROR MESSAGE +++++'
      write(UnitStop,*) Msg
      write(UnitStop,*) 
      write(UnitStop,*) 
      write(UnitStop,*) 
      close(UnitStop)
   endif

   return
end subroutine WriteErrorFile

!++++++++++++++++++++++++++++++++++ Utility_ErrorMsgChar +++++++++++++++++++++++++++
subroutine Utility_ErrorMsgChar(Msg,Name)     
!  write error message, String and terminate
implicit none

!  ----- arguments -----
   character(*), intent(in):: Msg     
   character(*), intent(in) :: Name
   logical(4) :: UnitOutOpen   
   logical(4) :: UnitLogOpen

   inquire(Unit=UnitLog,opened=UnitLogOpen)
   if(UnitLogOpen)write(UnitLog,*)Msg,Name

   inquire(Unit=UnitOut,opened=UnitOutOpen)
   if(UnitOutOpen)then
      write(UnitOut,*)
      write(UnitOut,*)'+++++ ERROR MESSAGE +++++'   
   write(UnitOut,*)Msg,Name
   endif


   write(StdOut,*)
   write(StdOut,*)'+++++ ERROR MESSAGE +++++'
   write(StdOut,*)Msg,Name

   call WriteErrorFile(Msg//Name)

   call Utility_PauseCode

   stop 
end subroutine Utility_ErrorMsgChar

!++++++++++++++++++++++++++++++++++ Utility_ErrorMsg +++++++++++++++++++++++++++
subroutine Utility_ErrorMsg(Msg)
!  write error message and terminate

!  ----- arguments -----
   character(*), intent(in):: Msg  ! error message
   logical(4) :: UnitOutOpen   
   logical(4) :: UnitLogOpen

   inquire(Unit=UnitLog,opened=UnitLogOpen)
   if(UnitLogOpen)write(UnitLog,*)Msg

   inquire(Unit=UnitOut,opened=UnitOutOpen)
   if(UnitOutOpen)then
      write(UnitOut,*)
      write(UnitOut,*)'+++++ ERROR MESSAGE +++++'   
      write(UnitOut,*)Msg
   endif


   write(StdOut,*)
   write(StdOut,*)'+++++ ERROR MESSAGE +++++'
   write(StdOut,*)Msg

   call WriteErrorFile(Msg)

   call Utility_PauseCode

   stop 
end subroutine Utility_ErrorMsg

!+++++++++++++++++++++++++++++++++ Utility_ReplaceChar +++++++++++++++++++++++
subroutine Utility_PauseCode

if(PauseOn)then
   write(StdOut,*)   
   write(StdOut,*)'Hit Any Key to Continue'
   read(StdIn,*)
endif

return
end subroutine Utility_PauseCode

!++++++++++++++++++++++++++++++++++ Utility_ErrorMsgChar +++++++++++++++++++++++++++
subroutine Utility_InputError(ios,routine,section,Line)
 
!  write error message, String and terminate
implicit none

!  ----- arguments -----
   integer(4), intent(in) :: ios
   character(*), intent(in):: routine     
   character(*), intent(in) :: section
   character(*), intent(in) :: Line
   logical(4) :: UnitOutOpen   
   logical(4) :: UnitLogOpen

   inquire(Unit=UnitLog,opened=UnitLogOpen)
   if(UnitLogOpen)then
      write(UnitLog,'(a,a,a,a,a,i5)')'Read error in routine ',routine,' for ',section,'  ios = ',ios
   endif

   inquire(Unit=UnitOut,opened=UnitOutOpen)
   if(UnitOutOpen)then   
      write(UnitOut,*)
      write(UnitOut,*)'+++++ ERROR MESSAGE +++++' 
      write(UnitOut,'(/,a,a,a,a,a)')'Read error ',section,'  (routine ',routine,')'
      write(UnitOut,*)'Bad Line:'
      write(UnitOut,'(1a80)')Line
   endif

   write(StdOut,*)
   write(StdOut,*)'+++++ ERROR MESSAGE +++++'    
   write(StdOut,'(/,a,a,a,a,a)')'Read error ',section,'  (routine ',routine,')'
   write(StdOut,*)'Bad Line:'
   write(StdOut,'(1a80)')Line

   call Utility_PauseCode

   stop 
end subroutine Utility_InputError

end module


