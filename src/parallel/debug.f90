!  debug.f90 
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
!  PROGRAM: Debug(Module)
!
!  PURPOSE: Contains debug routines
!
!  SUBROUTINES: Init_Main - Main Startup, subroutine calls
!               Init_Allocate - Allocates arrays for NEGF
!               Init_Deallocate - Deallocates arrays 
!
!****************************************************************************
!
module Debug

use Errors
use Utility

implicit none

integer(kind=4) :: debug_l = 10 !Default debug level is 10
character*19 :: f_name


contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Debug_Hermi - Check Hermicity of Hamiltonian
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Debug_Hermi(H,n)

  integer(kind=4) :: i, j, nh, err
  character*120 :: var
  real(kind=8), allocatable, dimension(:,:) :: Hh

  integer(kind=4), intent(in) :: n
  real(kind=8), dimension(:,:), intent(in) :: H

  allocate(Hh(n,n), stat=err) ! Allocate Matrix
  write(var,*) 'Hh(Np,Np)'
  call Errors_Allocate(err,var)

  Hh = H - transpose(H) !Check Hermiticity

  nh = 0 !set counter to zero

  do j = 1, n ! Spectral Matrix horizontal && vertical centerlines
    do i = 1, n
       if(Hh(i,j) .ne. 0) nh = nh + 1
    end do
  end do

  write(var,*) ' +++ Checking Hamiltonian Hermicity  +++'
  call Debug_Stdout(var)
  call Debug_Output(var)

  if(nh .gt. 0) then 
    write(var,*) '!! Warning Hamiltonian Non-Hermitian !!'
    call Debug_Stdout(var)
    call Debug_Output(var)
  end if

  deallocate(Hh)

end subroutine Debug_Hermi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Debug_DPCHIM - Check errors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Debug_DPCHIM(err_code)

  integer(kind=4) :: i, j, nh, err
  character*120 :: var

  integer(kind=4), intent(in) :: err_code

  if(err_code .gt. 0) then 
    write(var,*) '!! Warning: DPCHIM - IERR switches in the direction of monotonicity were detected. !!'
    !call Debug_Stdout(var)
    call Debug_Output(var)
  end if

  if(err_code .eq. -1) then 
    write(var,*) '!! Error: DPCHIM - N.LT.2 . !!'
    call Debug_Stdout(var)
    call Debug_Output(var)
  end if

  if(err_code .eq. -2) then 
    write(var,*) '!! Error: DPCHIM - INCFD.LT.1 . !!'
    call Debug_Stdout(var)
    call Debug_Output(var)
  end if

  if(err_code .eq. -3) then 
    write(var,*) '!! Error: DPCHIM - X-array is not strictly increasing. !!'
    call Debug_Stdout(var)
    call Debug_Output(var)
  end if

end subroutine Debug_DPCHIM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Debug_DPCHIM - Check errors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Debug_DPCHFE(err_code)

  integer(kind=4) :: i, j, nh, err
  character*120 :: var

  integer(kind=4), intent(in) :: err_code

  if(err_code .gt. 0) then 
    !write(var,*) '!! Warning: DPCHFE - means that extrapolation was performed at IERR points !!'
    !call Debug_Stdout(var)
    call Debug_Output(var)
  end if

  if(err_code .eq. -1) then 
    write(var,*) '!! Error: DPCHFE - N.LT.2 . !!'
    call Debug_Stdout(var)
    call Debug_Output(var)
  end if

  if(err_code .eq. -2) then 
    write(var,*) '!! Error: DPCHFE - INCFD.LT.1 . !!'
    call Debug_Stdout(var)
    call Debug_Output(var)
  end if

  if(err_code .eq. -3) then 
    write(var,*) '!! Error: DPCHFE - NE.LT.1 . !!'
    call Debug_Stdout(var)
    call Debug_Output(var)
  end if

end subroutine Debug_DPCHFE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Debug_Output - Print Debug Data to file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Debug_Output(var)

  integer(kind=4) :: err
  character*120, intent(in) :: var

  !open(unit=UnitLog, status='old', access='append', iostat=err)

  write(UnitLog,*) var

  flush(UnitLog)
  !close(UnitLog)

end subroutine Debug_Output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Debug_Output - Print Debug Data to file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Debug_DebugOutput(var)

  integer(kind=4) :: err
  character*120, intent(in) :: var

  !open(unit=UnitLog, status='old', access='append', iostat=err)

  write(UnitDebug,*) var

  flush(UnitDebug)
  !close(UnitLog)

end subroutine Debug_DebugOutput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Debug_SetTime - Set time to filename
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Debug_SetTime(ftname)

  character*120, intent(out) :: ftname
  integer,dimension(8) :: values
  character*10 :: time

  call date_and_time(TIME=time) ! Call current time for filename
  call date_and_time(VALUES=values) ! To timemark inside file
  write(ftname,'(A5,A6,A4)') 'slpf_',time,'.out' ! Construct filename

end subroutine Debug_SetTime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Debug_SetTimeFloat - Set time to float
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Debug_SetTimeFloat(ftname)

  integer(kind=4), intent(out) :: ftname
  integer,dimension(8) :: values
  character*10 :: time

  call date_and_time(TIME=time) ! Call current time for filename
  call date_and_time(VALUES=values) ! To timemark inside file
  !write(*,*) 'values',values
  ftname=int(values(5)*3600+values(6)*60+values(5)+values(6)) ! Construct filename

end subroutine Debug_SetTimeFloat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Debug_Output - Print Debug Data to file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Debug_ParmOut(var,mn)

  logical :: exist1
  integer(kind=4) :: err,n,rnum
  integer(kind=4), intent(in) :: mn
  integer,dimension(8) :: values
  character*9 :: dir_acc, rw_stat
  character*80,dimension(80), intent(in) :: var
  character*10 :: time
  integer(kind=4) :: timeArray(3)

  !call date_and_time(TIME=time) ! Call current time for filename
  !call date_and_time(VALUES=values) ! To timemark inside file
  !call itime(timeArray)
  !rnum = rand(timeArray(3)) !time
  !rnum = int(rand(0)*9)+1 !need a random number because on the cluster jobs start simulateously
  !                       and two jobs write into same file then it is hard to parse.
  !write(f_name,'(A5,A6,I1,A4)') 'slpf_',time,rnum,'.out' ! Construct filename

  !inquire(file=f_name, exist=exist1, direct=dir_acc, readwrite=rw_stat) !Check file status

  !if(exist1 .eqv. .False.) then
  !   open(unit=10, file=f_name, status='new', action='write', iostat=err)
  !else
  !   call date_and_time(TIME=time) ! Call current time for filename
  !   call date_and_time(VALUES=values) ! To timemark inside file
  !   rnum = int(rand(0)*9)+1 !need a random number because on the cluster jobs start simulateously
     !                       and two jobs write into same file then it is hard to parse.
  !   write(f_name,'(A5,A6,I1,A4)') 'slpf_',time,rnum,'.out' ! Construct filename
  !   open(unit=10, file=f_name, status='new', action='write', iostat=err)
  !end if

  do n=1,mn
     write(UnitLog,*) var(n) !Write to file
  end do


  close(10)

end subroutine Debug_ParmOut

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Debug_Stdout - Print Debug Data to screen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Debug_Stdout(var)

  character*120, intent(in) :: var

  if(debug_l > 2) then
     write(StdOut,*) var !write to screen
  end if

end subroutine Debug_Stdout


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Debug_Title - Print Title
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Debug_Title(ver)

  character*54, intent(in) :: ver

  100 format (1X,A,1X)
  101 format (1X,A,23X,A,23X,A,1X)
   
  !                 |1         |10       |         |        |40       |50       |60        |70       |80       |90       |100
  write(StdOut,100) ' '
  write(StdOut,100) '*****************************************************************************************************'
  write(StdOut,100) '*                                                                                                   *'
  write(StdOut,100) '*                         Superlattice (SLPF) Quantum Transport Calculator                          *'
  write(StdOut,100) '*                               Thermoelectric Electrical Properties                                *'
  write(StdOut,100) '*                                                                                                   *'
  write(StdOut,100) '*                                    Musho, T.D., Walker, D.G.                                      *'
  write(StdOut,100) '*                          Vanderbilt University - Thermal Engineering Lab                          *'
  write(StdOut,100) '*                                                                                                   *'
  write(StdOut,100) '*                                                                                                   *'
  write(StdOut,100) '*	      Copyright 2012 Terence Musho Licensed under the                                        *'
  write(StdOut,100) '*	      Educational Community License, Version 2.0 (the "License"); you may                    *'
  write(StdOut,100) '*	      not use this program except in compliance with the License. You may                    *'
  write(StdOut,100) '*	      obtain a copy of the License at                                                        *'
  write(StdOut,100) '*                                                                                                   *'
  write(StdOut,100) '*            http://www.osedu.org/licenses/ECL-2.0                                                  *'
  write(StdOut,100) '*                                                                                                   *'
  write(StdOut,100) '*	      Unless required by applicable law or agreed to in writing,                             *'
  write(StdOut,100) '*	      software distributed under the License is distributed on an "AS IS"                    *'
  write(StdOut,100) '*	      BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express                    *'
  write(StdOut,100) '*	      or implied. See the License for the specific language governing                        *'
  write(StdOut,100) '*	      permissions and limitations under the License.                                         *'
  write(StdOut,100) '*                                                                                                   *'
  write(StdOut,100) '*                                                                                                   *'
  write(StdOut,100) '*                                                                                                   *'
  write(StdOut,100) '*****************************************************************************************************'
  write(StdOut,101) '*',ver,'*'
  write(StdOut,100) '*****************************************************************************************************'
  write(StdOut,100) ' '

end subroutine Debug_Title
 

end module


