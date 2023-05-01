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

implicit none

integer(kind=4) :: debug_l = 10 !Default debug level is 10
character*16 :: f_name
character*40 :: version = 'sltc(serial) ver.1.2.2 (4/19/2010)'


contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Debug_Hermi - Check Hermicity of Hamiltonian
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Debug_Hermi(H,n)

  integer(kind=4) :: i, j, nh, err
  character*100 :: var
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
       if(Hh(i,j) .ne. 0.0) nh = nh + 1
    end do
  end do

  write(var,*) '** Checking Hamiltonian Hermicity'
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
! Debug_Output - Print Debug Data to file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Debug_Output(var)

  integer(kind=4) :: err
  character*100, intent(in) :: var

  open(unit=10, file=f_name, status='old', access='append', iostat=err)

  write(10,*) var

  close(10)

end subroutine Debug_Output

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

  call date_and_time(TIME=time) ! Call current time for filename
  call date_and_time(VALUES=values) ! To timemark inside file
  call itime(timeArray)
  rnum = rand(timeArray(3)) !time
  rnum = int(rand(0)*9)+1 !need a random number because on the cluster jobs start simulateously
  !                       and two jobs write into same file then it is hard to parse.
  write(f_name,'(A5,A6,I1,A4)') 'sltc_',time,rnum,'.out' ! Construct filename

  inquire(file=f_name, exist=exist1, direct=dir_acc, readwrite=rw_stat) !Check file status

  if(exist1 .eqv. .False.) then
     open(unit=10, file=f_name, status='new', action='write', iostat=err)
  else
     call date_and_time(TIME=time) ! Call current time for filename
     call date_and_time(VALUES=values) ! To timemark inside file
     rnum = int(rand(0)*9)+1 !need a random number because on the cluster jobs start simulateously
     !                       and two jobs write into same file then it is hard to parse.
     write(f_name,'(A5,A6,I1,A4)') 'sltc_',time,rnum,'.out' ! Construct filename
     open(unit=10, file=f_name, status='new', action='write', iostat=err)
  end if

  write(10,*) 'Superlattice Quantum Transport Calculator - Thermoelectric Electrical Properties'
  write(10,*) version
  write(10,'(I2,A1,I2,A1,I4,2X,I2,A1,I2,A1,I2)') values(2),'/',values(3),'/',values(1),values(5),':',values(6),':',values(7)
  write(10,*) '--------------------------------------------------'
  write(10,*)

  do n=1,mn
     write(10,*) var(n) !Write to file
  end do

  write(10,*) '--------------------------------------------------'
  write(10,*)

  close(10)

end subroutine Debug_ParmOut

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Debug_Stdout - Print Debug Data to screen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Debug_Stdout(var)

  character*100, intent(in) :: var

  if(debug_l > 2) then
     write(6,*) var !write to screen
  end if

end subroutine Debug_Stdout

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Debug_Comarg - Print Command Line Arguments
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Debug_Comarg(vin,Lb,Lw,Lc,Nl,Nd_Ge)

  integer(kind=4), intent(in) :: vin
  integer(kind=4), intent(in) :: Nl
  real(kind=8), intent(in) :: Nd_Ge,Lb,Lw,Lc
 

  if(vin > 8) then
  write(*,'(1X,A)') 'Command Line Arguments'
  write(*,'(1X,A)') '----------------------'
  write(*,'(1X,A,I3)') 'Debug Level - ',debug_l
  write(*,'(1X,A,E15.6)') 'Barrier Length - ',Lb
  write(*,'(1X,A,E15.6)') 'Well Length - ',Lw
  write(*,'(1X,A,E15.6)') 'Contact Length, - ',Lc
  write(*,'(1X,A,I3)') 'Number of Layers - ',Nl
  write(*,'(1X,A,E15.6)') 'Doping Concentration - ',Nd_Ge
  write(*,*)
  end if

end subroutine Debug_Comarg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Debug_Title - Print Title
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Debug_Title

  write(6,*) ' '
  write(6,*) '********************************************************************************'
  write(6,*) '*                                                                              *'
  write(6,*) '*              Superlattice (SLTC) Quantum Transport Calculator                *'
  write(6,*) '*                  Thermoelectric Thermal Properties                           *'
  write(6,*) '*                                                                              *'
  write(6,*) '*                         Musho, T.D., Walker, D.G.                            *'
  write(6,*) '*               Vanderbilt University - Thermal Engineering Lab                *'
  write(6,*) '*                                                                              *'
  write(6,*) '********************************************************************************'
  write(6,*) '*                   ',version,'                   *'
  write(6,*) '********************************************************************************'
  write(6,*) ' '

end subroutine Debug_Title
 

end module


