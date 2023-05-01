!  grad.f90 
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
!  PROGRAM: Grad(Module)
!
!  PURPOSE: Contains subroutines to generate grad matrix
! 
!  SUBROUTINES: Grad_Main - Main Startup
!               Grad_Potent - Calculates inital channel potential
!               Grad_Allocate - allocate device matrices
!
!****************************************************************************
module Grad

!Use global variables from module
use Startup
use Debug
use Errors

implicit none

!! Define Global Variables !!
real(kind=8), allocatable, dimension(:,:) :: dz

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Device_Main - Main Startup
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Grad_Main

   call Grad_Allocate ! allocate device matrix
   call Grad_Dz ! Input Channel Potential
   if(rank .eq. 0 .and. debug_l > 5)call Grad_Plot

end subroutine Grad_Main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Device_Potent - Calculates inital channel potential
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Grad_Dz

    integer(kind=4) :: r, n, k

    do r = 1, Np
      dz(r,r) = -2*1/Gdx(r)
      if(r .gt. 1) dz(r,r-1) = 1*1/Gdx(r)
      if(r .lt. Np) dz(r,r+1) = 1*1/Gdx(r)
    end do

end subroutine Grad_Dz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Gradplot - Plot Grad
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Grad_Plot

    integer(kind=4):: err, r
    character*100 :: var

    open(unit=13, file='Dz.dat', status='new', action='write', iostat=err)
    write(var,*) 'Dz.dat'
    call Errors_Fileopen(err,var)

    do r=1, Np
      write(13,*) dz(r,:) ! Write Diagonals of Matrix to File for plotting
    end do
    close(unit=13)

end subroutine Grad_Plot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Device_allocate - allocate Device Matrices
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Grad_Allocate

    integer(kind=4) :: err
    character*100 :: var

    allocate(dz(Np,Np), stat=err) ! Allocate Uo
    write(var,*) 'dz(Np)'
    call Errors_Allocate(err,var) !Check Errors
   
    dz(:,:) = 0.0 ! Setting initial potential in channel to be zero 

end subroutine Grad_Allocate

end module
