!  restart.f90
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
!  PURPOSE: Contains subroutine to output restart data to bin file
!
!  SUBROUTINES: 
!               
!               
!                
!
!****************************************************************************
!
module Restart

use Startup
implicit none

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Restart_Write - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Restart_Write(Its,Itot,V,Imax,Vmax,Vneg,Vpos,Ineg,Ipos,Imin,c1,c2,c3,c4,Gn,Gp,Rh)

   integer(kind=4) :: r, j, m
   integer(kind=4), intent(in) :: c1,c2,c3,c4
   real(kind=8), intent(in) :: Its,Itot,V,Imax,Vmax,Vneg,Vpos,Ineg,Ipos,Imin
   real(kind=8), intent(in), dimension(:) :: Rh
   real(kind=8), intent(in), dimension(:,:) :: Gn,Gp

   character*40 :: Fname

   if(debug_l>5)write(6,*)'+++ Restart; Writing Restart Data +++'

   Fname = Utility_MakeFileName('restart',NumFluidProcs,c4,0,'bin')
   open(UnitRestart,file=Fname,form='unformatted')

   !write header file
   write(UnitRestart) Input_Fname

   write(UnitRestart) NE,Np,c1,c2,c3,c4

   write(UnitRestart) Its,Itot,V,Imax,Vmax,Vneg,Vpos,Ineg,Ipos,Imin

   do r = 1, NE
     write(UnitRestart) Rh(r) 
   end do

   do r = 1, NE
     do j = 1, Np
         write(UnitRestart) Gn(j,r),Gp(j,r)
     end do
   end do


   !do r = 1, NE
   !  write(UnitRestart) f_e(r) 
   !end do

   close(UnitRestart)

end subroutine Restart_Write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Restart_Write - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Restart_Read(Its,Itot,V,Imax,Vmax,Vneg,Vpos,Ineg,Ipos,Imin,c1,c2,c3,c4,Gn,Gp,Rh)

   integer(kind=4) :: r, j, m
   integer(kind=4), intent(out) :: c1,c2,c3,c4
   real(kind=8), intent(out) :: Its,Itot,V,Imax,Vmax,Vneg,Vpos,Ineg,Ipos,Imin
   real(kind=8), intent(out), dimension(:) :: Rh
   real(kind=8), intent(out), dimension(:,:) :: Gn,Gp

   character*40 :: Fname

   if(debug_l>5)write(6,*)'+++ Restart; Reading Restart Data +++'

   open(UnitRestart,file=Restart_Fname,form='unformatted')

   !write header file
   read(UnitRestart) Input_Fname

   read(UnitRestart) NE,Np,c1,c2,c3,c4

   read(UnitRestart) Its,Itot,V,Imax,Vmax,Vneg,Vpos,Ineg,Ipos,Imin

   do r = 1, NE
     read(UnitRestart) Rh(r) 
   end do

   do r = 1, NE
     do j = 1, Np
         read(UnitRestart) Gn(j,r),Gp(j,r)
     end do
   end do

   !do r = 1, NE
   !  read(UnitRestart) f_e(r) 
   !end do

   close(UnitRestart)

end subroutine Restart_Read
end module
