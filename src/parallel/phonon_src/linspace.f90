!  linspace.f90 
!
!  Authors: G. Walker - Vanderbilt University
!           A. Bulusu - Vanderbilt University
!
!  Modifier: T. Musho
!
!  Program was originally written in MathLab(Octave/GNU) Converted to Fortran 90/95
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
!  PROGRAM: LINSPACE
!
!  PURPOSE:  Fills Arrays with Specified Gradient 
!            Same as Linspace command in MatLab
!
!****************************************************************************
!
subroutine linspace (Len,size,sval,fval) !Len = Incoming Array, size = Size of Incoming Array, sval = Starting Value of Gradient,  fval = Ending Value of Gradient
implicit none

integer(kind=4) :: n
integer(kind=4), intent(IN) :: size
real(kind=8), intent(IN) :: sval, fval
real(kind=8), intent(INOUT), dimension(size) :: Len

Len(1) = sval

do n = 2, size
   Len(n) = Len(n-1) + (fval-sval)/dfloat(size-1) ! Evenly Spaced Array Gradient
end do

return
end subroutine linspace
