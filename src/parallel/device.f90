!  device.f90 
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
!  PROGRAM: Device(Module)
!
!  PURPOSE: Contains subroutines to construct device layout
! 
!  SUBROUTINES: Device_Main - Main Startup
!               Device_Potent - Calculates inital channel potential
!               Device_Allocate - allocate device matrices
!
!****************************************************************************
module Device

!Use global variables from module
use Startup
use Debug
use Errors

implicit none

!! Define Global Variables !!
real(kind=8), allocatable, dimension(:) :: Uo, Ubdy
real(kind=8), allocatable, dimension(:,:) :: U

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Device_Main - Main Startup
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Device_Main

   call Device_Allocate ! allocate device matrix
   call Device_Potent ! Input Channel Potential
   call Device_Init_bdy
   if(debug_l>5 .and. rank .eq. 0) call Device_Pot_Plot

end subroutine Device_Main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Device_Potent - Calculates inital channel potential
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Device_Potent

    integer(kind=4) :: r, n, k

    !!! Permittivity

    k=0
    do r = 1, Nl
       if(r .eq. 1) then
          Uo(1) = Gdx(1)/epsilrS
          k = 0
          if(Nec .ne. 0) then
            do n = 2, Nec
              Uo(n) = Gdx(n)/epsilrS ! silicon permittivity
              k=n
            end do
            Uo(k+1) = Gdx(k+1)/epsilrS
          end if
       end if

       do n = k+2, k+Neb
          Uo(n) = Gdx(n)/epsilrS ! germanium permittivity
          k=n
       end do
       Uo(k+1) = 0.5*(Gdx(k)/epsilrS + Gdx(k)/epsilrG) ! interface permittivity

       do n = k+2, k+New 
          Uo(n) = Gdx(n)/(epsilrG) ! germanium permittivity
          k=n
       end do
       Uo(k+1) = 0.5*(Gdx(k)/epsilrS + Gdx(k)/epsilrG) ! interface permittivity

       if(r .eq. Nl) then
          do n = k+2, k+Neb
             Uo(n) = Gdx(n)/epsilrS ! silicon permittivity
             k=n
          end do
          if(Nec .ne. 0) then
            Uo(k+1) = Gdx(k)/epsilrS
            do n = k+2, k+Nec
              Uo(n) = Gdx(n)/epsilrS ! silicon permittivity
              k=n
            end do
          end if
          Uo(Np) = Gdx(Np)/epsilrS
       end if
    end do

    !Mod Hamiltonian for no leads
    if(Nec .eq. 0) then
       Uo(1) = Gdx(1)/epsilrS
       Uo(Np) = Gdx(Np)/epsilrS
    end if

    !!! Channel Potential 
    Uo = Uo * q/(2.0*pi*epsil0) * log(2.0)   ! Units of eV 

end subroutine Device_Potent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Device_Init_bdy - Set conduction band edge of body
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Device_Init_bdy

    integer(kind=4) :: r, n, k

    k=0
    do r = 1, Nl
       if(r .eq. 1) then
          Ubdy(1) = Ecs
          k = 0
          if(Nec .ne. 0) then
            do n = 2, Nec
              Ubdy(n) = Ecs 
              k=n
            end do
            Ubdy(k+1) = Ecs
          end if
       end if

       do n = k+2, k+Neb
          Ubdy(n) = Ecs 
          k=n
       end do
       Ubdy(k+1) = 0.5*(Ecs + Ecg) 

       do n = k+2, k+New 
          Ubdy(n) = Ecg 
          k=n
       end do
       Ubdy(k+1) = 0.5*(Ecs + Ecg)

       if(r .eq. Nl) then
          do n = k+2, k+Neb
             Ubdy(n) = Ecs 
             k=n
          end do
          if(Nec .ne. 0) then
            Ubdy(k+1) = Ecs
            do n = k+2, k+Nec
              Ubdy(n) = Ecs
              k=n
            end do
          end if
          Ubdy(Np) = Ecs
       end if
    end do

    !Mod Hamiltonian for no leads
    if(Nec .eq. 0) then
       Ubdy(1) = Ecs
       Ubdy(Np) = Ecs
    end if

end subroutine Device_Init_bdy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Device_Potplot - Plot potential
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Device_Pot_Plot

    integer(kind=4):: err, r
    real(kind=4) :: l
    character*120 :: var

       open(unit=13, file='Uo.dat', status='new', action='write', iostat=err)
       write(var,*) 'Uo.dat'
       call Errors_Fileopen(err,var)
       open(unit=14, file='Ua.dat', status='new', action='write', iostat=err)
       write(var,*) 'Ua.dat'
       call Errors_Fileopen(err,var)

       l = 0.0
       do r=1, Np
          write(13,*) l, Uo(r) ! Write Diagonals of Matrix to File for plotting
          write(14,*) l, Ubdy(r)
          l = Gdx(r) + l
       end do
       close(unit=13)
       close(unit=14)

end subroutine Device_Pot_Plot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Device_allocate - allocate Device Matrices
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Device_Allocate

    integer(kind=4) :: err
    character*120 :: var

    allocate(Uo(Np), Ubdy(Np), stat=err) ! Allocate Uo
    write(var,*) 'Uo(Np)'
    call Errors_Allocate(err,var) !Check Errors
    
    allocate(U(Np,Np), stat=err) ! Allocate U
    write(var,*) 'U(Np,Np)'
    call Errors_Allocate(err,var) !Check Errors
   
    U(:,:) = 0.0 ! Setting initial potential in channel to be zero 
    Uo(:) = 0.0; Ubdy(:) = 0.0 ! Specifying size of permittivity matrix

end subroutine Device_Allocate

end module
