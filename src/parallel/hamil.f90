!  hamil.f90 
!
!  Authors: G. Walker - Vanderbilt University
!           T. Musho - Vanderbilt University         
!
!  Created: 1/30/09
!  Last Modified: 1/30/09
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
!  PROGRAM: Hamil(Module)
!
!  PURPOSE: Contains subroutine to construct hamiltonian
!
!  SUBROUTINES: Hamil_Main - Main Startup, subroutine calls
!               Hamil_Create - Create Hamiltonian Matrix
!               Hamil_Allocate - Allocate hamiltonian, zero
!               Hamil_Plot - Plot Hamiltonian, full(raw), main diag 
!
!****************************************************************************
!
module Hamil

!Use global variables from module
use Startup
use Device
use Errors

implicit none

!! Define Global Variables !!
real(kind=8), allocatable, dimension(:,:) :: H

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Hamil_Main - Main Startup, subroutine calls
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Hamil_Main

   call Hamil_Allocate ! Allocate device matrix
   call Hamil_Create ! Input values into device matrix
   if(debug_l>5 .and. rank .eq. 0) call Hamil_Plot ! Plot hamiltonian diag and raw
   if(rank .eq. 0)call Debug_Hermi(H,Np) ! Check Hermicity

end subroutine Hamil_Main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Hamil_Create - Create Hamiltonian Matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Hamil_Create

    integer(kind=4) :: r, n, k
    !(Ec_con,Ec_bar,Ec_wel,t_c,t_b,t_w,Nd,Nb,Nw,Gd)
    !integer(kind=8), intent(in) :: Nd, Nb, Nw
    !real(kind=8), intent(in) :: Ec_con, Ec_bar, Ec_wel, t_c, t_b, t_w
    !real(kind=8), intent(in), dimension(Nd) :: Gd

    !!! diagonal Hamiltonian terms
    k=0
    do r = 1, Nl
       if(r .eq. 1) then  
         H(1,1) = Ecs + 2*ts*(1/Gdx(1)**2)           
         k = 0
         if(Nec .ne. 0) then
           do n = 2, Nec
             H(n,n) = Ecs + ts*(1/Gdx(n)**2+1/Gdx(n-1)**2)
             k=n
           end do
           H(k+1,k+1) = Ecs + ts*(1/Gdx(k+1)**2+1/Gdx(k)**2)
         end if
       end if

       do n = k+2, k+Neb
          H(n,n) = Ecs + ts*(1/Gdx(n)**2+1/Gdx(n-1)**2)
          k=n
       end do
       H(k+1,k+1) = Ecb + (tg/Gdx(k+1)**2+ts/Gdx(k)**2)

       do n = k+2, k+New
          H(n,n) = Ecg + (tg)*(1/Gdx(n)**2+1/Gdx(n-1)**2)
          k=n
       end do
       H(k+1,k+1) = Ecb + (ts/Gdx(k+1)**2+tg/Gdx(k)**2)

       if(r .eq. Nl) then
          do n = k+2, k+Neb
             H(n,n) = Ecs + (ts)*(1/Gdx(n)**2+1/Gdx(n-1)**2)
             k=n
          end do
          if(Nec .ne. 0) then
            H(k+1,k+1) = Ecs + (ts)*(1/Gdx(k+1)**2+1/Gdx(k)**2)
            do n = k+2, k+Nec
              H(n,n) = Ecs + (ts)*(1/Gdx(n)**2+1/Gdx(n-1)**2)
              k=n
            end do
          end if
          H(Np,Np) = Ecs + 2*ts*(1/Gdx(Np)**2)
       end if
    end do

    !!! off-diagonal Hamiltonian terms
    k=0
    do r = 1, Nl
       if(r .eq. 1) then
         H(2,1) = -ts*(1/Gdx(n)**2)
         H(1,2) = -ts*(1/Gdx(n)**2)  
         if(Nec .ne. 0) then
           do n = 2, Nec-1
             H(n,n+1) = -ts*(1/Gdx(n)**2)
             H(n+1,n) = -ts*(1/Gdx(n)**2)
             k=n
           end do
             H(k+1,k+2) = -ts*(1/Gdx(k+1)**2)
             H(k+2,k+1) = -ts*(1/Gdx(k+1)**2)
         end if
        end if

       do n = k+1, k+Neb
         H(n,n+1) = -ts*(1/Gdx(n)**2)
         H(n+1,n) = -ts*(1/Gdx(n)**2)
         k=n
       end do

       do n = k+1, k+New
         H(n,n+1) = -tg*(1/Gdx(n)**2)
         H(n+1,n) = -tg*(1/Gdx(n)**2)
         k=n
       end do

       if(r .eq. Nl) then
          do n = k+1, k+Neb
             H(n,n+1) = -ts*(1/Gdx(n)**2)
             H(n+1,n) = -ts*(1/Gdx(n)**2)
             k=n
          end do
          if(Nec .ne. 0) then
            H(k+1,k+2) = -ts*(1/Gdx(k+1)**2)
            H(k+2,k+1) = -ts*(1/Gdx(k+1)**2)
            do n = k+2, k+Nec+1
              H(n,n+1) = -ts*(1/Gdx(n)**2)
              H(n+1,n) = -ts*(1/Gdx(n)**2)
              k=n
            end do
          end if
       end if
    end do

end subroutine Hamil_Create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Hamil_Allocate - Allocate Hamiltonian Matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Hamil_Allocate

    integer(kind=4) :: err
    character*120 :: var

    !! Hamiltonian Matrix, Units of eV !! 
    allocate(H(Np,Np), stat=err) ! Allocate Matrix
    write(var,*) 'H(Np,np)'
    call Errors_Allocate(err,var) !Check Errors

    H(:,:) = 0.0 !Define Hamiltonian Matrix to Zero

end subroutine Hamil_Allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Hamil_Allocate - Allocate Hamiltonian Matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Hamil_Plot

    integer(kind=4) :: err, r
    real(kind=8) :: l
    character*120 :: var

    open(unit=14, file='Hraw.dat', status='new', action='write', iostat=err)
    write(var,*) 'Hraw.dat'
    call Errors_Fileopen(err,var)
    open(unit=15, file='H.dat', status='new', action='write', iostat=err)
    write(var,*) 'H.dat'
    call Errors_Fileopen(err,var)

    l = 0 
    do r=1, Np
      if(r .eq. 1) then
        write(15,*) l,Gdx(r),H(r,r),0.0,-H(r,r+1) ! Write Diagonals of Matrix to File for plotting
      else if(r .eq. Np) then
        write(15,*) l,Gdx(r),H(r,r),-H(r,r-1),0.0 ! Write Diagonals of Matrix to File for plotting
      else
        write(15,*) l,Gdx(r),H(r,r),-H(r,r-1),-H(r,r+1) ! Write Diagonals of Matrix to File for plotting
      end if
      l = l + Gdx(r)
    end do

    do r=1, Np
      write(14,*) H(r,:) ! Write Raw Hamiltonian to File
    end do    

    close(unit=14)
    close(unit=15)

end subroutine Hamil_Plot

end module
