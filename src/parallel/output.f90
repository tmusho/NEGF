!  output.f90
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
!  PURPOSE: Contains subroutine to output data to text file
!
!  SUBROUTINES: 
!               
!               
!                
!
!****************************************************************************
!
module Output

!Use global variables from module
use Startup
use Errors
use Debug
use Utility


implicit none

integer(kind=4), private :: Unitout = 10 !Fileout
integer(kind=4), private :: Unitout1 = 11 !Fileout
integer(kind=4), private :: Unitout2 = 12 !Fileout
integer(kind=4), private :: Unitout3 = 13 !Fileout

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Output_Green - Plot Green's Matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Output_Green(G)

    integer(kind=4) :: r, n, err
    character*120 :: var
    complex(kind=8), dimension(:,:), intent(in) :: G

    open (unit=Unitout, file='Graw.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'Graw.dat'
    call Errors_Fileopen(err,var)
    open (unit=Unitout1, file='G.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'Error Opening :: G*.dat'
    call Errors_Fileopen(err,var)
 
    100 FORMAT (100(ES12.5,1X))

    do r=1, Np
       write(Unitout1,100) (DBLE(G(n,n))) ! Write Diagonals of Matrix to file for plotting
    end do

    do r=1, Np
       write(Unitout,*) G(r,:) ! Write Raw Hamiltonian to file
    end do    

    close(unit=Unitout)
    close(unit=Unitout1)

end subroutine Output_Green

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Output_ldos - Plot ASCII output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Output_ldos(A1,Aall,Gn1, c)

    integer(kind=4) :: r, n, j, err
    real(kind=8) :: w, l
    character*120 :: var
    character*80 :: Fname

    integer(kind=4), intent(in) :: c
    real(kind=8), dimension(:,:), intent(in) :: A1, Gn1
    real(kind=8), dimension(:,:,:), intent(in) :: Aall

    Fname = Utility_MakeFileName('A',NumFluidProcs,c,FluidId,'dat')
    open (unit=Unitout, file=Fname, status='NEW', action='write', iostat=err)
    Fname = Utility_MakeFileName('Gn',NumFluidProcs,c,FluidId,'dat')
    open (unit=Unitout1, file=Fname, status='NEW', action='write', iostat=err)
    Fname = Utility_MakeFileName('Gp',NumFluidProcs,c,FluidId,'dat')
    open (unit=Unitout2, file=Fname, status='NEW', action='write', iostat=err)
    Fname = Utility_MakeFileName('Asub',NumFluidProcs,c,FluidId,'dat')
    open (unit=Unitout3, file=Fname, status='NEW', action='write', iostat=err)
    write(var,*) 'Error Opening :: A*.dat'
    call Errors_Fileopen(err,var)

    101 format (3(ES12.5,1X))
    104 format (1000(ES12.5,1X))

    !LDOS Contour Plot
    do r=1, NE
       w = 0; l = 0
       do n=1,Np
         write(Unitout,101) E(r),l,A1(n,r)/(2.0*pi) !LDOS Contour Plot
         write(Unitout3,104) E(r),l,(Aall(n,r,j)/(2.0*pi),j=1,Nb) !LDOS Contour Plot
         write(Unitout1,101) E(r),l,Gn1(n,r)/(2.0*pi) !Filled States Contour Plot
         write(Unitout2,101) E(r),l,(A1(n,r)-Gn1(n,r))/(2.0*pi) !Unfilled States Contour Plot
         l = Gdx(n) + l
       end do
       write(Unitout,*);write(Unitout1,*);write(Unitout2,*);write(Unitout3,*) (' ',j=1,Nb)
    end do

    close(Unitout);close(Unitout1);close(Unitout2);close(Unitout3)

end subroutine Output_ldos

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Output_Scin - Plot ASCII output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Output_Scin(Si1, Nb, c)

    integer(kind=4) :: r, n, j, err
    real(kind=8) :: w, l
    character*120 :: var
    character*80 :: Fname

    integer(kind=4), intent(in) :: c, Nb
    real(kind=8), dimension(:,:,:), intent(in) :: Si1

    Fname = Utility_MakeFileName('Scin',NumFluidProcs,c,FluidId,'dat')
    open (unit=Unitout, file=Fname, status='NEW', action='write', iostat=err)
    write(var,*) 'Error Opening :: Sc*.dat'
    call Errors_Fileopen(err,var)

    101 format (3(ES12.5,1X))
    104 format (1000(ES12.5,1X))

    !LDOS Contour Plot
    do r=1, NE
       w = 0; l = 0
       do n=1,Np
         write(Unitout,104) E(r),l,(abs(Si1(n,r,j))/(2.0*pi),j=1,4),abs(sum(Si1(n,r,:)))/(2.0*pi) !LDOS Contour Plot
         l = Gdx(n) + l
       end do
       write(Unitout,*) (' ',j=1,Nb)
    end do

    close(Unitout)

end subroutine Output_Scin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Output_Scin - Plot ASCII output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Output_tau(Si1, Nb, c)

    integer(kind=4) :: r, n, j, err
    real(kind=8) :: w, l
    character*120 :: var
    character*80 :: Fname

    integer(kind=4), intent(in) :: c, Nb
    real(kind=8), dimension(:,:,:), intent(in) :: Si1

    Fname = Utility_MakeFileName('Tau',NumFluidProcs,c,FluidId,'dat')
    open (unit=Unitout, file=Fname, status='NEW', action='write', iostat=err)
    write(var,*) 'Error Opening :: Tau*.dat'
    call Errors_Fileopen(err,var)

    101 format (3(ES12.5,1X))
    104 format (1000(ES12.5,1X))

    !LDOS Contour Plot
    do r=1, NE          !Take the average scattering rate spatially but use matthersons rule at single site
     write(Unitout,104) E(r),(abs(sum(Si1(:,r,j)/Np)/(2.0*pi)),j=1,4),abs(sum(Si1(:,r,:)/Np))/(2.0*pi) !LDOS Contour Plot
    end do

    close(Unitout)

end subroutine Output_tau

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Output_tr - Plot ASCII output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Output_tr(Tr, nsb, c)

    integer(kind=4) :: r, n, err
    real(kind=8) :: w, l
    character*120 :: var
    character*80 :: Fname

    integer(kind=4), intent(in) :: c, nsb
    real(kind=8), dimension(:,:), intent(in) :: Tr

    Fname = Utility_MakeFileName('Tr',NumFluidProcs,c,FluidId,'dat')
    open (unit=Unitout, file=Fname, status='NEW', action='write', iostat=err)
    write(var,*) 'Tr.dat'
    call Errors_Fileopen(err,var)

    102 format (1000(ES12.5,1X))

    !Transmission Plot
    do r=1, NE
      w = 0; l = 0
        write(Unitout,102) E(r),(Tr(r,n), n=1,Nb),sum(Tr(r,:))/(nsb-1)
    end do

    close(Unitout)

end subroutine Output_tr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Output_i - Plot ASCII output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Output_i(Id1, Id2, c)

    integer(kind=4) :: r, err
    character*120 :: var
    character*80 :: Fname

    integer(kind=4), intent(in) :: c
    real(kind=8), dimension(:,:), intent(in) :: Id1,Id2

    Fname = Utility_MakeFileName('I',NumFluidProcs,c,FluidId,'dat')
    open (unit=Unitout, file=Fname, status='NEW', action='write', iostat=err)
    write(var,*) 'I.dat'
    call Errors_Fileopen(err,var)

    101 format (100(ES12.5,1X))

    !Current Plot
    do r=1, NE
      write(Unitout,101) E(r),sum(Id1(r,:)),sum(Id2(r,:))
    end do

    close(Unitout)

end subroutine Output_i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Output_i - Plot ASCII output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Output_i_scat(Id1, Id2, Id3, c)

    integer(kind=4) :: r, err
    character*120 :: var
    character*80 :: Fname

    integer(kind=4), intent(in) :: c
    real(kind=8), dimension(:,:), intent(in) :: Id1,Id2,Id3

    Fname = Utility_MakeFileName('I',NumFluidProcs,c,FluidId,'dat')
    open (unit=Unitout, file=Fname, status='NEW', action='write', iostat=err)
    write(var,*) 'I.dat'
    call Errors_Fileopen(err,var)

    101 format (100(ES12.5,1X))

    !Current Plot
    do r=1, NE
      write(Unitout,101) E(r),sum(Id1(r,:)),sum(Id2(r,:)),sum(Id3(r,:))
    end do

    close(Unitout)

end subroutine Output_i_scat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Output_u - Plot ASCII output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Output_u(Uo, U1, Ubdy, Ua1, No, c)

    integer(kind=4) :: r, n, j, err
    real(kind=8) :: w, l, mx
    character*120 :: var
    character*80 :: Fname

    integer(kind=4), intent(in) :: c
    real(kind=8), dimension(:), intent(in) :: Uo, Ubdy
    real(kind=8), dimension(:,:), intent(in) :: U1, No, Ua1

    Fname = Utility_MakeFileName('U',NumFluidProcs,c,FluidId,'dat')
    open (unit=Unitout, file=Fname, status='NEW', action='write', iostat=err)
    write(var,*) 'Error Opening :: U*.dat'
    call Errors_Fileopen(err,var)

    103 format (10(ES12.5,1X))

    ! Plot U
    l = 0.0
    do n=1, Np
      mx = -Uo(n)*sum(No(n,:))
      write(Unitout,103) l,U1(n,1),U1(n,2),U1(n,3),U1(n,4),U1(n,5),sum(U1(n,:)),mx,Ua1(n,c),Ua1(n,c)+Ubdy(n)+mx+sum(U1(n,:))
      l = Gdx(n) + l
    end do 

    close(Unitout)

end subroutine Output_u

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Output_rho - Plot ASCII output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Output_rho(Rh1, No, c)

    integer(kind=4) :: r, n, j, err
    real(kind=8) :: w, l
    character*120 :: var
    character*80 :: Fname

    integer(kind=4), intent(in) :: c
    real(kind=8), dimension(:), intent(in) :: Rh1
    real(kind=8), dimension(:,:), intent(in) :: No

    Fname = Utility_MakeFileName('Rho',NumFluidProcs,c,FluidId,'dat')
    open (unit=Unitout, file=Fname, status='NEW', action='write', iostat=err)
    write(var,*) 'Error Opening :: Rho*.dat'
    call Errors_Fileopen(err,var)

    101 FORMAT (1000(ES12.5,1X))

    l = 0.0
    do n=1, Np
      write(Unitout,101) l,(Rh1(n))/an*1E-6,sum(No(n,:))/Gdx(n)*1E-6 !Electron Density [1/m^2]
      l = Gdx(n) + l
    end do 

    close(Unitout)

end subroutine Output_rho 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Output_iv - Plot ASCII output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Output_iv(IV,c)

    integer(kind=4) :: n, err
    character*120 :: var

    integer(kind=4), intent(in) :: c
    real(kind=8), dimension(:,:), intent(in) :: IV

    open (unit=Unitout, file='IV.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'IV.dat'
    call Errors_Fileopen(err,var)

    !IV Plot
    do n = 1, NV
      write(Unitout,*) IV(n,1), IV(n,2)
    end do 

    close(Unitout)

end subroutine Output_iv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Output_mob - Plot ASCII output, Mobility
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Output_mob(MV,c)

    integer(kind=4) :: n, err
    character*120 :: var

    integer(kind=4), intent(in) :: c
    real(kind=8), dimension(:,:), intent(in) :: MV

    open (unit=Unitout, file='Mob.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'Mob.dat'
    call Errors_Fileopen(err,var)

    !Mobility Plot - Electric field, Mobility, Drift Velocity
    do n = 1, NV
      write(Unitout,*) MV(n,1), MV(n,2), MV(n,2)*MV(n,1)
    end do 

    close(Unitout)

end subroutine Output_mob

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Output_Hamil - Plot Hamiltonian, full(raw), main diag
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Output_Hamil(H)

    integer(kind=4) :: r, n, j, err
    real(kind=8) :: w, l
    character*120 :: var
    real(kind=8), dimension(:,:), intent(in) :: H

    open (unit=Unitout, file='Hraw.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'Hraw.dat'
    call Errors_Fileopen(err,var)
    open (unit=Unitout1, file='H.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'H.dat'
    call Errors_Fileopen(err,var)
 
    100 FORMAT (1000(F6.3,1X))

    do r=1, Np
       write(Unitout1,100) H(n,n) ! Write Diagonals of Matrix to File for plotting
    end do

    !H = H - transpose(H) !Check Hermiticity
    do r=1, Np
       write(Unitout,*) H(r,:) ! Write Raw Hamiltonian to File
    end do    

    close(unit=Unitout)
    close(unit=Unitout1)


end subroutine Output_Hamil

end module
