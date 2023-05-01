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
use Errors

implicit none

!! Define Global Variables !!
real(kind=8), dimension(2) :: Le
real(kind=8), allocatable, dimension(:,:) :: H, M, Ks

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Hamil_Main - Main Startup, subroutine calls
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Hamil_Main
   
   call Hamil_Allocate ! Allocate device matrix
   call Hamil_Create ! Input values into device matrix
   call Hamil_Mass_Matrix
   call Hamil_Assemble
   if(debug_l>5 .and. rank .eq. 0) call Hamil_Plot ! Plot hamiltonian diag and raw
   if(rank .eq. 0)call Debug_Hermi(H,Np) ! Check Hermicity

end subroutine Hamil_Main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Hamil_Assemble - Assemble Hamiltonian Matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Hamil_Assemble

   Ks = H
   H = Ks/M !Hamiltonian equals omega^2
   Ef = 2*sqrt(Ks(1,1)/2/M(1,1))*hbar
   if(rank .eq. 0) then
     write(*,'(A,ES15.6,A,ES15.6,A)') '  ** Mat 1 Spring Constant, Mass = ',Ks(1,1)/2,' N/m',M(1,1),' kg'
     write(*,'(A,ES15.6,A,ES15.6,A)') '  ** Mat 2 Spring Constant, Mass = ',Ks(Neb_l(1)+2,Neb_l(1)+2)/2,' N/m', &
   & M(Neb_l(1)+2,Neb_l(1)+2),' kg'
     write(*,'(A,ES15.6,A,ES15.6,A)') '  ** Redefine Cut-off energy,freq = ',Ef/q,' eV',max(2*sqrt(Ks(1,1)/2/M(1,1)), &
   & 2*sqrt(Ks(Neb_l(1) + 2,Neb_l(1)+2)/2/M(Neb_l(1)+2,Neb_l(1)+2))),' Hz'
   end if

end subroutine Hamil_Assemble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Hamil_Create - Create Hamiltonian Matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Hamil_Create

    integer(kind=4) :: r, n, k
    real(kind=8) :: t1, t2
    !(Ec_con,Ec_bar,Ec_wel,t_c,t_b,t_w,Nd,Nb,Nw,Gd)
    !integer(kind=8), intent(in) :: Nd, Nb, Nw
    !real(kind=8), intent(in) :: Ec_con, Ec_bar, Ec_wel, t_c, t_b, t_w
    !real(kind=8), intent(in), dimension(Nd) :: Gd

    Le(:) = 0.0;
    !!! diagonal Hamiltonian terms
    k=0; n=0
    do r = 1, Nl
       if(r .eq. 1) then
         call Hamil_Couple(ms,Gdx(1),t1)  
         H(1,1) = 2*t1*(1/Gdx(1)**2)           
         k = 0
         if(Nec .ne. 0) then
           do n = 2, Nec
             call Hamil_Couple(ms,Gdx(n),t1)
             call Hamil_Couple(ms,Gdx(n-1),t2)
             H(n,n) = (t1/Gdx(n)**2+t2/Gdx(n-1)**2)
             k=n
           end do
           call Hamil_Couple(ms,Gdx(k+1),t1)
           call Hamil_Couple(ms,Gdx(k),t2)
           H(k+1,k+1) = (t1/Gdx(k+1)**2+t2/Gdx(k)**2)
         end if
       end if

       do n = k+2, k+Neb_l(r)
          call Hamil_Couple(ms,Gdx(n),t1)
          call Hamil_Couple(ms,Gdx(n-1),t2) 
          H(n,n) = (t1/Gdx(n)**2+t2/Gdx(n-1)**2)
          k=n
       end do   
       call Hamil_Couple(mg,Gdx(k+1),t1)
       call Hamil_Couple(ms,Gdx(k),t2)
       H(k+1,k+1) = (t1*ms*ab_s/Gdx(k+1)**2+t2*mg*ab_g/Gdx(k)**2)

       do n = k+2, k+New_l(r)
          call Hamil_Couple(mg,Gdx(n),t1)
          call Hamil_Couple(mg,Gdx(n-1),t2) 
          H(n,n) = (t1/Gdx(n)**2+t2/Gdx(n-1)**2)
          k=n
       end do
       call Hamil_Couple(ms,Gdx(k+1),t1)
       call Hamil_Couple(mg,Gdx(k),t2)
       H(k+1,k+1) = (t1*mg*ab_g/Gdx(k+1)**2+t2*ms*ab_s/Gdx(k)**2)

       if(r .eq. Nl) then
          do n = k+2, k+Neb_l(r+1)
             call Hamil_Couple(ms,Gdx(n),t1)
             call Hamil_Couple(ms,Gdx(n-1),t2)
             H(n,n) = (t1/Gdx(n)**2+t2/Gdx(n-1)**2)
             k=n
          end do
          if(Nec .ne. 0) then
            call Hamil_Couple(ms,Gdx(k+1),t1)
            call Hamil_Couple(ms,Gdx(k),t2)
            H(k+1,k+1) = (t1/Gdx(k+1)**2+t2/Gdx(k)**2)
            do n = k+2, k+Nec
              call Hamil_Couple(ms,Gdx(n),t1)
              call Hamil_Couple(ms,Gdx(n-1),t2)
              H(n,n) = (t1/Gdx(n)**2+t2/Gdx(n-1)**2)
              k=n
            end do
          end if
          call Hamil_Couple(ms,Gdx(Np),t1)  
          H(Np,Np) = 2*t1*(1/Gdx(Np)**2)
       end if
    end do
    if(rank .eq. 0) then
    write(*,'(A,ES15.6,A)') '  ** Total Lattice Energy = ',Le(1),' J'
    write(*,'(A,ES15.6,A)') '  ** Mean Lattice Energy Per Atom = ',Le(1)/Np,' J/Atom'
    write(*,'(A,ES15.6,A)') '  ** RMS Lattice Energy Per Atom = ',Le(2)*2,' J/Atom'
    end if

    !!! off-diagonal Hamiltonian terms
    k=0
    do r = 1, Nl
       if(r .eq. 1) then
         call Hamil_Couple(ms,Gdx(n),t1)  
         H(2,1) = -t1*(1/Gdx(n)**2)
         H(1,2) = -t1*(1/Gdx(n)**2)  
         if(Nec .ne. 0) then
           do n = 2, Nec-1
             call Hamil_Couple(ms,Gdx(n),t1)  
             H(n,n+1) = -t1*(1/Gdx(n)**2)
             H(n+1,n) = -t1*(1/Gdx(n)**2)
             k=n
           end do
             call Hamil_Couple(ms,Gdx(k+1),t1)  
             H(k+1,k+2) = -t1*(1/Gdx(k+1)**2)
             H(k+2,k+1) = -t1*(1/Gdx(k+1)**2)
         end if
        end if

       do n = k+1, k+Neb_l(r)
         call Hamil_Couple(ms,Gdx(n),t1)  
         H(n,n+1) = -t1*(1/Gdx(n)**2)
         H(n+1,n) = -t1*(1/Gdx(n)**2)
         k=n
       end do

       do n = k+1, k+New_l(r)
         call Hamil_Couple(mg,Gdx(n),t1)  
         H(n,n+1) = -t1*(1/Gdx(n)**2)
         H(n+1,n) = -t1*(1/Gdx(n)**2)
         k=n
       end do

       if(r .eq. Nl) then
          do n = k+1, k+Neb_l(r+1)
             call Hamil_Couple(ms,Gdx(n),t1)  
             H(n,n+1) = -t1*(1/Gdx(n)**2)
             H(n+1,n) = -t1*(1/Gdx(n)**2)
             k=n
          end do
          if(Nec .ne. 0) then
            call Hamil_Couple(ms,Gdx(k+1),t1)  
            H(k+1,k+2) = -t1*(1/Gdx(k+1)**2)
            H(k+2,k+1) = -t1*(1/Gdx(k+1)**2)
            do n = k+2, k+Nec+1
              call Hamil_Couple(ms,Gdx(n),t1)  
              H(n,n+1) = -t1*(1/Gdx(n)**2)
              H(n+1,n) = -t1*(1/Gdx(n)**2)
              k=n
            end do
          end if
       end if
    end do

end subroutine Hamil_Create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Hamil_Mass_Matrix - Create Hamiltonian Mass Matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Hamil_Mass_Matrix

    integer(kind=4) :: r, n, k
    real(kind=8) :: t1, t2

    !!! diagonal Hamiltonian terms
    !! ab_s number of atoms per basis
    k=0; n=0
    do r = 1, Nl
       if(r .eq. 1) then
         M(1,1) = ms*ab_s           
         k = 0
         if(Nec .ne. 0) then
           do n = 2, Nec
             M(n,n) = ms*ab_s
             k=n
           end do
           M(k+1,k+1) = ms*ab_s
         end if
       end if

       do n = k+2, k+Neb_l(r) 
          M(n,n) = ms*ab_s
          k=n
       end do
       M(k+1,k+1) = ms*ab_s*mg*ab_g
       !M(k+1,k+1) = sqrt(ms*ab_s*mg*ab_g)
       !M(k+1,k+1) = sqrt((ms*ab_s)**2 + (mg*ab_g)**2)
       !M(k+1,k+1) = (ms*ab_s+mg*ab_g)/2
       do n = k+2, k+New_l(r)
          M(n,n) = mg*ab_g
          k=n
       end do
       M(k+1,k+1) = ms*ab_s*mg*ab_g
       !M(k+1,k+1) = sqrt(ms*ab_s*mg*ab_g)
       !M(k+1,k+1) = sqrt((ms*ab_s)**2 + (mg*ab_g)**2)
       !M(k+1,k+1) = (ms*ab_s+mg*ab_g)/2
       if(r .eq. Nl) then
          do n = k+2, k+Neb_l(r+1)
             M(n,n) = ms*ab_s
             k=n
          end do
          if(Nec .ne. 0) then
            M(k+1,k+1) = ms
            do n = k+2, k+Nec
              M(n,n) = ms*ab_s
              k=n
            end do
          end if
          M(Np,Np) = ms*ab_s
       end if
    end do
    !!! off-diagonal Hamiltonian terms
    k=0
    do r = 1, Nl
       if(r .eq. 1) then
         M(2,1) = ms*ab_s
         M(1,2) = ms*ab_s  
         if(Nec .ne. 0) then
           do n = 2, Nec-1
             M(n,n+1) = ms*ab_s
             M(n+1,n) = ms*ab_s
             k=n
           end do
             M(k+1,k+2) = ms*ab_s
             M(k+2,k+1) = ms*ab_s
         end if
        end if

       do n = k+1, k+Neb_l(r)
         M(n,n+1) = ms*ab_s
         M(n+1,n) = ms*ab_s
         k=n
       end do

       do n = k+1, k+New_l(r)
         M(n,n+1) = mg*ab_g
         M(n+1,n) = mg*ab_g
         k=n
       end do

       if(r .eq. Nl) then
          do n = k+1, k+Neb_l(r+1)
             M(n,n+1) = ms*ab_s
             M(n+1,n) = ms*ab_s
             k=n
          end do
          if(Nec .ne. 0) then
            M(k+1,k+2) = ms*ab_s
            M(k+2,k+1) = ms*ab_s
            do n = k+2, k+Nec+1
              M(n,n+1) = ms*ab_s
              M(n+1,n) = ms*ab_s
              k=n
            end do
          end if
       end if
    end do
end subroutine Hamil_Mass_Matrix

subroutine Hamil_Couple(m,a,t)

    integer(kind=8) :: tt
    real(kind=8) :: U
    real(kind=8), intent(in) :: m, a
    real(kind=8), intent(out) :: t

    if(m .eq. ms_o) then
      tt = 1
       !call Atompot_LJ(a,U,tt) !Call interatomic potential
       call Atompot_Harrison(a,U,tt)
       !call Atompot_Harrison_Spring(a,U,tt)
       !call Atompot_Hooke(a,U,tt)
       !call Atompot_Stillinger_Weber(a,U,tt)
    else
      tt = 2
       !call Atompot_LJ(a,U,tt) !Call interatomic potential
       call Atompot_Harrison(a,U,tt)
       !call Atompot_Harrison_Spring(a,U,tt)
       !call Atompot_Hooke(a,U,tt)
       !call Atompot_Stillinger_Weber(a,U,tt)

    end if
    t = 2*U ! Inter-unit cell coupling energy for silicon 
    Le(1) = Le(1) + U
    Le(2) = sqrt(Le(2)**2 + U**2/Np)


end subroutine Hamil_Couple

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Atompot_LJ - Lennard Jones potential. Values from Kittel
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Atompot_LJ(r,U,t)

    real(kind=8), dimension(2) :: e,sig
    integer(kind=8), intent(in) :: t
    real(kind=8), intent(in) :: r
    real(kind=8), intent(out) :: U

    !e = 14e-16 * 6.24150974e11 ! [eV] LJ Parameters
    e(1) = 49.1e-16 * 1e+7 ! [J] LJ Parameters Silicon
    e(2) = 47.2e-16 * 1e+7 ! [J] LJ Parameters Germanium
    sig(1) = 4.753e-10 ! [m] Silicon
    sig(2) = 5.109e-10 ! [m] Germanium

    U = (4*e(t)*((sig(t)/r)**12-(sig(t)/r)**6)) ! [J] Lennar-Jones potential 

end subroutine Atompot_LJ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Atompot_Hooke - Hooke's law potential. Values from Kittel
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Atompot_Hooke(r,U,t)

    real(kind=8), dimension(2) :: e,sig,c
    integer(kind=8), intent(in) :: t
    real(kind=8), intent(in) :: r
    real(kind=8), intent(out) :: U

    !e = 14e-16 * 6.24150974e11 ! [eV] LJ Parameters
    e(1) = 1.66e11 ! [N/m^2] Parameters Silicon
    e(2) = 1.29e11 ! [N/m^2] Parameters Germanium
    sig(1) = 3.835e-10 ! [m] Silicon
    sig(2) = 3.995e-10 ! [m] Germanium
    !rho(1) = ms/e(1)/vsi**2
    !rho(2) = mg/e(2)/vge**2
    c = 3*sig**3/16*e

    U = 0.5*c(t)*(r-sig(t))**2/sig(t)**2 ! [J] Lennar-Jones potential 

end subroutine Atompot_Hooke

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Atompot_Stillinger_Weber - Stillinger Weber potential with just pair part
!                            This potential is for the a/4(x + y + z) direction of
!                            the diamond structure. Not correct for 1d chain
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Atompot_Stillinger_Weber(r,U,t)

    real(kind=8) :: f,rr
    real(kind=8), dimension(2) :: e, A,B,p,q,aa,sig
    integer(kind=8), intent(in) :: t
    real(kind=8), intent(in) :: r
    real(kind=8), intent(out) :: U

    e(1) = 3.4723e-19 ! [J] Silicon
    e(2) = 3.4723e-19 ! [J] Germanium
    A(1) = 1.66e11 
    A(2) = 1.29e11 
    B(1) = 0.6022245584
    B(2) = 0.6022245584
    p(1) = 4
    p(2) = 4
    q(1) = 0
    q(2) = 0
    aa(1) = 1.8
    aa(2) = 1.8
    sig(1) = 2.0951e-10 ![m]
    sig(2) = 2.0951e-10 ![m]

    rr = r/sig(t)
    f = A(t)*(B(t)*rr**(-p(t))-rr**(-q(t)))*exp(1/(rr-aa(t)))
    write(*,*) 'rr, f =',rr,f

     U = -e(t)*f ! [J] Lennar-Jones potential

end subroutine Atompot_Stillinger_Weber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Atompot_Harrison - Harrison potential. Values from Mingo phonon papers
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Atompot_Harrison(r,U,t)

    real(kind=8), dimension(2) :: e,sig
    integer(kind=8), intent(in) :: t
    real(kind=8), intent(in) :: r
    real(kind=8), intent(out) :: U

    !e = 14e-16 * 6.24150974e11 ! [eV] LJ Parameters
    e(1) = 49.1*q ! [J] LJ Parameters Silicon
    e(2) = 47.2*q ! [J] LJ Parameters Germanium
    sig(1) = 2.35e-10 ! [m] Silicon
    sig(2) = 2.44e-10 ! [m] Germanium

    U = (0.5*e(t)*((r-sig(t))/sig(t))**2) ! [J] Lennar-Jones potential 

end subroutine Atompot_Harrison

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Atompot_Harrison_Spring - Harrison potential. Backout spring
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Atompot_Harrison_Spring(r,U,t)

    real(kind=8), dimension(2) :: e,e1,sig
    integer(kind=8), intent(in) :: t
    real(kind=8), intent(in) :: r
    real(kind=8), intent(out) :: U

    !e = 14e-16 * 6.24150974e11 ! [eV] LJ Parameters
    e(1) = 49.1*q ! [J] LJ Parameters Silicon
    e(2) = 47.2*q ! [J] LJ Parameters Germanium
    e1(1) = 1.07*q ! [J] LJ Parameters Silicon
    e1(2) = 0.845*q ! [J] LJ Parameters Germanium
    sig(1) = 2.35e-10 ! [m] Silicon
    sig(2) = 2.44e-10 ! [m] Germanium

    U = 0.5*4/3*(e(t))*r**2 ! [J] Lennar-Jones potential 
    
end subroutine Atompot_Harrison_Spring

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Atompot_Tersoff - Tersoff potential. Doesn't work
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Atompot_Tersoff(r,U,t)

    real(kind=8) :: fc
    real(kind=8), dimension(2) :: e,sig
    integer(kind=8), intent(in) :: t
    real(kind=8), intent(in) :: r
    real(kind=8), intent(out) :: U

    !e = 14e-16 * 6.24150974e11 ! [eV] LJ Parameters
    e(1) = 3.4*q !16.9e-16 * 1e+7 ! [J] LJ Parameters Silicon
    e(2) = 13.7e-16 * 1e+7 ! [J] LJ Parameters Germanium
    sig(1) = 2.74e-10 ! [m]
    sig(2) = 2.74e-10 ! [m]
   
    fc = 0.5-0.5*sin(pi/2*(r-sig(t)/r))
    U = 0.5*fc

end subroutine Atompot_Tersoff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Hamil_Allocate - Allocate Hamiltonian Matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Hamil_Allocate

    integer(kind=4) :: err
    character*100 :: var

    !! Hamiltonian Matrix, Units of eV !! 
    allocate(H(Np,Np), M(Np,Np), Ks(Np,Np), stat=err) ! Allocate Matrix
    write(var,*) 'H(Np,np)'
    call Errors_Allocate(err,var) !Check Errors

    H(:,:) = 0.0; Ks(:,:) = 0.0; M(:,:) = 1.0 !Define Hamiltonian Matrix to Zero

end subroutine Hamil_Allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Hamil_Allocate - Allocate Hamiltonian Matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Hamil_Plot

    integer(kind=4) :: err, r
    real(kind=8) :: l
    character*100 :: var

    open(unit=14, file='Hraw_p.dat', status='new', action='write', iostat=err)
    write(var,*) 'Hraw_p.dat'
    call Errors_Fileopen(err,var)
    open(unit=15, file='H_p.dat', status='new', action='write', iostat=err)
    write(var,*) 'H_p.dat'
    open(unit=16, file='Mraw_p.dat', status='new', action='write', iostat=err)
    write(var,*) 'Mraw_p.dat'
    open(unit=17, file='Kraw_p.dat', status='new', action='write', iostat=err)
    write(var,*) 'Kraw_p.dat'
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
      write(16,*) M(r,:) ! Write Raw Hamiltonian to File
      write(17,*) Ks(r,:) ! Write Raw Hamiltonian to File
    end do    

    close(unit=14); close(unit=15); close(unit=16); close(unit=17)

end subroutine Hamil_Plot

end module
