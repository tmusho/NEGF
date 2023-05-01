!  startup.f90 
!
!  Authors: G. Walker - Vanderbilt University
!           T. Musho - Vanderbilt University         
!
!  Created: 3/30/09
!  Last Modified: 5/9/09
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
!******************************************************************************
!
!  PROGRAM: Startup(Module)
!
!  PURPOSE: Contains subroutine to initiate independent variables
! 
!  SUBROUTINES: Startup_Main - Main Startup, calls other subroutines in order
!               Startup_Energy - allocate and constructs energy array with gradient
!               Startup_Ident - Construct Identity Matrix
!               Startup_Var - Initiates Variables, contains dependent variables
!               Startup_Couple - Calculate Coupling Energy
!               Startup_Strain - Calculate Strain
!               Startup_Cband - Calculate Conduction Band
!               Startup_Discrz - Sets up discrization of matrix
!               
!******************************************************************************
!
module Startup

use Debug
use Errors

implicit none

!! Define Constants !!
real(kind=8), parameter :: pi = 3.14159265
real(kind=8), parameter :: hbar = 1.0545716e-34  ! Reduced Plancks constant, Units of Joule-sec 
real(kind=8), parameter :: hbar_eV = 6.58211899e-16 !Units eV-s
real(kind=8), parameter :: q = 1.6021765e-19  ! Charge of an electron in Coulombs 
real(kind=8), parameter :: kB = 8.6173422e-05  ! Boltzmann constant in eV/K
real(kind=8), parameter :: kBm = 1.3806503e-23  ! Boltzmann constant in J/K
real(kind=8), parameter :: me = 9.1093819e-31  ! Electron mass in Kg
complex(kind=8), parameter :: i = (0.0,1.0) ! Imaginary Number

!! Define Global Variables !!
integer(kind=4) :: RANK, NPROC, En, lde

integer(kind=4) :: amr_o, amr_d, amr_l, Nank, Nank_s, coe, cop, coc, sct, ab_s, ab_g, coc_s, csp, tst, psc, mwa, usr
integer(kind=4) :: Np, NE, NV, Nb, Nl, Neb, New, Nec, mat1, mat2, Lg_o, Neb_n, New_n, Neb_o, New_o, Naa, off, Ns
integer(kind=4), allocatable, dimension(:) :: Er, Neb_l, New_l
real(kind=8) :: Io, epsil0, epsilrS, epsilrG, epsilrJ, ms, mg, mu, Ecs, Ecg, Ecs_e, Ecg_e, x, am, Lt, ms_o, mg_o
real(kind=8) :: a0Si, a0Ge, D001Si, D001Ge, G001Si, G001Ge, EdUSi, EdUGe, a11, a_1Si, eps_parSi, eps_perSi, Tp, Lc
real(kind=8) :: dEcs, dEcs_s, a_1Ge, eps_parGe, eps_perGe, dEcg, dEcg_s, Ecb, kT, dE, S, sigma, Pf, Cv, Df
real(kind=8) :: kT1, kT2, an, ts, tsc, tg, No1, No2, Ds, V, chg, Eo, Ef, Vo, Vf, dT, dV, Isb, ac, Icb, rf, T1, T2
real(kind=8) :: Nd_Si, Nd_Ge, Nc_Si, Nc_Ge, ab, aw, Lb, Lw, Vneg, Vpos, Ineg, Ipos, Vmax, Vold, Imax, Imin, Itot, Qtot
real(kind=8) :: t_potent, smix, smix_s, wmix, wmix_s, sc_ip, sc_dv, sc_sb, sc_po, sc_sc, vsi, vge, Lg, Eph, So
real(kind=8), allocatable, dimension(:) :: E, Gdx, w, Gdx_p
complex(kind=8) :: zplus
complex(kind=8), allocatable, dimension(:,:) :: eye

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Main - Main Startup
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Main

   call Startup_Var
   call Startup_Material ! Select Material
   call Startup_Discrz ! Discretize Finite Element Mesh
   call Startup_Energy
   call Startup_Ident
   if(RANK .eq. 0) call Startup_Writeparm
   if(RANK .eq. 0 .and. debug_l > 5) call Startup_Nodeplot !plot nodes to file

   !! Specifying applied drain bias !! 
   Nb = Np ! Number of Subbands
   NV = Nb

end subroutine Startup_Main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Var - Initiates Variables
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Var

    character*100 :: var

    !! AMR Energy Inputs !!
    amr_o = 0 !AMR On = 1 AMR Off = 0
    amr_d = 10 !AMR Divisions

    !! Sparse Matrix Inversion !!
    ! Using a sparse matrix factorization to determine LU components
    ! over 3 times speed up.
    csp = 4 ! 1=Lapack (general solver), 2=WASP(not working) 3=SuperLU, 
            ! 4=Lapack (gaussian trisolve) 5=Lapack (direct trisolve)  

    !! Cut off Energy !!
    ! Will integrate energy range until no current contribution. Transmission plot cut off.
    coe = 0 ! Off = 0 On = 1

    !! Integration Sensing !!
    ! Stop integrating once current difference is zero
    ! Note full transmission range will not be captured
    coe = 0 ! 0=off 1=on

    !! Strain Effects !!
    sct = 1 ! 0=off 1=on

    !! Maxwell Approximation for Distributions
    ! 0 = no, 1=yes use maxwell
    mwa = 1

    !! Convergence Criteria Inputs !!
    sc_dv = 1E-6 !dV/V Voltage Change over Voltage
    sc_sc = 1e-2 !6E-2 !scattering convergence

    !! Parallel Scattering/Mixing
    psc = 1 ! 0 = Off, 1= On

    !! Type of Scattering
    !  1= Inelastic, 2= Elastic Scattering 3=Inelastic single phonon 4=Elastic single phonon 5=ballistic 
    tst = 1
    !! Specifying phonon energy
    So = 0.00005 !Phonon correlation factor
    ! Convergence method
    coc_s = 1 ! 1=simple mixing; 2=anderson mixing
    smix_s = 0.6 ! simple mixing variable
    wmix_s = 10 ! anderson mix variable
    Nank_s = 5 ! number of previous vectors to keep, anderson mix

    !! Grid Inputs !!
    !Because we are using a harmonic potential set these at equil position
    !then adjust the harmic potential to get the correct K
    !an = 2.715e-10  ! Target barrier width
    an = 5.43e-10
    !an = 1.17e-10
    !an = 3.379e-10 ! 5.43nm*0.6223 Special Points
    ab_s = 2 ! Atoms per basis
    !am = 2.825e-10  ! Target well width
    am = 5.65e-10
    !am = 3.51e-10 ! 5.65nm*0.6223 Special Points
    ab_g = 2 ! Atoms per basis
    !an = 3.6e-10
    !am = 3.6e-10
    !mat1 = 1; mat2 = 2 ! 1 = Silicon 2 = Germanium

    !! Specifying Material Inputs !! 
    ms = 1.6913e-26  ! Silicon effective mass along 001 direction 
    mg = 4.3744e-26  ! Germanium effective mass along 001 direction
    ms_o = ms; mg_o = mg 
    vsi = 9130 ! Silicon sound speed [m/s]
    vge = 5400 ! Germanium sound speed [m/s]

    !! Energy Inputs !!
    ! phonon excitation from contacts
    NE = 8201 ! Number of energy steps
    !User-specifed Energy range or dynamic range
    usr = 0 !0 = no 1 =yes; if yes use values below to define range
    !These values will be redefined once eigenvalues are calculated if dynamic used
    Eo = 1e-4*q ! Minimux Potential Range, units in joules
    Ef = 45e-3*q ! Maximum Potential Range [J]
    !dynamic factor below and above lowest eigenvalue
    rf = 0.537 !Range Extent above and below min and max eigenvalue [0-1]

    !! Specify Quantum Inputs !!
    Io = 1/(2.0*pi)**3  ! Quantized conductance term, units of Amp/eV

    !! Specifying temperature of source and drain !!
    !Tp = 300 ! Temperature
    !dT = 10 ! Temperature difference across device 
    kT = kBm*Tp ! General room temperature in J (Boltzmann constant k * Temperature T) 0.0259
    kT1 = kBm*(Tp-dT/2) ! Source temperature at 300K 
    kT2 = kBm*(Tp+dT/2) ! Drain temperature at 300K
    T1 = Tp-dT/2
    T2 = Tp+dT/2

    zplus = i*1e-12! 0.0001*(1-Eo/Ef) ! Incremental term for energy for broadening during coupling

end subroutine Startup_Var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_DOS - Calculate Fermi Function Source/Drain
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_DOS(v,kbT,Eph,No)

    real(kind=8), intent(in) :: v, kbT, Eph
    real(kind=8), intent(out) :: No

    ! 1d Phonon density of states
    No = kbT*Eph/(2*pi*v**2*hbar**2) !(Ep/hbar)**2/(2*pi**2*vsi**3) ! Constant used in Fermi function (1/m2)

end subroutine Startup_DOS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Material - Select Material\Swap
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Material

    real(kind=8) :: f1, f2, f3, f4, f5

    if(mat1 .eq. 1 .and. mat2 .eq. 2) then !Silicon Material 1
      Ecs = Ecs; Ecg = Ecg; Ecb = Ecb; epsilrS = epsilrS
    else if(mat1 .eq. 2 .and. mat2 .eq. 1) then
      f1 = Ecs; f2 = an; f3 = ms; f4 = epsilrS; f5 = vsi;
      ms = mg; vsi = vge
      Ecs = Ecg; epsilrS = epsilrG; Ecg = f1; an = am; am = f2
      mg = f3; epsilrG = f4; Ecb = Ecb; vsi = f5
    else if(mat1 .eq. 1 .and. mat2 .eq. 1) then
      Ecs = Ecs; Ecs_e = Ecs_e; Ecg = Ecs; Ecg_e = Ecs_e
      mg = ms; epsilrG = epsilrS; Ecb = Ecs; am = an; vsi = vsi
    else if(mat1 .eq. 2 .and. mat2 .eq. 2) then
      Ecs = Ecg; Ecs_e = Ecg_e; Ecg = Ecg; Ecg_e = Ecg_e
      ms = mg; epsilrS = epsilrG; Ecb = Ecs; an = am; vsi = vge
    end if

end subroutine Startup_Material

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Gsize - Sets up grid size
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Gsize(Lm,am,Nm)

   real(kind=8), intent(in) :: Lm, am
   integer(kind=4), intent(out) :: Nm

   !Round up a grid or down based on remainder
   if((Lm/am-floor(Lm/am)) .ge. 0.4) then 
     Nm = ceiling(Lm/am)
   else
     Nm = floor(Lm/am)
   end if

   if(Nm .eq. 0 .and. Lm .ne. 0) Nm = 1 !atleast one cell

   return

end subroutine Startup_Gsize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Discrz - Sets up discrization of matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Discrz

    integer(kind=4) :: r, k, n, err
    character*100 :: var

    !Neb = ceiling(Lb / an) ! Number of cells barrier
    call Startup_Gsize(Lb,an,Neb)
    Neb_o = Neb
    ab = Lb / Neb
    !New = ceiling(Lw / an) ! Number of cells well
    call Startup_Gsize(Lw,an,New)
    New_o = New
    aw = Lw / New
    if(Lc .ne. 0) then
       Nec = ceiling(Lc / an) ! Contact Length
       ac = Lc / Nec
    else
       Nec = 0
       ac = ab
       Lc = 0
    end if
    
    Np = (Nl+1)*(Neb) + Nl*(New) + 2*Nec + 1

    !! Construct Identity Matrix (1's Along Main Diagonal)!!
    allocate(Gdx(Np), Neb_l(Nl+1), New_l(Nl+1), stat=err) ! allocate Identity Matrix
    write(var,*) 'Gdx(Np),Neb_l(Nl),New_l(Nl)'
    call Errors_Allocate(err,var)
    Gdx(:) = 0.0;

    !!! diagonal Hamiltonian terms
    k=0
    do r = 1, Nl
       if(r .eq. 1) then  
          Gdx(1) = ac
          k = 0
          if(Nec .ne. 0) then           
            do n = 2, Nec
              Gdx(n) = ac
              k=n
            end do
            Gdx(k+1) = ac
          end if          
       end if

       do n = k+1, k+Neb
          Gdx(n) = Lb/Neb
          k=n
       end do

       do n = k+1, k+New
          Gdx(n) = Lw/New
          k=n
       end do

       if(Lg_o .ne. 0) then
          if(Lg_o .eq. 1) then
            Lb = Lb + Lg
            call Startup_Gsize(Lb,an,Neb_n)
            New_n = New_o
          else if(Lg_o .eq. 2) then
            Lw = Lw + Lg
            call Startup_Gsize(Lw,an,New_n)
            Neb_n = Neb_o
          else
            Lb = Lb + Lg
            Lw = Lw + Lg
            call Startup_Gsize(Lb,an,Neb_n)
            call Startup_Gsize(Lw,an,New_n)
          end if

          if(r .eq. NL) then 
             Np = Np+(Neb_n-Neb_o) !increase Np
          else
             Np = Np+(Neb_n-Neb_o)+(New_n-New_o)
          end if

          allocate(Gdx_p(Np), stat=err) ! define temp array
          Gdx_p(:) = 0.0;
          do n = 1, k
             Gdx_p(n) = Gdx(n)
          end do
          deallocate(Gdx)
          allocate(Gdx(Np), stat=err) ! increase array size
          Gdx = Gdx_p;
          deallocate(Gdx_p)
       
          write(var,*) 'Reallocation: Gdx(Np)'
          call Errors_Allocate(err,var)

          Neb_l(r) = Neb; New_l(r) = New

          if(r .eq. NL) then
            Neb_l(r+1) = Neb_n; New_l(r+1) = New_n
          else
            Neb_l(r+1) = Neb_n; New_l(r+1) = New_n
          end if

          Neb = Neb_n; New = New_n
       end if

       if(r .eq. Nl) then
          do n = k+1, k+Neb+1
             Gdx(n) = Lb/Neb
             k=n
          end do
          if(Nec .ne. 0) then
            Gdx(k+1) = ac
            do n = k+2, k+Nec
              Gdx(n) = Lb/Neb
              k=n
            end do
          end if
          Gdx(Np) = Lb/Neb
       end if
    end do

    !Compute total length
    Lt = 0
    do n = 1, Np-1
      Lt = Lt + Gdx(n)
    end do

    if(rank .eq. 0) write(*,*) 'Total Length, Lt = ',Lt

end subroutine Startup_Discrz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Nodeplot - Plot Node file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Nodeplot

    integer(kind=4) :: n, err
    character*100 :: var

    open(unit=14, file='Node_p.dat', status='new', action='write', iostat=err)
    write(var,*) 'Node_p.dat'
    open(unit=15, file='Neb_l_p.dat', status='new', action='write', iostat=err)
    write(var,*) 'Node.dat'
    open(unit=16, file='New_l_p.dat', status='new', action='write', iostat=err)
    write(var,*) 'Node_p.dat'
    call Errors_Fileopen(err,var)

    do n = 1, Np
       write(14,*) Gdx(n)
    end do

    do n = 1, Nl
       write(15,*) Neb_l(n)
       write(16,*) New_l(n)
    end do

    close(unit=14);close(unit=15);close(unit=16)

end subroutine Startup_Nodeplot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Ident - Startup Identity Matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Ident

    integer(kind=4) :: r, err
    character*100 :: var

    !! Construct Identity Matrix (1's Along Main Diagonal)!!
    allocate(eye(Np,Np), stat=err) ! Allocate Identity Matrix
    write(var,*) 'eye(Np,Np)'
    call Errors_Allocate(err,var)

    eye(:,:) = (0.0,0.0) ! Set Matrix to Zero
    do r = 1, Np
       eye(r,r) = (1.0,0.0) ! 1's Along diagonal
    end do

end subroutine Startup_Ident

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Energy_Redefine - redefines energy array with gradient
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Energy_Redefine(DH)

    real(8), dimension(:), intent(in) :: DH
    integer(kind=4) :: err
    real(kind=8) :: rfac
    character*100 :: var

    Eo = sqrt(DH(1))*hbar*(1.0-rf) !Set start energy below lowest eigenvalue
    rfac = sqrt(DH(1))*hbar-Eo !Set difference above eigenvalues same as below eigenvalues
    Ef = sqrt(DH(Np))*hbar+rfac !Set end energy above highest eigenvalue

    E(:) = 0.0 ! Zero
    
    call linspace(E,NE,Eo,Ef) ! Specifying incoming electron range
    Er(:) = 1 ! Refinement Vector 1=no division

    dE = E(3)-E(2) ! Difference between each energy step
    Ns = Nint(Eo/dE)

    if(RANK .eq. 0) then
      111 format (A,I6,E15.6,E15.6,E15.6)
      write(var,111) 'Redefine [J]: NE, Eo, Ef, dE = ',NE,Eo,Ef,dE
      call Debug_Stdout(var)
      write(var,111) 'Redefine [eV]: NE, Eo, Ef, dE = ',NE,Eo/q,Ef/q,dE/q
      call Debug_Stdout(var)
    end if

end subroutine Startup_Energy_Redefine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Energy - allocate and constructs energy array with gradient
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Energy

    integer(kind=4) :: err
    character*100 :: var

    En = ceiling(real(NE/NPROC)) ! Calculate energys per task, integer value
    NE = NPROC * En ! Recalc number of energy levels
    allocate(E(NE), stat=err)
    write(var,*) 'E(NE)'
    call Errors_Allocate(err,var)

    allocate(Er(NE), stat=err) ! Allocation AMR vector
    write(var,*) 'Er(NE)'
    call Errors_Allocate(err,var)

    E(:) = 0.0 ! Zero

    call linspace(E,NE,Eo,Ef) ! Specifying incoming electron range
    Er(:) = 1 ! Refinement Vector 1=no division

    dE = E(2)-E(1) ! Difference between each energy step
 
    if(RANK .eq. 0) then
      111 format (A,I6,E15.6,E15.6,E15.6)
      write(var,111) 'NE, Eo, Ef, dE = ',NE,Eo,Ef,dE
      call Debug_Stdout(var)
    end if

end subroutine Startup_Energy 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Write - Write parameters to file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Writeparm

    integer(kind=4) :: n
    character*80, dimension(80) :: var


    112 format (A,I6)
    113 format (A,E15.6)
    114 format (A)

    n = 1
    write(var(n),114) 'Command Line Inputs';n = n + 1
    write(var(n),114) '----------------------';n = n + 1
    write(var(n),112) 'Number of Processors = ',NPROC;n = n + 1
    write(var(n),112) 'Debug Level = ',debug_l;n = n + 1
    write(var(n),113) 'Barrier Length = ',Lb;n = n + 1
    write(var(n),113) 'Well Length = ',Lw;n = n + 1
    write(var(n),113) 'Contact Length, = ',Lc;n = n + 1
    write(var(n),113) 'Total Length, = ',Lt;n = n + 1
    write(var(n),112) 'Number of Layers = ',Nl;n = n + 1
    write(var(n),113) 'Silicon Lattice Constant = ',an;n = n + 1
    write(var(n),113) 'Germanium Lattice Constant = ',an;n = n + 1
    write(var(n),114);n = n + 1
    write(var(n),114) 'Quantum Inputs';n = n + 1
    write(var(n),114) '---------------------';n = n + 1
    write(var(n),113) 'hbar = ', hbar;n = n + 1
    write(var(n),113) 'q = ', q;n = n + 1
    write(var(n),113) 'Io = ', Io;n = n + 1
    write(var(n),113) 'kB = ', kB;n = n + 1
    write(var(n),113) 'me = ', me;n = n + 1
    write(var(n),114);n = n + 1
    write(var(n),114) 'Temperature Inputs';n = n + 1
    write(var(n),114) '---------------------';n = n + 1
    write(var(n),113) 'Tp = ', Tp;n = n + 1
    write(var(n),113) 'dT = ', dT;n = n + 1
    write(var(n),113) 'kT = ', kT;n = n + 1
    write(var(n),113) 'Source Temp, kT1 = ', kT1;n = n + 1
    write(var(n),113) 'Drain Temp, kT2 = ', kT2;n = n + 1
    write(var(n),114);n = n + 1
    write(var(n),114) 'Material Inputs';n = n + 1
    write(var(n),114) '---------------------';n = n + 1
    write(var(n),113) 'Silicon effective mass, ms = ', ms;n = n + 1
    write(var(n),113) 'Germanium effective mass, ms = ', mg;n = n + 1
    write(var(n),114);n = n + 1
    write(var(n),114) 'Bias Parameter Inputs';n = n + 1
    write(var(n),114) '---------------------';n = n + 1
    write(var(n),112) 'Maximum Number of Subbands, Nb = ', Nb;n = n + 1
    write(var(n),112) 'Number of Energy Steps, NE = ', NE;n = n + 1
    write(var(n),113) 'Minimum Energy Range, Eo = ', Eo;n = n + 1
    write(var(n),113) 'Maximum Energy Range, Ef = ', Ef;n = n + 1
    write(var(n),114);n = n + 1
    write(var(n),114) 'Scattering Inputs';n = n + 1
    write(var(n),114) '---------------------';n = n + 1
    write(var(n),113) 'Scattering Parameter = ', So;n = n + 1
    write(var(n),114)
    call Debug_ParmOut(var,n)

end subroutine Startup_Writeparm

end module Startup
