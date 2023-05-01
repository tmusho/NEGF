!  startup.f90 
!
!  Authors: G. Walker - Vanderbilt University
!           T. Musho - Vanderbilt University         
!
!  Created: 3/30/09
!  Last Modified: 2/1/10
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
use Utility

implicit none

include 'mpif.h'

!! Define Constants !!
real(kind=8), parameter :: pi = 3.14159265
complex(kind=8), parameter :: i = (0.0,1.0) ! Imaginary Number
real(kind=8), parameter :: hbar = 1.0545716e-34  ! Reduced Plancks constant, Units of Joule-sec 
real(kind=8), parameter :: hbar_eV = 6.58211899e-16 !Units eV-s
real(kind=8), parameter :: q = 1.6021765e-19  ! Charge of an electron in Coulombs 
real(kind=8), parameter :: Io = 2.0*(q**2.0)/(2.0*pi*hbar)  ! Quantized conductance term, units of Amp/eV
real(kind=8), parameter :: kB = 8.6173422e-05  ! Boltzmann constant in eV/K
real(kind=8), parameter :: kBm = 1.3806503e-23  ! Boltzmann constant in J/K 
real(kind=8), parameter :: me = 9.1093819e-31  ! Electron mass in Kg

!! Define Global Variables !!
integer(kind=4) :: RANK, NPROC, En !Mpi variables
integer(kind=4) :: itype, nstartOld = 1, cse
integer(kind=4) :: amr_o, amr_d, amr_l, Nank, Nank_s, coe, cop, coc, coc_s, sct, csp, vst, tst, rst=0, psc
integer(kind=4) :: Np, NE, NV, Nb, Nl, Neb, New, Nec, mat1, mat2, off, mwa, Ns, Lg_o, rss, tor, sc_num
integer(kind=4) :: Np_p, NE_p, Nb_p, Ep_a = 0 !phonon
integer(kind=4), allocatable, dimension(:) :: Er, mtyp
real(kind=8) :: Eo_p, Ef_p, dE_p !phonon
real(kind=8) :: op_cut_si, op_cut_ge
real(kind=8) :: epsil0, epsilrS, epsilrG, epsilrJ, ms, mg, mu, Ecs, Ecg, Ecs_e, Ecg_e, x, op_cut, att, pcs, rsr
real(kind=8) :: a0Si, a0Ge, D001Si, D001Ge, G001Si, G001Ge, EdUSi, EdUGe, a11, a_1Si, eps_parSi, eps_perSi, Tp, Lc, Lg
real(kind=8) :: dEcs, dEcs_s, a_1Ge, eps_parGe, eps_perGe, dEcg, dEcg_s, Ecb, kT, dE, S, sigma, Pf, ms_b, mg_b
real(kind=8) :: kT1, kT2, an, ts, tsc, tg, No1, No2, Ds, V, chg, Eo, Ef, Vo, Vf, dT, dV, Isb, ac, Icb, sc_tp
real(kind=8) :: Nd_Si, Nd_Ge, Nc_Si, Nc_Ge, ab, aw, Lb, Lw, Vneg, Vpos, Ineg, Ipos, Vmax, Vold, Imax, Imin, Itot, Lt, ke, Itot_ke
real(kind=8) :: t_potent, smix, smix_s, wmix, wmix_s, sc_ip, sc_dv, sc_sb, sc_po, sc_sc, Nw, chng, Eph_si, Eph_ge, So, mob
real(kind=8), allocatable, dimension(:) :: E, Ep, Gdx
complex(kind=8) :: zplus
complex(kind=8), allocatable, dimension(:,:) :: eye

character*80 :: Rname
character(53) :: release='+++++ slpf(parallel) ver.3.0.0 (8/25/2011) +++++'

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Main - Main Startup
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Main

   integer(kind=4) :: err

   !call Startup_Var
   if(itype .eq. 1) then
     call Startup_Default
   else
     call Startup_Default
     if(rst .eq. 1) call Utility_GetFname(Input_fname)
     call Startup_Input
   end if
   call MPI_BARRIER(MPI_COMM_WORLD,err) 
   call Startup_UnivC !universal variables not matter input type
   call Startup_Material ! Select Material
   call Startup_Discrz ! Discretize Finite Element Mesh
   call Startup_Couple(ms,ab,ts) ! Calculate Inter-unit cell coupling energy - Si @ boundary to use in NEGF BC
   if(Nec .eq. 0) then
     tsc = ts
   else
     call Startup_Couple(ms,ac,tsc)
   end if
   call Startup_Couple(mg,aw,tg) ! Calculate Inter-unit cell coupling energy - Ge @ boundary to use in NEGF BC
   call Startup_Fermi(ms,kT1,No1) ! Constant used in Fermi function (source) (1/m2)
   call Startup_Fermi(ms,kT2,No2) ! Constant used in Fermi function (drain) (1/m2)
   call Startup_Energy
   call Startup_Ident
   if(RANK .eq. 0) call Startup_Writeparm
   if(RANK .eq. 0 .and. debug_l > 5 .and. rst .ne. 1) call Startup_Nodeplot !plot nodes to file
end subroutine Startup_Main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Input
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Input

   use input
   use readinput
   use inputformat

   real(kind=8) :: ms_e, mg_e

   character(len=LINE_LENGTH) :: headt

   call Utility_Load(release, NPROC, RANK)

   ! open input file
   call Input_Load(Input_fname, 0, 1)
   call InputFormat_Check('slpf')

   if(RANK .eq. 0) then
     call Debug_Title(release)
     call Startup_Output
   end if

   call ReadInput_Case(headt,cse)
   vst = cse

   call ReadInput_solver(csp)

   call ReadInput_Integration(amr_o,amr_d,cop,coe,sc_dv,sc_sb,sc_po,sc_ip,NE,Eo,Ef,an,sc_sc,sc_tp,sc_num)
   !write(*,*) 'Integration =',amr_o,amr_d,cop,coe,sc_dv,sc_sb,sc_po,sc_ip,NE,Eo,Ef,an

   call ReadInput_potential(coc,smix,wmix,nank)
   !write(*,*) 'Potential =',coc,smix,wmix,nank
   
   call ReadInput_strain(sct)
   !write(*,*) 'Strain =',sct

   call ReadInput_temp(Tp,dT)
   !write(*,*) 'Temp = ',Tp,dT

   call ReadInput_restart(tor,rsr,rss)
   !write(*,*) 'Restart = ',tor,rsr,rss

   call ReadInput_geometry(Nl,Lg,Lg_o,Lc)
   !write(*,*) 'Geo =',Nl,Lg,Lg_o,Lc

   call ReadInput_scatter(tst,psc,coc_s,smix_s,wmix_s,Nank_s)
   !write(*,*) 'Scat =',tst,psc,coc_s,smix_s,wmix_s,Nank_s

   call ReadInput_material1(mat1,Lb,Nd_Si,Nc_Si,epsilrS,ms_e,a0Si,D001Si,G001Si,EdUSi,Eph_si,op_cut_si)
   !write(*,*) 'Mat1 = ',mat1,Lb,Nd_Si,Nc_Si,epsilrS,ms_e,a0Si,D001Si,G001Si
   ms = me*ms_e

   call ReadInput_material2(mat2,Lw,Nd_Ge,Nc_Ge,epsilrG,mg_e,a0Ge,D001Ge,G001Ge,EdUGe,Eph_ge,op_cut_ge)
   !write(*,*) 'Mat2 = ',mat2,Lb,Nd_Ge,Nc_Ge,epsilrG,ms_e,a0Ge,D001Ge,G001Ge
   mg = me*mg_e
 
   call ReadInput_voltage(NV,Vo,Vf,Vneg,Vpos,Ineg,Ipos,Vmax,Imax,Imin,Nb)
   !stop
   !call Startup_Writeparm
   
   if(FluidID==0)then
      write(StdOut,*)'+++ Done Reading Input +++'
      write(StdOut,*)
   endif 

   close(UnitIn); 
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Output

   use input
   use readinput

   character(FILENAME_LENGTH) :: Outfile

   !  open Output File
   call Debug_SetTimeFloat(nstartOld)
   Outfile=Utility_MakeFileName('input',NPROC,nstartOld,RANK,'out')   
   open(UnitOut,file=OutFile)
   write(UnitOut,'(30x,a)')release        
   call Input_Echo(Input_fname,RANK)

   Outfile=Utility_MakeFileName('log',NPROC,nstartOld,RANK,'out')  
   open(UnitLog,file=OutFile)
   write(UnitLog,'(30x,a)')release

   call ReadInput_Debug(NPROC,nstartOld,RANK)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Var - Default Variables if Si/Ge default input
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Default

    character*120 :: var

    !! AMR Energy Inputs !!
    amr_o = 0 !AMR On = 1 AMR Off = 0
    amr_d = 10 !AMR Divisions
    
    !! Voltage selections
    ! Determine IV curve or search for zero current
    vst = 2 ! 1=Voltage Search (Seebeck) 2= IV Curve

    !! Sparse Matrix Inversion !!
    ! Using a sparse matrix factorization to determine LU components
    ! over 3 times speed up.
    csp = 4 ! 1=Lapack (dense solver), 2=WASP(not working) 3=SuperLU, 
            ! 4=Lapack (trisolve gaussian)

    !! Cut off Energy !!
    ! Will integrate energy range until no current contribution. Transmission plot cut off.
    coe = 0 ! Off = 0 On = 1

    !! Potential Calc
    cop = 0 ! 0 = hartree-fock approx, 1=full potential,

    !! Maxwell Approximation for Distributions
    ! 0 = no, 1=yes use maxwell
    mwa = 0

    !! Strain Effects !!
    sct = 0 ! 0=off 1=on

    !! Restart Parameter !!
    tor = 0 !type of restart 0=time interval, 1=voltage steps
    rsr = 30 !seconds between restarts
    rss = 5 !voltage steps between restarts
    
    !! Convergence Criteria Inputs !!
    sc_tp = 1e-4 !Temperature change criteria norm(dT)/norm(T)
    sc_dv = 1E-4 !dV/V Voltage Change over Voltage
    sc_sb = 1E-4 !Isb/Icb Subband Current Ratio
    sc_po = 1E-3 !chg Change in potential
    sc_ip = 1E-3 !inital chg, force higher convergence of inital calc
    sc_sc = 6E-1 !scattering convergence, this is a norm criterion not factor
                 !because it is norm it doesn't have to be set as low
    sc_num = 20  !force number of subbands @ V=0 to collect init charge dist

    !! Potential Poisson Inputs !!
    ! Convergence method
    coc = 1 ! 1=simple mixing; 2=anderson mixing
    smix = 0.3 ! simple mixing variable
    wmix = 10 ! anderson mix variable
    Nank = 5 ! number of previous vectors to keep, anderson mix

    !! Type of Scattering
    !  1= Fully coupled one-way, inelatic, 
    !  2= Multiple phonon Scalar DOS Inelastic Scattering 
    !  3= Multiple phonon Scalar DOS Elastic Scattering 
    !  4= Multiple phonon BE Dist Non-coupled Inelastic
    !  5= Multiple phonon BE Dist Non-coupled Elastic
    !  6= Single phonon BE Dist Non-coupled Inelastic
    !  7= Single phonon BE Dist Non-coupled  Elastic
    !  8=No scattering
    tst = 8

    !! Parallel Scattering/Mixing - Calculate scattering matricies in parallel
    !  Only scales up to 10procs
    if(tst .gt. 5) Nb_p = 1
    psc = 0 ! 0 = Off, 1= On

    !! Specifying phonon energy
    Eph_si = 0.0612 ! Longitudinal Intervalley Optical phonon energy of Silicon
    Eph_ge = 0.03704 ! Longitudinal Intervalley Optical phonon energy of Germanium
    !So = 0.009475 ! Phonon correlation factor in silicon obtained from LO intervalley def potential
    op_cut = 0.03 !Acoustic cut off freq, anything > is considered optical
    ! Convergence method
    coc_s = 1 ! 1=simple mixing; 2=anderson mixing
    smix_s = 0.3 ! simple mixing variable
    wmix_s = 10 ! anderson mix variable
    Nank_s = 7 ! number of previous vectors to keep, anderson mix

    !! Energy Inputs !!
    NE = 5201 ! Number of energy steps 
    Eo = -1.0 ! Minimux Potential Range
    Ef = 1.0 ! Maximum Potential Range

    !! Specifying temperature of source and drain !!
    Tp = 300 ! Temperature
    dT = 10 ! Temperature difference across device 
    
    !! Specifying Material Inputs !! 
    Nd_Si = Nd_Ge ! Assume doping concentration is the same for silicon and germanium
    Nc_Si = 6.5e25 ! Effective density of states silicon
    Nc_Ge = 3.4e25 ! Effective density of states germanium
    epsil0 = 8.8541878e-12  ! Permittivity of free space, Units of Farad/m = Coulomb/Volt-m 
    epsilrS = 11.7  ! Relative permittivity for Silicon 
    epsilrG = 16.0  ! Relative permittivity for Germanium 
    epsilrJ = 0.5*(epsilrS + epsilrG)  ! Relative permittivity of the Si/Ge junction  
    ms = 0.91*me  ! Silicon effective mass along 001 direction 
    mg = 0.95*me  ! Germanium effective mass along 001 direction 
    mu = 0.1  ! Fermi energy 

    ! Specifying Dopant Concentration
    Nd_Si = Nd_Ge ! Command Line input of dopant conc.
    Nc_Si = 2.8e19  ! Effective 3d density of states in silicon
    Nc_Ge = 1.04e19  ! Effective 3d density of states in germanium
    !Calculate effective density of states - bulk values, assume bulk contacts
    ms_b = 1.08*me !bulk effective mass
    mg_b = 0.55*me !bulk effective mass
    !This is temperature dependent, will modify conduction band edge
    Nc_Si = 2*(2*pi*ms_b*kBm*Tp/(hbar*2*pi)**2.0)**(3.0/2.0)/100**3.0  ! Effective 3d density of states in emitter
    Nc_Ge = 2*(2*pi*mg_b*kBm*Tp/(hbar*2*pi)**2.0)**(3.0/2.0)/100**3.0  ! Effective 3d density of states in medium

    call Startup_Cband(Nd_Si,Nc_Si,Ecs) ! Calculate Conduction Band - Si
    call Startup_Cband(Nd_Ge,Nc_Ge,Ecg) ! Calculate Conduction Band - Ge

    !! Grid Inputs !!
    call Startup_Size(ms,Eo,att) ! Plate
    !if(rank .eq. 0) write(*,*) 'Recommended Grid Size =',att
    an = att !1e-10  ! Target cell width
    !mat1 = 1; mat2 = 2 ! 1 = Silicon 2 = Germanium

    !! Input Parameters to Calculate Strain !! 
    a0Si = 5.4307  ! Silicon Lattice constant parallel to substrate 
    a0Ge = 5.6579  ! Germanium Lattice constant parallel to substrate 
    D001Si = 0.7713  ! Silicon deformation potential 
    D001Ge = 0.7509  ! Germanium deformation potential
    G001Si = 3.641 ! Silicon shear modulus
    G001Ge = 2.876 ! Germanium shear modulus

    EdUSi = 9.16  !Deformation potential
    EdUGe = 9.42  !Deformation potential

    !! Specifying applied drain bias !! 
    NV = 20 ! Number of voltage steps 
    Vo = 0.0 ! Minimum voltage range
    Vf = 0.1 ! Maximum voltage range
    Vneg = 0.0
    Vpos = -1
    Ineg = 0.0
    Ipos = 0.0
    Vmax = 0.0
    Imax = 0.0
    Imin = 0.0
    Nb = 31 ! Number of Subbands
    dV = (Vf-Vo)/NV

end subroutine Startup_Default

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_UnivC - Default parameters that are universal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_UnivC

    character*120 :: var

    !universal constants
    epsil0 = 8.8541878e-12  ! Permittivity of free space, Units of Farad/m = Coulomb/Volt-m 
    mu = 0.1  ! Fermi energy
    zplus = i*1E-12 ! Incremental term for energy for broadening during coupling

    call Startup_Checklayer

    !Temperature variables
    kT = kB*Tp ! General room temperature in eV (Boltzmann constant k * Temperature T) 0.0259
    kT1 = kB*(Tp-dT/2) ! Source temperature at 300K 
    kT2 = kB*(Tp+dT/2) ! Drain temperature at 300K

    !Calculate effective density of states in contacts
    Nc_Si = 2*(2*pi*ms*kBm*Tp/(hbar*2*pi)**2.0)**(3.0/2.0)/100**3.0
    Nc_Ge = 2*(2*pi*mg*kBm*Tp/(hbar*2*pi)**2.0)**(3.0/2.0)/100**3.0

    !Calculate conduction band edge
    call Startup_Cband(Nd_Si,Nc_Si,Ecs) ! Calculate Conduction Band - Si
    call Startup_Cband(Nd_Ge,Nc_Ge,Ecg) ! Calculate Conduction Band - Ge
 
    !Calculate strain
    if(sct .ne. 0) call Startup_StrainParam

    !! New conduction band edges due to shift in band edges with strain !!
    if(sct .eq. 1) then
      Ecs_e = Ecs + dEcs ! New conduction band edge in silicon substrate 
      Ecg_e = Ecg + dEcg ! New conduction band edge in germanium substrate
      Ecs = Ecs + dEcs_s ! New conduction band edge in silicon 
      Ecg = Ecg + dEcg_s ! New conduction band edge in germanium 
    else !No straining effects
      Ecs_e = Ecs ! New conduction band edge in silicon substrate 
      Ecg_e = Ecg ! New conduction band edge in germanium substrate
      Ecs = Ecs ! New conduction band edge in silicon 
      Ecg = Ecg ! New conduction band edge in germanium 
    end if
    Ecb = 0.5*(Ecs + Ecg) ! Conduction band at interface of Si and Ge

    111 format (A,F9.5,F9.5,F9.5,F9.5,F9.5,F9.5)
    !write(var,111) 'Ecs, Ecg, dEcs, dEcg, dEcs_s, dEcg_s = ', Ecs, Ecg, dEcs, dEcg, dEcs_s, dEcg_s
    !call Debug_Stdout(var)

    !Voltage constants
    dV = (Vf - Vo)/NV
    V = Vneg

end subroutine Startup_UnivC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Checklayer - Check layer size for graded case
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Startup_Checklayer

    character*120 :: var

    if(Nl*Lg .ge. Lb .or. Nl*Lg .ge. Lw) then
       write(var,*) 'Increase Layers or Decrease Gradient - Zero thickness layer'
       call Errors_Allocate(0,var)
    end if

end subroutine Startup_Checklayer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_StrainParam - Calculate Strain Offset
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_StrainParam

    integer(4) :: x
    character*120 :: var

    !! Calculation of substrate lattice constant !! 
    ! Concentration of silicon in Si(x)Ge(1-x) substrate
    if (mat1 .eq. 1) then
      x = 1.0 ! Si 001 substrate
    else 
      x = 0.0
    end if
    ! x = 0.5 ! Si(0.5)Ge(0.5) substrate 
    a11 = x*(a0Si) + (1-x)*(a0Ge) ! Lattice constant parallel to the superlattice layers  
    call Startup_Strain(a0Si,D001Si,a11,EdUSi,dEcs) ! Calculate Strain - Si
    call Startup_Strain(a0Ge,D001Ge,a11,EdUGe,dEcg) ! Calculate Strain - Ge

    a11 = (a0Si*Lb*G001Si + a0Ge*Lw*G001Ge)/(G001Si*Lb + G001Ge*Lw) ! Lattice constant parallel to the superlattice layers
    call Startup_Strain(a0Si,D001Si,a11,EdUSi,dEcs_s) ! Calculate Strain - Si
    call Startup_Strain(a0Ge,D001Ge,a11,EdUGe,dEcg_s) ! Calculate Strain - Ge  

    111 format (A,F9.5,F9.5,F9.5,F9.5,F9.5,F9.5)
    if(RANK .eq. 0) then
       write(var,111) 'Ecs, Ecg, dEcs, dEcg, dEcs_s, dEcg_s = ', Ecs, Ecg, dEcs, dEcg, dEcs_s, dEcg_s
       call Debug_Stdout(var)
    end if

    !! New conduction band edges due to shift in band edges with strain !!
    if(sct .eq. 1) then
      Ecs_e = Ecs + dEcs*pcs ! New conduction band edge in silicon substrate 
      Ecg_e = Ecg + dEcg*pcs ! New conduction band edge in germanium substrate
      Ecs = Ecs + dEcs_s*pcs ! New conduction band edge in silicon 
      Ecg = Ecg + dEcg_s*pcs ! New conduction band edge in germanium 
    else !No straining effects
      Ecs_e = Ecs ! New conduction band edge in silicon substrate 
      Ecg_e = Ecg ! New conduction band edge in germanium substrate
      Ecs = Ecs ! New conduction band edge in silicon 
      Ecg = Ecg ! New conduction band edge in germanium 
    end if
    Ecb = 0.5*(Ecs + Ecg) ! Conduction band at interface of Si and Ge

    zplus = i*1E-12 ! Incremental term for energy for broadening during coupling

end subroutine Startup_StrainParam

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Material - Select Material\Swap
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Material

    real(kind=8) :: f1, f2, f3, f4

    if(mat1 .eq. 1 .and. mat2 .eq. 2) then !Silicon Material 1
      Ecs = Ecs
      ts = ts
      Ecg = Ecg
      tg = tg
      Ecb = Ecb
      epsilrS = epsilrS
      tsc = ts
    else if(mat1 .eq. 2 .and. mat2 .eq. 1) then
      f1 = Ecs
      f2 = ts
      f3 = ms
      f4 = epsilrS
      ms = mg
      Ecs = Ecg
      epsilrS = epsilrG
      ts = tg
      Ecg = f1
      tg = f2
      mg = f3
      epsilrG = f4
      Ecb = Ecb
      tsc = tg
    else if(mat1 .eq. 1 .and. mat2 .eq. 1) then
      Ecs = Ecs; Ecs_e = Ecs_e
      ts = ts
      Ecg = Ecs; Ecg_e = Ecs_e
      tg = ts
      mg = ms
      epsilrG = epsilrS
      Ecb = Ecs
      tsc = ts
    else if(mat1 .eq. 2 .and. mat2 .eq. 2) then
      Ecs = Ecg; Ecs_e = Ecg_e
      ts = tg
      Ecg = Ecg; Ecg_e = Ecg_e
      tg = tg
      ms = mg
      epsilrS = epsilrG
      Ecb = Ecs
      tsc = tg
    end if

end subroutine Startup_Material

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Couple - Calculate Coupling Energy
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Couple(m,a,t)

    real(kind=8), intent(in) :: m, a
    real(kind=8), intent(out) :: t

    t = (hbar**2.0)/(2.0*m*q) ! Inter-unit cell coupling energy for silicon 

end subroutine Startup_Couple

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Couple - Calculate Required Grid Spacing
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Size(m,Ec,a)

    real(kind=8), intent(in) :: m, Ec
    real(kind=8), intent(out) :: a

    if(Ef .gt. Ec) then
      a = sqrt((hbar**2.0)/(2.0*m*abs(Ef-Ec)*q)) ! Inter-unit cell coupling energy for silicon
    else
      a = sqrt((hbar**2.0)/(2.0*m*abs(Ef-Eo)*q)) ! Inter-unit cell coupling energy for silicon
    end if

end subroutine Startup_Size


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Strain - Calculate Strain
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Strain(a0,D001,a11,EdU,dEc)

    real(kind=8), intent(in) :: a0, D001, a11, EdU
    real(kind=8), intent(out) :: dEc
    
    real(kind=8) :: a_1, eps_par, eps_per ! Local Varibles

    !!!!!!!!! Strain Parameters !!!!!!!!!! 
    a_1 = a0*(1.0-D001*((a11/a0)-1.0)) ! Lattice constant perpendicular to the layers 
    eps_par = (a11 - a0)/a0 ! Strain parallel to the layers 
    eps_per = (a_1 - a0)/a0 ! Strain perpendicular to the layers 
    dEc = (2.0/3.0)*(EdU)*(eps_per - eps_par) ! Shift in the conduction band edge due to strain

end subroutine Startup_Strain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Cband - Calculate Conduction Band
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Cband(Nd_Sc,Nc_Sc,Ed)

   real(kind=8), intent(in) :: Nd_Sc, Nc_Sc
   real(kind=8), intent(out) :: Ed
   character*120 :: var

   if(mwa .eq. 1 .or. (Nc_Sc/Nd_Sc-1) .le. 1) then
     !maxwell distribution approximation in this
     if(RANK .eq. 0) write(var,*) '**Using FD not MB distribution to determine Ec -> Ec-mu<kT'
     if(RANK .eq. 0) call Debug_Stdout(var)
     Ed = mu-kT*dlog(Nd_Sc/Nc_Sc)  ! Conduction band edge in silicon for doping level of 1018/cm3
   else
     Ed = mu+kT*dlog(Nc_Sc/Nd_Sc-1)  ! Conduction band edge in silicon for doping level of 1018/cm3
   end if

end subroutine Startup_Cband

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
    character*120 :: var

    !Neb = ceiling(Lb / an) ! Number of cells barrier
    call Startup_Gsize(Lb,an,Neb)
    ab = Lb / Neb
    !New = ceiling(Lw / an) ! Number of cells well
    call Startup_Gsize(Lw,an,New)
    aw = Lw / New
    if(Lc .ne. 0) then
       Nec = ceiling(Lc / an) ! Contact Length
       ac = Lc / Nec
    else
       Nec = 0
       ac = ab
    end if
    
    Np = (Nl+1)*(Neb) + Nl*(New) + 2*Nec + 1
    Nb = Np
    !! Construct Identity Matrix (1's Along Main Diagonal)!!
    allocate(Gdx(Np), mtyp(Np), stat=err) ! allocate Identity Matrix
    write(var,*) 'Gdx(Nd)'
    call Errors_Allocate(err,var)

    !!! diagonal Hamiltonian terms
    k=0
    do r = 1, Nl
       if(r .eq. 1) then  
          Gdx(1) = ac
          mtyp(1) = mat1
          k = 0
          if(Nec .ne. 0) then           
            do n = 2, Nec
              Gdx(n) = ac
              k=n
            end do
            Gdx(k+1) = ac
            mtyp(k+1) = mat1
          end if          
       end if

       do n = k+1, k+Neb
          Gdx(n) = ab
          mtyp(n) = mat1
          k=n
       end do

       do n = k+1, k+New
          Gdx(n) = aw
          mtyp(n) = mat2
          k=n
       end do


       if(r .eq. Nl) then
          do n = k+1, k+Neb
             Gdx(n) = ab
             mtyp(n) = mat1
             k=n
          end do
          if(Nec .ne. 0) then
            Gdx(k+1) = ac
            do n = k+2, k+Nec
              Gdx(n) = ac
              mtyp(n) = mat1
              k=n
            end do
          end if
          Gdx(Np) = ac
          mtyp(Np) = mat1
       end if
    end do

    Nb = Np

    Lt = 0
    do n = 1, Np-1
      Lt = Lt + Gdx(n)
    end do
end subroutine Startup_Discrz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Nodeplot - Plot Node file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Nodeplot

    integer(kind=4) :: n, err
    character*100 :: var

    open(unit=14, file='Node.dat', status='new', action='write', iostat=err)
    write(var,*) 'Node.dat'
    call Errors_Fileopen(err,var)

    do n = 1, Np
       write(14,*) Gdx(n)
    end do

    close(unit=14)

end subroutine Startup_Nodeplot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Fermi - Calculate Fermi Function Source/Drain
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Fermi(m,kT,No)

    real(kind=8), intent(in) :: m, kT
    real(kind=8), intent(out) :: No

    !Actually the 2d effective density of states
    No = (m*kT*q)/(2.0*pi*(hbar**2.0)) ! Constant used in Fermi function (1/m2)

end subroutine Startup_Fermi

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
! Startup_Energy_ph - allocate and constructs phonon energy array with gradient
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Energy_ph

    integer(kind=4) :: err
    character*100 :: var

    if(Ep_a .eq. 0) then
      allocate(Ep(NE_p), stat=err)
      write(var,*) 'Ep(NE)'
      call Errors_Allocate(err,var)
      Ep_a = 1
    end if

    Ep(:) = 0.0 ! Zero
    !write(*,*) 'NE_p,Eo_p,Ef_p =',NE_p,Eo_p,Ef_p 
    call linspace(Ep,NE_p,Eo_p,Ef_p) ! Specifying incoming electron range

    dE_p = Ep(2)-Ep(1) ! Difference between each energy step 

end subroutine Startup_Energy_ph

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Startup_Energy - allocate and constructs energy array with gradient
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Startup_Energy

    integer(kind=4) :: err
    character*120 :: var

    En = ceiling(real(NE/NPROC)) ! Calculate energys per task, integer value
    NE = NPROC * En ! Recalc number of energy levels

    allocate(E(NE), stat=err)
    write(var,*) 'E(NE)'
    call Errors_Allocate(err,var)

    allocate(Er(NE), stat=err) ! Allocation AMR vector
    write(var,*) 'Er(NE)'
    call Errors_Allocate(err,var)

    E(:) = 0.0 ! Zero

    !! Redfine Start of range slightly below Ec!! 
    !if(Ecs<Ecg .or. Ecs_e<Ecg_e) then
    !   if(Ecs<Ecs_e) then
    !      Eo = Ecs - 0.1*Ecs
    !   else
    !      Eo = Ecs_e - 0.1*Ecs_e
    !   end if
    !else
    !   if(Ecg<Ecg_e) then
    !      Eo = Ecg - 0.1*Ecg
    !   else
    !      Eo = Ecg_e - 0.1*Ecg_e
    !   end if
    !end if

    !Eo = minval([Ecs,Ecg,Ecs_e,Ecg_e])*(1-2)
    !Eo = -.4
    call linspace(E,NE,Eo,Ef) ! Specifying incoming electron range
    Er(:) = 1 ! Refinement Vector 1=no division

    Ns = Nint(Eo/dE) !start bin
    dE = abs(E(3)-E(2)) ! Difference between each energy step
 
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
    write(var(n),113) 'Doping Concentration = ',Nd_Ge;n = n + 1
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
    write(var(n),113) 'Permittivity of Free Space, epsil0 = ', epsil0;n = n + 1
    write(var(n),113) 'Relative Permittivity for Silicon, epsilrS = ', epsilrS;n = n + 1
    write(var(n),113) 'Relative Permittivity for Germanium, epsilrG = ', epsilrG;n = n + 1
    write(var(n),113) 'Relative Permittivity of the Si/Ge junction, epsilrJ = ', epsilrJ;n = n + 1
    write(var(n),113) 'Silicon effective mass, ms = ', ms;n = n + 1
    write(var(n),113) 'Germanium effective mass, ms = ', mg;n = n + 1
    write(var(n),113) 'Fermi Energy, mu = ', mu;n = n + 1
    write(var(n),113) 'Doping Concentration in Silicon, Nd_Si = ', Nd_Si;n = n + 1
    write(var(n),113) 'Doping Concentration in Germanium, Nd_Ge = ', Nd_Ge;n = n + 1
    write(var(n),113) 'Effective 3d Density of States in Silicon, Nc_Si = ', Nc_Si;n = n + 1
    write(var(n),113) 'Effective 3d Density of States in Germanium, Nc_Si = ', Nc_Si;n = n + 1
    write(var(n),114);n = n + 1
    write(var(n),114) 'Strain Parameter Inputs';n = n + 1
    write(var(n),114) '---------------------';n = n + 1
    write(var(n),113) 'Material 1 Lattice constant parallel to substrate, a0Si = ', a0Si;n = n + 1
    write(var(n),113) 'Material 2 Lattice constant parallel to substrate, a0Ge = ', a0Ge;n = n + 1
    write(var(n),113) 'Material 1  deformation potential, D001Si = ', D001Si;n = n + 1
    write(var(n),113) 'Material 2 deformation potential, D001Ge = ', D001Ge;n = n + 1
    write(var(n),113) 'Lattice constant parallel to the superlattice, a11 = ', a11;n = n + 1
    write(var(n),114);n = n + 1
    write(var(n),114) 'Bias Parameter Inputs';n = n + 1
    write(var(n),114) '---------------------';n = n + 1
    write(var(n),113) 'Conduction Band Edge Mat 1, Ecs = ', Ecs;n = n + 1
    write(var(n),113) 'Conduction Band Edge Mat 2, Ecg = ', Ecg;n = n + 1
    write(var(n),112) 'Number of Voltage Steps, NV = ', NV;n = n + 1
    write(var(n),113) 'Minimum Voltage Range, Vo = ', Vo;n = n + 1
    write(var(n),113) 'Maximum Voltage Range, Vf = ', Vf;n = n + 1
    write(var(n),112) 'Maximum Number of Subbands, Nb = ', Nb;n = n + 1
    write(var(n),112) 'Number of Energy Steps, NE = ', NE;n = n + 1
    write(var(n),113) 'Minimum Energy Range, Eo = ', Eo;n = n + 1
    write(var(n),113) 'Maximum Energy Range, Ef = ', Ef;n = n + 1
    write(var(n),113) 'Silicon Coupling Energy, ts = ', ts;n = n + 1
    write(var(n),113) 'Germanium Coupling Energy, tg = ', tg;n = n + 1
    write(var(n),114)
    call Debug_ParmOut(var,n)

end subroutine Startup_Writeparm

end module Startup
