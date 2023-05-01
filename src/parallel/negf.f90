!  negf.f90 
!
!  Authors: G. Walker - Vanderbilt University
!           T. Musho - Vanderbilt University
!
!  Created: 12/09/08
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
!  PROGRAM: Negf(Module)
!
!  PURPOSE: Contains subroutine to solve NEGF
! 
!  SUBROUTINES: Negf_Main - Allocates lapack varibles and begins main loop
!               Negf_Ivchar - Calculate Seebeck, Electrical Conduct, PowerFactor
!               Negf_Voltsrch - Interation of voltage to find zero current
!               Negf_Spotent - Calculate self potential
!               Negf_Matmul - Matrix multiplication with matrix algebra simplification
!               Negf_Invers - Lapack Green's function inversion
!               Negf_Bound - Boundary conditions
!               
!******************************************************************************
!
module Negf

!Use global variables from module
use Startup
use Hamil
use Eigen
use Potent
use Scatter
use Device, ONLY : U, Uo
use Energy_amr
use Utility
use Interpolate
use Sparse
use Mpinterface
use Negf_Scatter
use Restart

implicit none

  !! Define Global Variables !!
  integer(kind=4) :: c0, c1, c2, c3, c4, c2t, lde, k, nsb, stp, cfr = 0, csi = 0
  real(kind=8) :: Its, Itr, Eold, Imv, dEo, Icb_o, chng_tprt = 1, tot_end=0.0
  real(kind=8), allocatable, dimension(:) :: IT, I1, I2, I3, Ud, NN, N1, Vr, Tint, Rh, Rh1, f_p, f_e, Ua
  real(kind=8), allocatable, dimension(:,:) :: A, A1, No, Tr, Tr1, Tr2, U1, Id1, Id2, Id3, IV, Gf, Gn1, MV
  complex(kind=8) :: ka1, ka2
  complex(kind=8), allocatable, dimension(:,:) :: Sig1, Sig2, Gam1, Gam2, SigIn1, SigIn2, SigOut1, SigOut2
  complex(kind=8), allocatable, dimension(:,:) :: G, Gtc, G1, G2, G3, G4, T
  real(kind=8), allocatable, dimension(:,:) :: G5
  complex(kind=8), allocatable, dimension(:) :: Gamp 
  real(kind=8), allocatable, dimension(:,:) :: Gn, Gp, SigInp, SigOutp, Ua1
  real(kind=8), allocatable, dimension(:,:,:) :: Aall, Gn_p, Gp_p, Gn_pi, Gp_pi, Si, Si1, Sit

  !! Define Private Variables
  real(kind=8), private :: f1, f2, f3, f4, f5, f5_o
  integer(kind=4), allocatable, dimension(:) :: rcount, disp, disp_mpi, rcountv, disp_mpiv, disp_np,rcount_np
  integer(kind=4), allocatable, dimension(:) :: rcount_scat, disp_scat
  integer(kind=8), private, allocatable, dimension(:) :: IPIV
  complex(kind=8), allocatable, dimension(:) :: WORK
  complex(kind=8), private, allocatable, dimension(:) :: DU,DL,D,DU2
  real(kind=8), private, allocatable, dimension(:) :: Trt1, Trt2, I1t, I2t, I3t, f_et
  real(kind=8), allocatable, dimension(:,:) :: At, Gnt, Gpt, Gft, Sigt

 real(kind=8), private :: t_rend, t_rstart

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Main - Allocates lapack varibles and begins main loop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Main

  integer(kind=4) :: ierr, n, j
  real(kind=8) ::  rn, tn

  call Negf_Lapack_Alloc

  LDA = Np
  c1 = 0; c2 = 0; c3 = 0
  nsb = 0; N1(:) = 0.0; Isb = 1; Icb = 0

  !This is for splitting up the scattering and mixing routines
  !Profiling says the mixing consumes 64% of the comp time
  if(NPROC .lt. Nb_p) then
    tn = floor(Nb_p/DBLE(NPROC))
    disp_np(:) = 0; rcount_np(:) = 0
    !assign 1 mix to each proc
    do n = 1, NPROC
      disp_np(n) = Nint(tn)*(n-1)
      rcount_np(n) = Nint(tn)
    end do
    !remainder 
    rn = DBLE(Nb_p) - floor(tn*DBLE(NPROC))
    do n = 1, Nint(rn)
      disp_np(n) = disp_np(n) + 1
      rcount_np(n) = rcount_np(n) + 1
    end do
    !increase remaining displacements
    disp_np(Nint(rn)+1:NPROC) = disp_np(int(rn)+1:NPROC) + Nint(rn)
    disp_np(1) = 0 !always zero
  else
    !if more procs then length
    disp_np(:) = -1; rcount_np(:) = -1
    do n = 1, NPROC
      disp_np(n) = n
      rcount_np(n) = 1
    end do
    disp_np(1) = 0 !always zero
  end if

  !Determine Leading Dimension of Energy Array
  do n = 1, NPROC !Setup MPI send array
    !each client will do a subset of the energy array
    disp(n) = (En * (n-1)) + 1 !for energy array split - k variable
    disp_mpi(n) = Np*(En * (n-1)) !for mpiscatterv - 2d array
    rcount(n) = Np*En !for mpiscatterv - 2d array
    disp_mpiv(n) = (En * (n-1)) !for mpiscatterv - 1d array
    rcountv(n) = En !for mpiscatterv - 1d array
    if(disp_np(n) .ge. 0) then
      disp_scat(n) = disp_np(n)*Np*NE !for mpiscatterv - scatter array
      rcount_scat(n) = rcount_np(n)*Np*NE !for mpiscatterv - scatter array
    end if
  end do
  disp_scat(1) = 0 !always zero

  call Negf_Work_Opt !Optimize Workspace

  !Need to send Eigenvalues to all proc, were only calculated on root
  if(psc .eq. 1) call MPI_Bcast(DH,Np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  !electron

  if(psc .eq. 1 .and. tst .lt. 8) then
    if(rank .ne. 0) allocate(DH_p(Nb_p), stat=ierr) !was never allocated
    call MPI_Bcast(dE_p,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) !phonon
    call MPI_Bcast(DH_p,Nb_p,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) !phonon
  end if

  if(tst .lt. 8) then
    call MPI_Bcast(NE_p,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) !phonon
    call MPI_Bcast(Eo_p,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) !phonon
    call MPI_Bcast(Ef_p,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) !phonon
    call Startup_Energy_ph !Creates linespace of phonon, didn't bring over so have to construct
    if(rank .eq. 0) then
      call Negf_Inter_Trans
      call MPI_Barrier(MPI_COMM_WORLD,ierr)  !Sync
    else
      call MPI_Barrier(MPI_COMM_WORLD,ierr)  !Sync
    end if
    if(psc .eq. 1) call Negf_Scatp_Init
  end if


  !this is where the root and client split
  if(RANK .eq. 0) then
    call Negf_Root
  else
    call Negf_Client
  end if

  call Negf_Lapack_Dealloc

end subroutine Negf_Main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Scatp_Init - Calculate the inital phonon scattering matrices
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Scatp_Init

  integer(kind=4) :: ierr, n, j

  !calculate phonon scat matricies outside scattering routine
  !do this upfront because it is constant and wasteful to do it each iter
  if(tst .eq. 2 .or. tst .eq. 3) then
  if(rank .eq. 0) then !only root because I don't want to send phonon DOS to all
    phe(:,:) =0.0;pha(:,:)=0.0
    do j=1,Nb_p
    do n=1,Np
      !The phonon density of states has units of 1/omega^2 we need to convert to a number
      !without an energy association
      !Because the energy is interpolated use dE of the Energy spectrum
      phe(n,j) = sum(Gn_pi(n,:,j))*(dE_p/hbar_eV)**2/(2.0*pi) !Convert to a number FIXME: Noninteger values?
      pha(n,j) = sum(Gn_pi(n,:,j)-Gp_pi(n,:,j))*(dE_p/hbar_eV)**2/(2.0*pi) !Convert to a number N+1
    end do
    end do
    call MPI_Barrier(MPI_COMM_WORLD,ierr)  !Sync
    call MPI_Bcast(pse,Np*NE,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(pha,Np*NE,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  else
    phe(:,:) =0.0;pha(:,:)=0.0
    call MPI_Barrier(MPI_COMM_WORLD,ierr)  !Sync
    call MPI_Bcast(pse,Np*NE,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(pha,Np*NE,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  end if
  end if

end subroutine Negf_Scatp_Init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Inter_Trans - Need to do some shifting
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Inter_Trans

  integer(4) :: err, i, j

  real(kind=8), allocatable, dimension(:,:,:) :: Gn_pi_t, Gp_pi_t ! temp matrices

  !Interpolate phonon spectral data
  write(6,*)'** Interpolating Gn Phonon Spectra'
  ! Select type of interpolation in Interpolate Mod
  call Interpolate_ph2el(Gn_pi,Gn_p,NE,NE_p,Np,Np_p,Nb,Nb_p,Ef,Ef_p,Eo,Eo_p,dE,dE_p,Lt,debug_l) 
  !if(debug_l>5)call Interpolate_plot(Gn_pi,E,NE,Nb_p,Np,Gdx,'Gni_p.dat') ! plot interpolated phonon spectra to check against actual
  !if(debug_l>5)call Interpolate_plot(Gn_p,Ep,NE_p,Nb_p,Np_p,Gdx,'Gno_p.dat')! plot read-in phonon spectra
  write(6,*)'** Interpolating Gp Phonon Spectra'
  call Interpolate_ph2el(Gp_pi,Gp_p,NE,NE_p,Np,Np_p,Nb,Nb_p,Ef,Ef_p,Eo,Eo_p,dE,dE_p,Lt,debug_l) 
  !if(debug_l>5)call Interpolate_plot(Gp_pi,E,NE,Nb_p,Np,Gdx,'Gpi_p.dat')
  !if(debug_l>5)call Interpolate_plot(Gp_p,Ep,NE_p,Nb_p,Np_p,Gdx,'Gpo_p.dat')

  write(6,*)'** Transpose Gn and Gp Phonon Spectra'
  !FIXME:A little blunder, phonon spectra and electron spectra have rows and col interchanged.
  ! This reallocates the phonon spectra so we can simple multiply time electron spectra in code
  allocate(Gp_pi_t(NE,Np,Nb_p),Gn_pi_t(NE,Np,Nb_p), stat=err) !Lots of memory used
  Gp_pi_t = Gp_pi; Gn_pi_t = Gn_pi
  deallocate(Gp_pi,Gn_pi)
  allocate(Gp_pi(Np,NE,Nb_p),Gn_pi(Np,NE,Nb_p), stat=err)
  Gp_pi(:,:,:) = 0.0; Gn_pi(:,:,:) = 0.0
  do i = 1, Np
     do j = 1, NE
     Gp_pi(i,j,:)=Gp_pi_t(j,i,:); Gn_pi(i,j,:)=Gn_pi_t(j,i,:)
     end do
  end do
  deallocate(Gp_pi_t,Gn_pi_t) !Recover memory chunk

end subroutine Negf_Inter_Trans

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Lapack_Alloc - Allocate Lapack space
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Lapack_Alloc

  integer(kind=4) :: err
  character*120 :: var

  allocate(Trt1(NE), Trt2(NE), I1t(NE), I2t(NE), I3t(NE), Sigt(Np,En), stat=err)
  allocate(At(Np,NE),Gnt(Np,NE),Gpt(Np,NE), Gft(Np,NE), stat=err)
  allocate(Vr(5), N1(Np), f_et(NE), stat=err)
  Trt1(:)=0.0; Trt2(:)=0.0; I1t(:)=0.0; I2t(:)=0.0; I3t(:)=0.0; Sigt(:,:)=0.0;
  At(:,:)=0.0; Gnt(:,:)=0.0; Gpt(:,:)=0.0; Gft(:,:) = 0.0
  Vr(:)=0.0; N1(:)=0.0; f_et(:)=0.0

  allocate(WORK(LWORK), stat=err) !dlange needs work

  if(csp .eq. 1) then 
    allocate(IPIV(Np), stat=err) !Allocate LAPACK Variables
  end if

  if(csp .eq. 4) then 
    allocate(IPIV(Np), stat=err) !Allocate LAPACK Variables
    allocate(D(Np),DU(Np-1),DL(Np-1),DU2(Np-2), stat=err) !Allocate LAPACK Variables
    D(:)=0.0;DU(:)=0.0;DL(:)=0.0;DU2(:)=0.0
  end if

  allocate(rcount(NPROC), disp(NPROC), disp_mpi(NPROC), stat=err) 
  allocate(rcountv(NPROC), disp_mpiv(NPROC), stat=err)
  allocate(rcount_np(NPROC), disp_np(NPROC), stat=err)
  allocate(rcount_scat(NPROC), disp_scat(NPROC), stat=err)
  rcount(:)=0.0; disp(:)=0.0; disp_mpi(:)=0.0
  rcountv(:)=0.0; disp_mpiv(:)=0.0
  rcount_np(:)=0.0; disp_np(:)=0.0
  rcount_scat(:)=0.0; disp_scat(:)=0.0

  write(var,*) 'WORK(LWORK),IPIV(Np)'
  call Errors_Allocate(err,var) !Check Errors

  if(RANK .eq. 0) then
    call Potent_init ! Allocate potential temp vectors
    call Scatter_init ! Allocate scatter temp vectors
  end if

end subroutine Negf_Lapack_Alloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Lapack_Dealloc - Deallocate Lapack space
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Lapack_Dealloc

  integer(kind=4) :: err

  deallocate(Trt1, Trt2, I1t, I2t, I3t, stat=err)
  deallocate(At,Gnt,Gpt, Gft, stat=err)
  deallocate(Vr, N1, f_et, stat=err)
  deallocate(rcount, disp, disp_mpi, stat=err)
  deallocate(rcountv, disp_mpiv, stat=err)
  deallocate(rcount_np, disp_np, stat=err)
  deallocate(rcount_scat, disp_scat, stat=err)
  deallocate(WORK, stat=err)
  if(csp .eq. 1) deallocate(IPIV, stat=err)
  if(csp .eq. 4) deallocate(IPIV,D,DU,DL,DU2, stat=err)

  if(RANK .eq. 0) call Potent_deallocate ! Allocate potential temp vectors

end subroutine Negf_Lapack_Dealloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Client - Client Proc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Client

  integer(kind=4) :: i, ierr, n, c1_o = 1

  lde = (En * rank) + 1 !determine stop value
  i = 1 
  Its = 0; Itr = 0; f_e(:) = 0.0

  do while(i .ne. 0)

    call MPI_Bcast(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    if(i .eq. 4) then
      if(psc .eq. 1) call Negf_Scatp_Init
    end if
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    if(i .ne. 0) then !messy logic
      call MPI_Bcast(nsb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(Vr,5,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
      call MPI_Bcast(Ud,Np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(Ua,Np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      if(tst .lt. 8) then
      call MPI_Scatterv(SigInp,rcount,disp_mpi,MPI_DOUBLE_PRECISION,SigInp(:,lde:(lde + En - 1)),&
                        Np*En,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_Scatterv(SigOutp,rcount,disp_mpi,MPI_DOUBLE_PRECISION,SigOutp(:,lde:(lde + En - 1)),&
                        Np*En,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      end if
      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      Gn(:,:) = 0.0; Gp(:,:) = 0.0; 
      Rh(:) = 0.0; A(:,:) = 0.0; Gf(:,:) = 0.0  ! Initial electron density along length of channel
      IT(:) = 0.0; I1(:) = 0.0; I2(:) = 0.0; Tr(:,:) = 0.0 ! Zero for +=
      if(debug_l>0) Si(:,:,:) = 0.0; !monitor all scattering rates
      Sig1(:,:) = (0.0,0.0); Sig2(:,:) = (0.0,0.0); SigIn1(:,:) = (0.0,0.0); SigIn2(:,:) = (0.0,0.0)
      Ds = Vr(1);V = Vr(2); c1 = Nint(Vr(3)); c0 = Nint(Vr(4)); c4 = Nint(Vr(5))

      if(c4 .eq. 1) then
        f_e(lde:(lde + En - 1)) = 1e-12
        if(tst .ge. 8) f_e(lde:(lde + En - 1)) = 0.0
      end if

      do n = 1, Np
       U(n,n) = Ud(n) ! Place Array Along Diagonal
      end do

      k = lde; stp = 1

        do while(k .le. (lde + En - 1) .and. stp .ne. 0) ! Starting energy loop for incoming electron energy
             if(Er(k) .eq. 1 .or. amr_o .eq. 0) then !No refinement of energy
               amr_l = 0 ! Non-Amr Loop
               call Negf_Bound ! Boundary Conditions
               call Negf_Invers ! Matrix Inversion [1/eV]
               call Negf_Matmul ! Sparse Matrix Multiplication
             else !Refinement of energy
               amr_l = 1 ! Amr Loop
               Eold = E(k)
               dEo = dE
               dE = dE/Er(k)
               do n = 1, Er(k)
                 if(n .eq. 1) then
                   E(k) = Eold
                 else
                   E(k) = Eold*(1 + real((n-1)/Er(k)))
                 end if
                 call Negf_Bound ! Boundary Conditions
                 call Negf_Invers ! Matrix Inversion [1/eV]
                 call Negf_Matmul ! Sparse Matrix Multiplication
               end do
               E(k) = Eold
               dE = dEo
             end if

         !Check condition of fermi functions - No current contribition when equal
         !if(abs(f1-f2) .eq. 0 .and. c1 .ne. 1 .and. coe .eq. 1) stp = 0

         f_e(k) = f5
         k = k + 1 !increment loop
         cfr = 1 !first run

       end do !Energy loop

      ! Calculate energy refinement for AMR
      if(amr_o .ne. 0) then
        call energy_amr_refine_log(IT,Er)
      end if
   
      Its = sum(IT)
      call MPI_Reduce(Rh,NN,Np,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(Its,Itr,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if(tst .lt. 8) then
      call MPI_Gatherv(f_e(lde),En,MPI_DOUBLE_PRECISION,f_et,rcountv,disp_mpiv,MPI_DOUBLE_PRECISION, &
                       0,MPI_COMM_WORLD,ierr)
      call MPI_Gatherv(Gn(:,lde:(lde + En - 1)),Np*En,MPI_DOUBLE_PRECISION,Gnt,rcount,disp_mpi,MPI_DOUBLE_PRECISION, &
                       0,MPI_COMM_WORLD,ierr)
      call MPI_Gatherv(Gp(:,lde:(lde + En - 1)),Np*En,MPI_DOUBLE_PRECISION,Gpt,rcount,disp_mpi,MPI_DOUBLE_PRECISION, &
                       0,MPI_COMM_WORLD,ierr)
      end if
      if(debug_l>5) then
        call MPI_Gatherv(Tr(lde,1),En,MPI_DOUBLE_PRECISION,Trt1,rcountv,disp_mpiv,MPI_DOUBLE_PRECISION, &
                         0,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(Tr(lde,2),En,MPI_DOUBLE_PRECISION,Trt2,rcountv,disp_mpiv,MPI_DOUBLE_PRECISION, &
                         0,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(I1(lde),En,MPI_DOUBLE_PRECISION,I1t,rcountv,disp_mpiv,MPI_DOUBLE_PRECISION, &
                               0,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(I2(lde),En,MPI_DOUBLE_PRECISION,I2t,rcountv,disp_mpiv,MPI_DOUBLE_PRECISION, &
                         0,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(I3(lde),En,MPI_DOUBLE_PRECISION,I3t,rcountv,disp_mpiv,MPI_DOUBLE_PRECISION, &
                         0,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(A(:,lde:(lde + En - 1)),Np*En,MPI_DOUBLE_PRECISION,At,rcount,disp_mpi,MPI_DOUBLE_PRECISION, &
                       0,MPI_COMM_WORLD,ierr)
        call MPI_Gatherv(Gf(:,lde:(lde + En - 1)),Np*En,MPI_DOUBLE_PRECISION,Gft,rcount,disp_mpi,MPI_DOUBLE_PRECISION, &
                       0,MPI_COMM_WORLD,ierr)
      end if

      if(psc .eq. 1 .and. tst .lt. 8 .and. c3 .ne. 1) then
      select case (tst)
        case (1) ! Fully coupled one-way
          !call Negf_Scatter_Couple 
        case (2) ! Scalar DOS Inelastic Scattering
          if(psc .eq. 0) then
            call Negf_Scatter_Couple_Indv(DH,DH_p,Gnt,Gpt,Gn_pi,Gp_pi,Si,SigInp,SigOutp,chng,LDA,WORK)
          else
            !parallel mixing - faster
            call Negf_Scatter_Couple_Indv_parallel(DH,DH_p,Gnt,Gpt,Gn_pi,Gp_pi,Si,Sit,SigInp,SigOutp,chng,LDA,WORK, &
                                                   disp_scat,rcount_scat,rcount_np,disp_np)
          end if
        case default
       end select
      end if

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

    end if
  end do

end subroutine Negf_Client

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Root - Root Proc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Root

   integer(kind=4) :: n, ierr, i, j, sc_ee=1, c_scat=0, rstp
   character*120 :: op
   real(kind=4) :: t_rend, t_rendt=0.0
   real(kind=4), dimension(2) :: t_rstart, t_rstarto

   i=1 !MPI End-Client Flag
   t_rstart = MPI_Wtime() ! get time for restart dumps
   call dtime(t_rstart,t_rend)

   !! START OF MAIN CURRENT-VOLTAGE CALCULATIONS !! 
   write(6,*)
   if(debug_l > 0) write(6,*) 'START OF MAIN CURRENT-VOLTAGE CALCULATION' ! Write to screen
   write(op,'(A)') 'START OF MAIN CURRENT-VOLTAGE CALCULATION' !Write to file
   call debug_output(op)
   call flush(6)

   No(:,:) = 0.0; Itot_ke = 0.0; IV(:,:) = 0.0;
   c0 = 0; c1 = 0; c2 = 0; c3 = 0; c4 = 0; c2t=0
   if(vst .eq. 2) V = Vo
   f_e(:) = 0.0; 
   do while((dabs(dV/V) > sc_dv) .and. csi .ne. 1) ! Voltage Loop 

     c0 = c0 + 1
     c1 = 0; chng_tprt = 1; c2t = 0
     do while(chng_tprt .gt. sc_tp) !coupled scattered
     c1 = c1 + 1 ! Loop Counter
     c2 = 0 !Zero subband counter
     Itot = 0.0; Isb = 0.0; IT(:) = 0.0; Icb = 0.0; Rh(:) = 0.0; Rh1(:) = 0.0  ! Zero current total
     if(debug_l>5) then
       Tr1(:,:) = 0.0; Id1(:,:) = 0.0; Id2(:,:) = 0.0; A1(:,:) = 0.0; U1(:,:)= 0.0; Gn1(:,:) = 0.0
       Si1(:,:,:) = 0.0
     end if
     nsb = 0 ! Zero subband count

     do while(nsb < Nb .and. dabs(Isb/Icb) > sc_sb .or. c2 .eq. 0) ! Subband Loop 
      
       c2 = c2 + 1 ! Subband loop couter
       c2t = c2t + 1 !Total Subband
       c3 = 0 ! Zero potential counter
       nsb = nsb + 1  ! Next Subband
       if(debug_l > 0) write(6,'(A,I4,1X,ES12.5)') 'Subband #, E =',nsb,DH(nsb) ! Write to screen
       write(op,'(A,I4,1X,ES12.5)') 'Subband #, E =',nsb,DH(nsb)
       call debug_output(op)
       call flush(6)
       Ds = DH(nsb) ! Additional iterations over subbands
       chg = 1.0 ! Initial guess for potential change
       if(debug_l>5) then
         Tr1(:,nsb) = 0.0; Id1(:,nsb) = 0.0; Id2(:,nsb) = 0.0; Id3(:,nsb) = 0.0; U1(:,nsb) = 0.0
       end if
       U(:,:) = 0e-8; Ud(:) = 0 !initial guess at potential
       call Potent_reset !reset potential variables
       do n=1, Np
         Ua(n) = -(-V/(Np-1)*n + (V+V/(Np-1)))
       end do

       do while(chg > sc_po)! Checks potential change 

         c3 = c3 + 1 ! Potential loop counter
         c4 = 0
         chng = 1.0
         f_e(1:En) = 1e-12 !slight non-zero to allow pertubation into scatter
         if(tst .ge. 8) f_e(1:En) = 0.0
         if(tst .lt. 8) call Scatter_reset !reset potential variables
         !SigInp(:,:) = 1e-30; SigOutp(:,:) = 1e-30 !initial guess at scattering rate  
         SigInp(:,:) = 0e-30; SigOutp(:,:) = 0e-30 !initial guess at scattering rate  
         do while(chng > sc_sc) ! Scattering convergence loop

           Gn(:,:) = 0.0; Gp(:,:) = 0.0; 
           Rh(:) = 0.0; A(:,:) = 0.0; Gf(:,:) = 0.0  ! Initial electron density along length of channel
           IT(:) = 0.0; I1(:) = 0.0; I2(:) = 0.0; Tr(:,:) = 0.0 ! Zero for +=
           Sig1(:,:) = (0.0,0.0); Sig2(:,:) = (0.0,0.0); SigIn1(:,:) = (0.0,0.0); SigIn2(:,:) = (0.0,0.0)
           if(debug_l>0) Si(:,:,:) = 0.0; !monitor all scattering rates
           c4 = c4 + 1
           stp = 1; k = 1
           Vr(1) = Ds; Vr(2) = V; Vr(3) = DBLE(c1); Vr(4) = DBLE(c0); Vr(5) = DBLE(c4) !Setup variables to MPIsend

           !read restart state
           if(rst .eq. 1) then
             call Restart_read(Its,Itot,V,Imax,Vmax,Vneg,Vpos,Ineg,Ipos,Imin,c1,c2,c3,c4,A,Gf,Rh)
             rst = 0
           end if

           !t_rend = MPI_Wtime()
           t_rstarto = t_rstart
           call dtime(t_rstart,t_rend)
           tot_end = tot_end + t_rend
           !Check to see if restart should be written
           if(tor .eq. 0) then
             if(t_rendt .ge. rsr) then
               rstp = 1
             else
               t_rendt = t_rendt + t_rend
             end if
           else if(tor .eq. 1) then
             if(mod(c1,rss) .eq. 0 .and. c2 .eq. 1 .and. c1 .ne. 1) rstp = 1
           else
             rstp = 0 !don't stop for restart this time
           end if
           if(rstp .eq. 1) then
             c4 = c4 + 1 !increment restart counter
             call Restart_write(Its,Itot,V,Imax,Vmax,Vneg,Vpos,Ineg,Ipos,Imin,c1,c2,c3,c4,A,Gf,Rh)
             call dtime(t_rstart,t_rend) !restart time reset
             t_rendt = 0.0 !restart
             rstp = 0
           end if

           ! Mpi Send
           call MPI_Bcast(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
           if(i .eq. 4) then
             if(psc .eq. 1) call Negf_Scatp_Init
             i = 1 !reset
           end if
           call MPI_Barrier(MPI_COMM_WORLD,ierr)
           call MPI_Bcast(nsb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
           call MPI_Bcast(Vr,5,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
           call MPI_Bcast(Ud,Np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
           call MPI_Bcast(Ua,Np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
           if(tst .lt. 8) then
           call MPI_Scatterv(SigInp,rcount,disp_mpi,MPI_DOUBLE_PRECISION,Sigt, &
                             Np*En,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
           call MPI_Scatterv(SigOutp,rcount,disp_mpi,MPI_DOUBLE_PRECISION,Sigt, &
                             Np*En,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
           end if
           call MPI_Barrier(MPI_COMM_WORLD,ierr)
           do while(k .le. En .and. stp .ne. 0) ! Starting energy loop for incoming electron energy
             if(Er(k) .eq. 1 .or. amr_o .eq. 0) then ! No refinement of energy
               amr_l = 0 ! Non-Amr Loop
               call Negf_Bound ! Boundary Conditions
               call Negf_Invers ! Matrix Inversion [1/eV]
               call Negf_Matmul ! Sparse Matrix Multiplication
             else !Refinement of energy
               amr_l = 1 ! Amr Loop
               Eold = E(k)
               dEo = dE
               dE = dE/Er(k)
               do n = 1, Er(k)
                 if(n .eq. 1) then
                   E(k) = Eold
                 else
                   E(k) = Eold*(1 + real((n-1)/Er(k)))
                 end if
                 call Negf_Bound ! Boundary Conditions
                 call Negf_Invers ! Matrix Inversion [1/eV]
                 call Negf_Matmul ! Sparse Matrix Multiplication
               end do
               E(k) = Eold
               dE = dEo
             end if

             !Check condition of fermi functions - No current contribition when equal
             !if(abs(f1-f2) .eq. 0 .and. c1 .ne. 1 .and. coe .eq. 1) stp = 0

             f_e(k) = f5
             k = k + 1 !increment loop
             cfr = 1 !first run
           end do !Energy loop

         ! Calculate energy refinement for AMR
         if(amr_o .ne. 0) then
           call energy_amr_refine_log(IT,Er)
         end if

         Itr = 0; NN(:) = 0; Trt1(:) = 0; Trt2(:) = 0; I1t(:) = 0; I2t(:) = 0; I3t(:) = 0; At(:,:) = 0
         Gnt(:,:) = 0; Gpt(:,:) = 0.0; f_et(:) = 0.0; Gft(:,:) = 0
         Its = sum(IT)
         !Reduce = summation from all processors. Each processor has a subset of the total energy range
         if(NPROC .ne. 1) then
            call MPI_Reduce(Rh,NN,Np,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_Reduce(Its,Itr,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            if(tst .lt. 8) then
            call MPI_Gatherv(f_e,En,MPI_DOUBLE_PRECISION,f_et,rcountv,disp_mpiv,MPI_DOUBLE_PRECISION, &
                             0,MPI_COMM_WORLD,ierr)
            call MPI_Gatherv(Gn,Np*En,MPI_DOUBLE_PRECISION,Gnt,rcount,disp_mpi,MPI_DOUBLE_PRECISION, &
                             0,MPI_COMM_WORLD,ierr)
            call MPI_Gatherv(Gp,Np*En,MPI_DOUBLE_PRECISION,Gpt,rcount,disp_mpi,MPI_DOUBLE_PRECISION, &
                             0,MPI_COMM_WORLD,ierr)
            end if
            if(debug_l>5) then
              call MPI_Gatherv(Tr(:,1),En,MPI_DOUBLE_PRECISION,Trt1,rcountv,disp_mpiv,MPI_DOUBLE_PRECISION, &
                               0,MPI_COMM_WORLD,ierr)
              call MPI_Gatherv(Tr(:,2),En,MPI_DOUBLE_PRECISION,Trt2,rcountv,disp_mpiv,MPI_DOUBLE_PRECISION, &
                               0,MPI_COMM_WORLD,ierr)
              call MPI_Gatherv(I1,En,MPI_DOUBLE_PRECISION,I1t,rcountv,disp_mpiv,MPI_DOUBLE_PRECISION, &
                               0,MPI_COMM_WORLD,ierr)
              call MPI_Gatherv(I2,En,MPI_DOUBLE_PRECISION,I2t,rcountv,disp_mpiv,MPI_DOUBLE_PRECISION, &
                               0,MPI_COMM_WORLD,ierr)
              call MPI_Gatherv(I3,En,MPI_DOUBLE_PRECISION,I3t,rcountv,disp_mpiv,MPI_DOUBLE_PRECISION, &
                               0,MPI_COMM_WORLD,ierr)
              call MPI_Gatherv(A,Np*En,MPI_DOUBLE_PRECISION,At,rcount,disp_mpi,MPI_DOUBLE_PRECISION, &
                               0,MPI_COMM_WORLD,ierr)
              call MPI_Gatherv(Gf,Np*En,MPI_DOUBLE_PRECISION,Gft,rcount,disp_mpi,MPI_DOUBLE_PRECISION, &
                               0,MPI_COMM_WORLD,ierr)
            end if
            f_e = f_et !dump it back
         else
            NN = Rh; Itr = Its; Gnt = Gn; Gpt = Gp; f_et = f_e
            if(debug_l>5) then
              Trt1 = Tr(:,1); Trt2 = Tr(:,2); I1t = I1; I2t = I2; I3t = I3; At = A; Gft=Gf
            end if
         end if

           !Select type of scattering, Construct scattering matricies
           if(c3 .gt. 1) then !skip scattering first time to get potential solution
           select case (tst)
             case (1) ! Fully coupled one-way - Inelastic
               call Negf_Scatter_Couple(DH_p,Gnt,Gpt,Gn_pi,Gp_pi,Si,SigInp,SigOutp,chng,LDA,WORK)
             case (2) ! Scalar DOS Inelastic Scattering - Inelastic
               if(psc .eq. 0) then
                 call Negf_Scatter_Couple_Indv(DH,DH_p,Gnt,Gpt,Gn_pi,Gp_pi,Si,SigInp,SigOutp,chng,LDA,WORK)
               else
                 call Negf_Scatter_Couple_Indv_parallel(DH,DH_p,Gnt,Gpt,Gn_pi,Gp_pi,Si,Sit,SigInp,SigOutp,chng,LDA,WORK, &
                                                        disp_scat,rcount_scat,rcount_np,disp_np)
               end if
             case (3) ! Scalar DOS Elastic Scattering - Elastic
               call Negf_Scatter_Couple_Indv_Elast(DH_p,Gnt,Gpt,Gn_pi,Gp_pi,Si,SigInp,SigOutp,chng,LDA,WORK)
             case (4) ! Multiple phonon BE Dist Non-coupled - Inelastic
               call Negf_Scatter_Multiple(DH_p,Gnt,Gpt,Si,SigInp,SigOutp,chng,LDA,WORK)
             case (5) ! Multiple phonon BE Dist Non-coupled - Elastic
               call Negf_Scatter_Multiple_Elast(DH_p,Gnt,Gpt,Si,SigInp,SigOutp,chng,LDA,WORK)
             case (6) ! Single phonon BE Dist Non-coupled - Inelastic
               call Negf_Scatter_Single(Gnt,Gpt,Si,SigInp,SigOutp,chng,LDA,WORK)
             case (7) ! Single phonon BE Dist Non-coupled - Elastic
               call Negf_Scatter_Single_Elast(Gnt,Gpt,Si,SigInp,SigOutp,chng,LDA,WORK)
             case default
               !Do nothing, No Scattering, Ballistic
               chng = 0.0 !no scattering
           end select
           else
            chng = 0.0
           end if
         call MPI_Barrier(MPI_COMM_WORLD,ierr) !Align Processes
         end do ! End of scattering loop

         ! Calculate energy refinement for AMR
         if(amr_o .ne. 0) then
           call energy_amr_refine_log(IT,Er)
         end if

         ! Potential Calculation
         call Matxopt_trace(U,Np,Ud) ! Stores diagonal elements in Array (eV)

         ! Poisson convergence method
         select case (coc)
           case (1) ! Simple mixing
             call Potent_lin_smix(NN,chg,Ud,No(:,nsb))
           case (2) ! Anderson mixing
             call Potent_anderson(NN,chg,Ud,No(:,nsb))
         end select
         if(debug_l>2) write(6,'(/,A,ES11.4,1X,ES11.4)') 'chg, t/sb = ',chg,tot_end/dble(c2t) ! Write change of charge to screen
         call Matxopt_diag(Ud,Np,U) ! Place Array Along Diagonal
         if(c3 .lt. 2) chg = 1.0
       end do ! End of potential

       !Save data
       if(debug_l>5) then !don't need to save if not plotting
         do n = 1, Np ! Add ionized donors across subbands
           A1(n,:) = A1(n,:) + At(n,:) !Local density of states
           Gn1(n,:) = Gn1(n,:) + Gft(n,:) !Local density of states
           Aall(n,:,nsb) = At(n,:)
         end do
         do n = 1, NE !need to divide by fermi function to get scattering rate in 1/s
           do j = 1, Nb_p 
             Si1(:,n,j) = Si1(:,n,j) + Si(:,n,j)/(f_et(n))
           end do
         end do
         Tr1(:,nsb) = Trt1 ! Transmission for each subband
         Tr2(:,nsb) = Trt2 ! Transmission for each subband
         U1(:,nsb) = Ud !Electrostatic potential
         Id1(:,nsb) = dE*Io*I1t; Id2(:,nsb) = dE*Io*I2t !Total current/subband
         Id3(:,nsb) = dE*Io*I3t;
       end if

       Rh1 = Rh1 + NN !Electron density
       Isb = dE*Io*Itr !Subband current
       Itot = Itot + Isb !Total current
       if(c1 .eq. 1) Itot_ke = Itot_ke + dE*Io*Itr !Total energy
       if (c1 .eq. 1 .and. nsb .le. sc_num) then
         Icb = Icb + Isb
         Isb = Icb
         !write(*,*) 'in loop', Isb/Icb
       else
         Icb = Icb + Isb
       end if

     end do ! End Subband
  
     if(tst .lt. 8) then
       !prepare to send back to phonon
       Sit(:,1,:) = 0.0
       do n = 1, Np
         do j = 1, Nb_p
           do k = 1, NE
            Sit(n,1,j) = Sit(n,1,j) + Si1(n,k,j)/NE !1/s average scattered energy
           end do
         end do
       end do

   
     !two-way coupling.

       call Interpolate_ph1d(Sit,Gn_p,Np,Np_p,Nb_p,Lt) !interpolate spatial mesh
       call Mpinterface_Write_Scph(Gn_p,Np_p,Nb_p,sc_ee,debug_l) !pass scattered energy
       call Mpinterface_Read_Gn2(Gn_p,NE_p,Eo_p,Ef_p,Np_p,Nb_p,debug_l,chng_tprt)
       call Mpinterface_Read_Gp2(Gp_p,NE_p,Np_p,Nb_p,debug_l,chng_tprt)
       call Negf_Inter_Trans
       i = 4
     else
       chng_tprt = 0.0
     end if

     end do ! Coupled Scattered

     if(debug_l>5 .and. c0 .lt. NV+2) Ua1(:,c0) = Ua
     if(debug_l>2) write(StdOut,'(A,ES12.5,1X,ES12.5,1X,ES12.5,1X,I4,1X,ES12.5)') 'V, Itot, Isb/Icb, nsb, t/dV = ',V &
                   & ,Itot,DABS(Isb/Icb),nsb,tot_end/dble(c1)
     call Loop_Plot
     call Restart_write(Its,Itot,V,Imax,Vmax,Vneg,Vpos,Ineg,Ipos,Imin,c1,c2,c3,c4,A,Gf,Rh);c4 = c4 + 1 !increment restart counter
     if(debug_l>5) call State_Plot ! Zero current LDOS plot
     call flush(6)
     if(debug_l>5 .and. c1 .lt. NV+2) then
       IV(c0,1) = V; IV(c0,2) = Itot; MV(c0,1) = V/Lt;
       if(c1 .gt. 1) then 
         MV(c0,2) = (IV(c0,2)-IV(c0-1,2))/(IV(c0,1)-IV(c0-1,1))/(q*sum(Rh1))
       else
         MV(c0,2) = 0.0
       end if
     end if
     select case (vst)
       case (1) ! Voltage Search
         call Negf_Voltsrch
       case (2) ! IV Char
         V = V + dV
         !V = (1.0e-11)*10.0**dble(c1)
         !dV = V
         if(V .gt. Vf) dV = 0
         if(V .gt. Vf) csi = 1  

     end select
   end do ! End of volage loop

   i = 0 ! End client processes
   call MPI_Bcast(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   call Negf_Ivchar

   sc_ee = 0 !exit phonon code
   if(tst .lt. 8)call Mpinterface_Write_Scph(Gn_p,Np_p,Nb_p,sc_ee,debug_l) !send phonon exit flag

   if(debug_l>0) write(6,*) 'END OF MAIN CURRENT-VOLTAGE CALCULATION' !Write to screen
   write(op,'(A)') 'END OF MAIN CURRENT-VOLTAGE CALCULATION' !Write to file
   call debug_output(op)
   call flush(6)

end subroutine Negf_Root

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Ivchar - Calculate Seebeck, Electrical Conduct, PowerFactor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Ivchar

    use Output

    character*120 :: op

    select case (vst)
      case (1) ! Voltage Search
        S =  V/dT * 1E+6
        sigma = Imax/(Vmax-V) * Lt
        PF = S**2 * sigma * 1E-12
        ke = -Itot_ke*Lt/dT
        mob = sigma/(q*sum(Rh1))*Lt
      case (2) ! IV Char
        S =  0
        if(debug_l>0) write(6,*)  'Warning: Calculating Sigma based on last two values of IV curve'
        write(op,'(A)') 'Warning: Calculating Sigma based on last two values of IV curve' !Write to file
        call debug_output(op)
        sigma = (IV(c1,2)-IV(c1-1,2))/(IV(c1,1)-IV(c1-1,1))*Lt! Imax/(Vmax-V) * Lt
        PF = S**2 * sigma * 1E-12
        ke = -Itot_ke*Lt/dT
        mob = sigma/(q*sum(Rh1))
    end select

    if(debug_l>5) call Output_iv(IV,c0)
    if(debug_l>5) call Output_mob(MV,c0)

end subroutine Negf_Ivchar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Voltsrch - Interation of voltage
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Voltsrch

   real(kind=8) :: Vold

   !! Save current for 0 and dV
   if (V .eq. 0) then
     Imin = Itot
   else if (V > Vmax) then
     Vmax = V
     Imax = Itot
   end if
  
   !! Use a Newton method to bracket and pick a new voltage guess
   Vold = V
   if (Vpos > Vneg) then  ! solution already bracketed
     if (Itot > 0.0) then
       Vpos = V
       Ipos = Itot
     else
       Vneg = V
       Ineg = Itot
     end if

     V = Vneg - Ineg/(Ipos-Ineg) * (Vpos-Vneg)

   else  ! keep increasing until solution is bracketed
     if (Itot < 0.0) then
       Vneg = V
       Ineg = Itot
       V = V + dV
     else
       Vpos = V
       Ipos = Itot
       V = Vneg - Ineg/(Ipos-Ineg) * (Vpos-Vneg)
     end if
   end if

   dV = V - Vold

end subroutine Negf_Voltsrch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Matmul - Matrix multiplication with matrix algebra simplification
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Matmul
   
   integer(kind=4) :: n, r
   real(kind=8) :: G1_i, G2_i, G3_i, G4_i

   G1(:,:) = (0,0); G2(:,:) = (0,0); G3(:,:) = (0,0); G4(:,:) = (0,0); G5(:,:) = 0.0

   Gtc = transpose(conjg(G)) ! Take Transpose Conjg once and save for multiple uses [1/eV]

   ! Divergent AMR statement
   if(amr_l .eq. 0) then

     G1(1,1) =  Gam1(1,1)
     G2(Np,Np) = Gam2(Np,Np)
     !G3 = matmul(matmul(G1,G),matmul(G2,Gtc)) !FIXME: so faster multiply T12
     call Matxopt_zsff(G1,1,G,Np,G3)
     call Matxopt_zsff(G2,Np,Gtc,Np,G4)
     call Matxopt_zffd2(G3,G4,Np,G5)
     call Matxopt_tracesum(G5,Np,G3_i)

     !G4 = matmul(matmul(G2,G),matmul(G1,Gtc)) !FIXME: so faster multiply T21
     call Matxopt_zsff(G2,Np,G,Np,G3)
     call Matxopt_zsff(G1,1,Gtc,Np,G4)
     call Matxopt_zffd2(G3,G4,Np,G5)
     call Matxopt_tracesum(G5,Np,G4_i)

     I1(k) = I1(k) + G3_i*(f1-f2)
     I2(k) = I2(k) + G4_i*(f2-f1)
     Tr(k,1) = Tr(k,1) + G3_i
     Tr(k,2) = Tr(k,2) + G4_i
     G2(:,:) = 0.0
     do n = 1, Np
        !G2(n,n) = SigInp(n,k)/f1  !inital guess
        if(f5 .ne. 0.0) G2(n,n) = SigInp(n,k)/f5 !determine gamma_in_phonon
     end do

     !Butticker Probe 
     G1(:,:) = (0,0);
     G1(1,1) =  Gam1(1,1)
     !G3 = matmul(matmul(G2,G),matmul(G1,Gtc)) !FIXME: so faster multiply T31
     call Matxopt_zdff(G2,G,Np,G3)
     call Matxopt_zsff(G1,1,Gtc,Np,G4)
     call Matxopt_zffd2(G3,G4,Np,G5)
     call Matxopt_tracesum(G5,Np,G3_i)

     if(tst .lt. 8) then
     G1(1,1) = (0,0);
     G1(Np,Np) = Gam2(Np,Np)
     !G4 = matmul(matmul(G2,G),matmul(G1,Gtc)) !FIXME: so faster multiply T32
     call Matxopt_zdff(G2,G,Np,G3)
     call Matxopt_zsff(G1,Np,Gtc,Np,G4)
     call Matxopt_zffd2(G3,G4,Np,G5)
     call Matxopt_tracesum(G5,Np,G4_i)

     !Trying to determine fermi function of contact
     !Summation of current at all contacts equals zero
     f5_o = f5
     if(G3_i+G4_i .ne. 0.0) then
       f5 = (f1*G3_i+f2*G4_i-(I1(k)+I2(k)))/(G3_i+G4_i)
     end if

     I3(k) = 0.0
     !back calculate current at scattering contact
     !conservation of energy
     if(G3_i+G4_i .ne. 0.0) then
        I3(k) = -G3_i*f1 - G4_i*f2 + f5*(G3_i+G4_i)
     end if 

     !Gn
     G1(:,:) = (0,0);
     G1(1,1) =  DBLE(SigIn1(1,1)) !(eV)
     G1(Np,Np) = DBLE(SigIn2(Np,Np)) !(eV)

     !Sigoutp from phonon scattering already contains fermi function
     do n = 1, Np
       G1(n,n) = G1(n,n) + SigInp(n,k)*f5/f5_o !(1/m^2) simplify multiplication full matrix times diag 
     end do

     else
       G1(:,:) = (0,0);
       G1(1,1) =  SigIn1(1,1) !(eV/m2)
       G1(Np,Np) = SigIn2(Np,Np) !(eV/m2)

       f5 = 0; I3(:) = 0.0
     end if

     !Total current is from T12
     IT(k) = I1(k)

     !G2 = matmul(matmul(G,G1),Gtc) !FIXME: so faster multiply
     call Matxopt_zfdf(G,G1,Np,G3) ! (1/eV*eV/m2)
     call Matxopt_zffd(G3,Gtc,Np,G2) !(1/eV-m2)

     !Gn = correlation function
     do n = 1, Np
       Gn(n,k) = (G2(n,n)) !(1/eV-m2)
     end do

     do n = 1, Np
       A(n,k) = A(n,k) + DBLE(i*(G(n,n) - Gtc(n,n))) ! Channel avaliable local density of states matrix (1/eV)
       Gf(n,k) = Gf(n,k) + DBLE(Gn(n,k))/(f1+f2+f5_o) ! Channel density of states matrix (1/eV)
     end do

     do n = 1, Np ! Calculating electron density in the channel (1/m2)
       Gp(n,k) = Gp(n,k) + (A(n,k)*(f1+f2+f5_o) - Gn(n,k)) ! Unfilled states, take view point of contact 1 (1/eV-m^2)
       Rh(n) = Rh(n) + (dE/(2.0*pi))*Gn(n,k) ! Calculating electron density in the channel (1/m^2)
     end do

   else

     !FIXME copy above and divide by amr division

   end if

end subroutine Negf_Matmul

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Work_Opt - Optimize scaLAPACK Workspace
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Work_Opt

   integer(kind=4) :: err, INFO
   character*120 :: var

   if(RANK .eq. 0) write(*,*) '** Optimizing Lapack Workspace'

   !! Determine LAPACK Optimal WORK() size !!
   LWORK = -1 ! Determine Optimal LAPACK Workspace (Not sure this is necessary each iteration?)
   call ZGETRI(Np,G,LDA,IPIV,WORK,LWORK,INFO)
   LWORK=int(WORK(1)) ! In Work(1) the optimal LWORK is Returned
   deallocate(WORK) ! Deallocation
   allocate(WORK(LWORK), stat=err) ! Allocation Optimal WORK
   write(var,*) 'Optimizing Lapack Workspace'
   call Errors_Allocate(err,var)

end subroutine Negf_Work_Opt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Invers - Lapack Green's function inversion
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Invers

   integer :: i
   integer(kind=4) :: INFO, Np_t, ca = 0
   real(kind=4) :: t_end
   real(kind=4), dimension(2) :: t_start

   if(csp .eq. 1) then
     !lapack general solver
     !call dtime(t_start,t_end)
     !! This is the same LAPACK calls Matlab makes
     call ZGETRF(Np,Np,G,LDA,IPIV,INFO) !LAPACK ROUTINE - LU Factorization of Matrix
     if (INFO <  0) then   
       write (6,*) 'the i-th argument had an illegal value'
     stop 'ABORT'
     else if(INFO > 0) then
       write (6,*) 'U(i,i) is exactly zero.' ! The factorization has been completed, but the factor U is exactly singular, and division by zero will occur if it is used to solve a system of equations.'
       stop 'ABORT'
     end if
     call ZGETRI(Np,G,LDA,IPIV,WORK,LWORK,INFO) !LAPACK ROUTINE - Computes Inverse Using LU Results
     if (INFO <  0) then   
       write (6,*) 'the i-th argument had an illegal value'
       stop 'ABORT'
     else if(INFO > 0) then
       write (6,*) 'U(i,i) is exactly zero; the matrix is singular and its inverse could not be computed'
       stop 'ABORT'
     end if
     !call dtime(t_start,t_end)
     !write(*,*) 'Lapack Time (time) = ',t_end,' sec' !Write to file
   else if(csp .eq. 2) then
     Np_t = Np
     ! Not working, not open source. Get seg fault in closed library
     ! Suppose to be faster then superLU
     !call Sparse_zgsmp_init(Np_t,G)
     !call Sparse_zgsmp(Np_t,ca,G)
     stop
   else if(csp .eq. 3) then
     Np_t = Np
     !call dtime(t_start,t_end)
     if(cfr .eq. 0) call Sparse_superlu_init(Np_t,G) ! only initalize once
     call Sparse_superlu(Np_t,ca,G)
     !call dtime(t_start,t_end)
     !write(*,*) 'superLU Time (time) = ',t_end,' sec' !Write to file
   else if(csp .eq. 4) then
     !lapack tridiagonal solver
     call Matxopt_ztrace_sub(G,1,Np,DL)
     call Matxopt_ztrace_sup(G,1,Np,DU)
     call Matxopt_ztrace(G,Np,D)
     call ZGTTRF(Np,DL,D,DU,DU2,IPIV,INFO)
     if (INFO <  0) then   
       write (6,*) 'the i-th argument had an illegal value'
       stop 'ABORT'
     else if(INFO > 0) then
       write (6,*) 'U(i,i) is exactly zero; the matrix is singular'
       stop 'ABORT'
     end if
     !call Matxopt_zput_sub(DL,1,Np,G)
     !call Matxopt_zput_sup(DU,1,Np,G)
     !call Matxopt_zput_sup(DU2,2,Np,G)
     !call Matxopt_zput(D,Np,G)
     G1 = eye
     call ZGTTRS('N',Np,Np,DL,D,DU,DU2,IPIV,G1,Np,INFO)
     G = G1
     !call ZGETRI(Np,G,LDA,IPIV,WORK,LWORK,INFO) !LAPACK ROUTINE - Computes Inverse Using LU Results
     if (INFO <  0) then   
       write (6,*) 'the i-th argument had an illegal value'
       stop 'ABORT'
     else if(INFO > 0) then
       write (6,*) 'U(i,i) is exactly zero; the matrix is singular and its inverse could not be computed'
       stop 'ABORT'
     end if
   end if

end subroutine Negf_Invers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Bound - Boundary conditions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Bound

   integer(kind=4) :: n
   real(kind=8) :: tm_ss, tm_se

   tm_ss = tsc*(1/Gdx(1)**2)
   tm_se = tsc*(1/Gdx(Np)**2)

   !! Calculating Fermi functions of electrons entering through Silicon Source and Drain !! 
   f1 = No1*dlog(1+dexp((-(E(k)-mu + Ua(1)))/kT1)) ! 2D Source fermi function. Units of 1/m^2 
   f2 = No2*dlog(1+dexp((-(E(k)-mu + Ua(Np)))/kT2)) ! 2D drain fermi function. Units of 1/m^2
   f3 = No1*(1-dlog(1+dexp((-(E(k)-mu + Ua(1)))/kT1))) ! 2D Source fermi function. Units of 1/m^2
   f4 = No2*(1-dlog(1+dexp((-(E(k)-mu + Ua(Np)))/kT2))) ! 2D drain fermi function. Units of 1/m^2
   !f1 = 1/(1+dexp(((E(k)-mu - V/2))/kT1))
   !f2 = 1/(1+dexp(((E(k)-mu + V/2))/kT2))

   f5 = f_e(k)

   ka1 =  i*(1-(1 - ((dcmplx(E(k)) + zplus - dcmplx(Ds) - dcmplx(Ua(1)))/(2.0*dcmplx(tm_ss))))**2)**0.5
   ka1 = -i*cdlog((1 - ((dcmplx(E(k)) + zplus - dcmplx(Ds) - dcmplx(Ua(1)))/(2.0*dcmplx(tm_ss)))) + ka1) ! Real Wave vector of broadened waves at drain (ACOS can be written in LOG Format See Matlab Manual) 		
   Sig1(1,1) = dcmplx(-tm_ss)*zexp(i*ka1) ! Self energy of source (eV), CMPLX=Converts to Complex Data Type
   Gam1(1,1) = i*(Sig1(1,1) - conjg(Sig1(1,1))) ! Source broadening matrix (eV) 
   SigIn1(1,1) = Gam1(1,1)*dcmplx(f1) ! Inscattering term for source (eV/m^2)
   SigOut1(1,1) = Gam1(1,1)*dcmplx(f3)                                                    

   ka2 =  i*(1-(1 - ((dcmplx(E(k)) + zplus - dcmplx(Ds) - dcmplx(Ua(Np)))/(2.0*dcmplx(tm_se))))**2)**0.5 
   ka2 = -i*cdlog((1 - ((dcmplx(E(k)) + zplus - dcmplx(Ds) - dcmplx(Ua(Np)))/(2.0*dcmplx(tm_se)))) + ka2) ! Real Wave vector of broadened waves at drain (ACOS can be written in LOG Format See Matlab Manual)  
   Sig2(Np,Np) = dcmplx(-tm_se)*zexp(i*ka2) ! Self energy of drain (eV) dcmplx=Converts Real to Complex Data Type
   Gam2(Np,Np) = i*(Sig2(Np,Np) - conjg(Sig2(Np,Np))) ! Drain broadening matrix (eV) 
   SigIn2(Np,Np) = Gam2(Np,Np)*dcmplx(f2) ! Inscattering term for drain (eV/m^2)
   SigOut2(Np,Np) = Gam2(Np,Np)*dcmplx(f4)

   Gamp(:) = 0.0
   !Inscatter and Outscattering
   !To just include Inscattering only SigInp
   !write(*,*) 'SigInp(:,k) - SigOutp(:,k) = ',sum(SigInp(:,k))/f5,sum(SigOutp(:,k))/f5
   !Gamp = i*(SigInp(:,k) - SigOutp(:,k))/f1
   !Hilbert transformation to get real part
   !Gamp = 2*pi*-0.5*(SigInp(:,k)/f1 + SigOutp(:,k)/f1) + 0.5*i*(SigInp(:,k)/f1 + SigOutp(:,k)/f1)
   !if(f5 .ne. 0) Gamp = -SigInp(:,k)/f5 - 0.5*i*(SigOutp(:,k))/f5
   if(f5 .ne. 0.0 .and. tst .lt. 8) Gamp = -SigInp(:,k)/f5-0.5*i*(SigInp(:,k)+SigOutp(:,k))/f5 !eV
   !if(f5 .ne. 0.0 .and. tst .lt. 8) Gamp = -SigInp(:,k)-0.5*i*(SigInp(:,k)+SigOutp(:,k)) !eV
   !if(f5 .ne. 0.0 .and. tst .lt. 8) Gamp = 0.5*i*(SigOutp(:,k))/f5 !eV
   !if(f5 .ne. 0) Gamp = SigInp(:,k) - 0.5*i*(SigOutp(:,k))
   !if(f5 .ne. 0) Gamp = -0.5*i*(SigOutp(:,k))/f5
   !write(*,*) 'sum(Gamp) =',sum(Gamp)

   G(:,:) = (0,0); 
   do n = 1, Np
     G1(n,n) = (0,0); G2(n,n) = (0,0)
   end do
   G1(1,1) = Sig1(1,1); G1(Np,Np) = Sig2(Np,Np) !(eV)

   do n = 1, Np
     G2(n,n) = ((dcmplx(E(k)) + zplus)*eye(n,n) - dcmplx(Ua(n)) - dcmplx(U(n,n))) !(eV) Simplify Addition Because Only Diagonal Terms 
   end do

   call Matxopt_ztcpy(dcmplx(H),G,Np)
   call Matxopt_sub_zdtf(G2,G,Np) 
   call Matxopt_sub_zdtf2(G1,G,Np)
   call Matxopt_sub_zvtf2(Gamp,G,Np)
   !G = (G2 - dcmplx(H) - G1) ! Green's function (1/eV) 

end subroutine Negf_Bound


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! State_Plot - Plot Density of State
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine State_Plot

    use Output

    call Output_ldos(A1,Aall,Gn1, c0)
    if(tst .lt. 8) call Output_Scin(Si1, Nb_p, c0)
    if(tst .lt. 8) call Output_tau(Si1, Nb_p, c0)
    call Output_tr(Tr1, nsb, c0)
    call Output_i_scat(Id1, Id2, Id3, c0)
    call Output_u(Uo, U1, Ubdy, Ua1, No, c0)
    call Output_rho(Rh1,No,c0)

end subroutine State_Plot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Loop_Plot - Appends to file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Loop_Plot

    character*120 :: op

    write(op,'(A,ES12.5,1X,ES12.5,1X,ES12.5,1X,I4,1X,ES12.5)') 'V, Itot, Isb/Icb, nsb, t/dV = ',V &
            & ,Itot,DABS(Isb/Icb),nsb,tot_end/dble(c1)
    call debug_output(op)

end subroutine Loop_Plot

end module
