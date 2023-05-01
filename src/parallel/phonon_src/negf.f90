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
use Matxopt
use Scatter
use Energy_amr
use Sparse
use Grad
use Mpinterface
use Utility

implicit none
include 'mpif.h'

  !! Define Global Variables !!
  integer(kind=4) :: c0, c1, c2, c3, k, f, nsb, stp, cfr = 0, csi=0
  real(kind=8) :: Its, Itr, Eold, Imv, dEo, dw, chng, chng_tprt
  real(kind=8), allocatable, dimension(:) :: IT, I1, I2, I3, Rh, Tprt_old
  real(kind=8), allocatable, dimension(:,:) :: A, A1, Gf, Gf1, Gp1, Tr, Tr1, Id1, Id2, Id3, Rh1
  real(kind=8), allocatable, dimension(:,:) :: Gn, Gp, SigInp, SigOutp, f_e, Tprt, Tprtt, Tprt_NE
  real(kind=8), allocatable, dimension(:,:) :: nse, nsa, pse, psa, SigInpNew, SigOutpNew, SigInpNewt, SigOutpNewt
  real(kind=8), allocatable, dimension(:,:,:) :: Aall, Gnall
 
  complex(kind=8) :: ka1, ka2
  complex(kind=8), allocatable, dimension(:) :: Gamp 
  complex(kind=8), allocatable, dimension(:,:) :: Sig1, Sig2, Gam1, Gam2, SigIn1, SigIn2, SigOut1, SigOut2
  complex(kind=8), allocatable, dimension(:,:) :: G, Gtc, G1, G2, G3, G4, T
  real(kind=8), allocatable, dimension(:,:) :: G5, Gn_ee

  !! Define Private Variables
  real(kind=8), private :: f1, f2, f3, f4, f5, f5_o
  real(kind=8), private, dimension(4) :: norm
  integer(kind=4), private, allocatable, dimension(:) :: rcount, disp, disp_mpi, rcountv, disp_mpiv, disp_np, rcount_np
  integer(kind=4), private, allocatable, dimension(:) :: rcount_scat, disp_scat
  integer(kind=8), private, allocatable, dimension(:) :: IPIV
  complex(kind=8), private, allocatable, dimension(:) :: WORK
  complex(kind=8), private, allocatable, dimension(:) :: DU,DL,D,DU2
  real(kind=8), allocatable, dimension(:,:,:) :: Si, Si1, Sit

  real(kind=8), private, allocatable, dimension(:) :: Trt1, Trt2, I1t, I2t, I3t, Vr, f_et, NN
  real(kind=8), private, allocatable, dimension(:,:) :: At, Gnt, Gpt, Gft, Sigt, Tprt1
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Main - Allocates lapack varibles and begins main loop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Main

  integer(kind=4) :: ierr, n
  real(kind=8) ::  rn, tn, DH_inv
  character*100 :: var

  call Negf_Lapack_Alloc

  LDA = Np
  c0 = 0; c1 = 0; c2 = 0; c3 = 0
  nsb = 0; Isb = 1; Icb = 0

  !This is for splitting up the scattering and mixing routines
  !Profiling says the mixing consumes 64% of the comp time
  if(NPROC .lt. Np) then
    tn = floor(Np/DBLE(NPROC))
    disp_np(:) = 0; rcount_np(:) = 0
    !assign 1 mix to each proc
    do n = 1, NPROC
      disp_np(n) = Nint(tn)*(n-1)
      rcount_np(n) = Nint(tn)
    end do
    !remainder 
    rn = DBLE(Np) - floor(tn*DBLE(NPROC))
    do n = 1, Nint(rn)
      disp_np(n) = disp_np(n) + 1
      rcount_np(n) = rcount_np(n) + 1
    end do
    !increase remaining displacements
    disp_np(int(rn)+1:NPROC) = disp_np(int(rn)+1:NPROC) + Nint(rn)
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
  if(psc .eq. 1) call MPI_Bcast(DH,Np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(E,NE,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Ns,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(dE,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Eo,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Ef,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  dw = dE/hbar   
  call MPI_Barrier(MPI_COMM_WORLD,ierr)  !Sync

  !Do this upfront for scattering
  if(psc .eq. 1 .or. rank .eq. 0) then 
  SigInpNewt = 0.0; SigOutpNewt = 0.0
  do n=1, Np
    !unroll inv
    DH_inv = DH(n)
    !note: there are Np number of eigenvalues: DH(j), this is just a summation over all eigenvalues
    SigInpNewt  = SigInpNewt  + DH_inv !omega^2
    SigOutpNewt = SigOutpNewt + DH_inv !omega^2
  end do
  end if

  if(RANK .eq. 0) then
    call Negf_Root
    if(debug_l>5) call State_Plot ! Zero current LDOS plot
  else
    call Negf_Client
  end if

  call Negf_Lapack_Dealloc

end subroutine Negf_Main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Lapack_Alloc - Allocate Lapack space
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Lapack_Alloc

  integer(kind=4) :: err
  character*100 :: var

  allocate(Trt1(NE), Trt2(NE), I1t(NE), I2t(NE), I3t(NE), Sigt(Np,En), stat=err)
  allocate(At(Np,NE),Gnt(Np,NE),Gpt(Np,NE), Gft(Np,NE), Tprt1(Np,NE), Tprt_old(Np), stat=err)
  allocate(Vr(4), NN(Np), f_et(NE), stat=err)
  Trt1(:)=0.0; Trt2(:)=0.0; I1t(:)=0.0; I2t(:)=0.0; I3t(:)=0.0; Sigt(:,:)=0.0;
  At(:,:)=0.0; Gnt(:,:)=0.0; Gpt(:,:)=0.0; Gft(:,:) = 0.0
  Vr(:)=0.0; NN(:)=0.0; f_et(:)=0.0; Tprt1(:,:) = 0.0; Tprt_old(:) = 0.0

  allocate(WORK(LWORK), stat=err) !dlange needs work

  if(csp .eq. 1) then 
    allocate(IPIV(Np), stat=err) !Allocate LAPACK Variables
  end if

  if(csp .eq. 4) then 
    allocate(IPIV(Np), stat=err) !Allocate LAPACK Variables
    allocate(D(Np),DU(Np-1),DL(Np-1),DU2(Np-2), stat=err) !Allocate LAPACK Variables
    D(:)=0.0;DU(:)=0.0;DL(:)=0.0;DU2(:)=0.0
  end if

  if(csp .eq. 5) then 
    allocate(IPIV(Np), stat=err) !Allocate LAPACK Variables
    allocate(D(Np),DU(Np-1),DL(Np-1), stat=err) !Allocate LAPACK Variables
    D(:)=0.0;DU(:)=0.0;DL(:)=0.0;
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
  deallocate(At, Gnt, Gpt, Gft, stat=err)
  deallocate(Vr, NN, f_et, stat=err)
  deallocate(rcount, disp, disp_mpi, stat=err)
  deallocate(rcountv, disp_mpiv, stat=err)
  deallocate(rcount_np, disp_np, stat=err)
  deallocate(rcount_scat, disp_scat, stat=err)
  deallocate(WORK, stat=err)
  if(csp .eq. 1) deallocate(IPIV, stat=err)
  if(csp .eq. 4) deallocate(IPIV,D,DU,DL,DU2, stat=err)
  if(csp .eq. 5) deallocate(IPIV,D,DU,DL, stat=err)

end subroutine Negf_Lapack_Dealloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Client - Client Proc
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Client

  integer(kind=4) :: i, ierr, n

  lde = (En * rank) + 1 !determine stop value
  i = 1 
  Its = 0; Itr = 0
  f_e(:,:) = 0.0

  do while(i .ne. 0)

    call MPI_Bcast(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    if(i .ne. 0) then !messy logic
      call MPI_Bcast(Vr,4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
      if(tst .lt. 5) then
      call MPI_Scatterv(SigInp,rcount,disp_mpi,MPI_DOUBLE_PRECISION,SigInp(:,lde:(lde + En - 1)),&
                        Np*En,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_Scatterv(SigOutp,rcount,disp_mpi,MPI_DOUBLE_PRECISION,SigOutp(:,lde:(lde + En - 1)),&
                        Np*En,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      end if
      if(Vr(3) .eq. 1 .and. Vr(4) .eq. 1) call MPI_Bcast(Gn_ee,Np*Nb,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      if(i .eq. 1) then
        c2 = 1
        f_e(lde:(lde + En - 1),3) = 1e-12
        SigInp(:,:) = 0.0; SigOutp(:,:) = 0.0 
      end if

      IT(:) = 0.0; I1(:) = 0.0; I2(:) = 0.0; I3(:) = 0.0; Tr(:,:) = 0.0; 
      A(:,:) = 0.0; Gf(:,:) = 0.0; Gn(:,:) = 0.0; Gp(:,:) = 0.0
      Ds = Vr(1); V = Vr(2); c1=Vr(3); c2=Vr(4)

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
         !if(abs(f1-f2) .eq. 0 .and. c1 .ne. 1) stp = 0

         f_e(k,1) = f1; f_e(k,2) = f2; f_e(k,3) = f5 
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
      call MPI_Gatherv(f_e(lde,3),En,MPI_DOUBLE_PRECISION,f_et,rcountv,disp_mpiv,MPI_DOUBLE_PRECISION, &
                       0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(Tprt,Tprt1,Np*NE,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if(tst .lt. 5) then
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

      if(psc .eq. 1 .and. tst .lt. 5) then
      select case (tst)
        case (1) ! Simple mixing
          call Negf_Scatter_parallel !Construct scattering matricies - Inelastic
        case (2) ! Simple mixing
          call Negf_Scatter_Elast !Construct scattering matricies - Elastic
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

   integer(kind=4) :: n, ierr, i, sc_ee=1
   character*100 :: op

   i=1 !MPI End-Client Flag

   !! START OF MAIN CURRENT-VOLTAGE CALCULATIONS !! 
   if(debug_l > 0) write(6,*)
   if(debug_l > 0) write(6,*) 'START OF MAIN CURRENT-VOLTAGE CALCULATION' ! Write to screen
   write(op,'(A)') 'START OF MAIN CURRENT-VOLTAGE CALCULATION' !Write to file
   call debug_output(op)
   call flush(6)

   !Zero everything
   c0 = 0
   do while(sc_ee .gt. 0)
     if(debug_l > 5) then 
       A1(:,:) = 0.0; Gf1(:,:) = 0.0; Gp1(:,:) = 0.0
       Tr1(:,:) = 0.0; Id1(:,:) = 0.0; Id2(:,:) = 0.0; Id3(:,:) = 0.0 
     end if
     Isb = 0.0; IT(:) = 0.0; Icb = 0.0; Itot = 0.0; Qtot = 0.0
     csi = 0; f_e(:,:) = 0.0; Tprtt(:,:) = 0.0
     c0 = c0 + 1 ! Loop Counter
     c1 = 0
   do while((c1+1) .le. Nb .and. csi .lt. 2) ! Frequency Loop 

     c1 = c1 + 1 ! Loop Counter
     Ds = DH(c1) !Eigenvalues
     chng = 1; f_e(1:En,3) = 1e-12 !slightly non-zero
     if(tst .lt. 5) call Scatter_reset !reset potential variables
     SigInp(:,:) = 0.0; SigOutp(:,:) = 0.0 
     c2 = 0

     do while(chng > sc_sc) ! Scattering convergence loop

         c2 = c2 + 1 ! Loop Counter
         k = 1
         IT(:) = 0.0; I1(:) = 0.0; I2(:) = 0.0; I3(:) = 0.0; Tr(:,:) = 0.0; 
         A(:,:) = 0.0; Gf(:,:) = 0.0; Gn(:,:) = 0.0; Gp(:,:) = 0.0
         Vr(1) = Ds; Vr(2) = V; Vr(3) = c1; Vr(4) = c2; Tprt(:,:) = 0  !Setup variables to MPIsend

         ! Mpi Send
         if(c2 .eq. 1) then
            i = 1
         else
            i = 3 !do nothing
         end if
         call MPI_Bcast(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         call MPI_Bcast(Vr,4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         if(tst .lt. 5) then
         call MPI_Scatterv(SigInp,rcount,disp_mpi,MPI_DOUBLE_PRECISION,Sigt, &
                           Np*En,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         call MPI_Scatterv(SigOutp,rcount,disp_mpi,MPI_DOUBLE_PRECISION,Sigt, &
                           Np*En,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         end if
         if(c1 .eq. 1 .and. c2 .eq. 1) call MPI_Bcast(Gn_ee,Np*Nb,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
         call MPI_Barrier(MPI_COMM_WORLD,ierr)

         do while(k .le. En) ! Starting energy loop for incoming electron energy
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

           f_e(k,1) = f1; f_e(k,2) = f2; f_e(k,3) = f5
           k = k + 1 !increment loop
           cfr = 1 !first run

         end do !Energy loop

         ! Calculate energy refinement for AMR
         if(amr_o .ne. 0) then
           call energy_amr_refine_log(IT,Er)
         end if

         Its = sum(IT);
         Itr = 0; Trt1(:) = 0; Trt2(:) = 0; I1t(:) = 0; I2t(:) = 0; I3t(:) = 0; At(:,:) = 0
         Gft(:,:) = 0; Gpt(:,:) = 0; Gnt(:,:) = 0; f_et(:) = 0.0; NN=0; Tprt1(:,:) = 0.0
         if(NPROC .ne. 1) then
            call MPI_Reduce(Rh,NN,Np,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_Reduce(Its,Itr,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            call MPI_Gatherv(f_e,En,MPI_DOUBLE_PRECISION,f_et,rcountv,disp_mpiv,MPI_DOUBLE_PRECISION, &
                             0,MPI_COMM_WORLD,ierr)
            call MPI_Reduce(Tprt,Tprt1,Np*NE,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            if(tst .lt. 5) then
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
            f_e(:,3) = f_et
         else
            Itr = Its; Gnt = Gn; Gpt = Gp; f_et = f_e(:,3) 
            if(debug_l>5) then
              Trt1 = Tr(:,1); Trt2 = Tr(:,2); I1t = I1; I2t = I2; I3t = I3; At = A; Gft=Gf; NN = Rh
            end if
         end if

         select case (tst)
           case (1) ! Simple mixing
             if(psc .eq. 0) then
               call Negf_Scatter !Construct scattering matricies - Inelastic
             else
               call Negf_Scatter_parallel !Construct scattering matricies - Inelastic
             end if
           case (2) ! Simple mixing
             call Negf_Scatter_Elast !Construct scattering matricies - Elastic
           case (3) ! Simple mixing
             call Negf_Scatter_Single !Construct scattering matricies - Inelastic, single phonon
           case (4) ! Simple mixing
             call Negf_Scatter_Elast_Single !Construct scattering matricies - Inelastic, single phonon
           case default
             chng = 0.0 !no scattering
         end select
         call MPI_Barrier(MPI_COMM_WORLD,ierr) !Align Processes
       end do ! End of scattering loop 

       if(debug_l>5) then
         do n = 1, NE !need to divide by fermi function to get scattering rate in 1/s
           Si1(:,n,:) = Si1(:,n,:) + Si(:,n,:)/f_et(n)
         end do
         A1 = A1 + At !Local density of states
         Gf1 = Gf1 + Gft !Electron density
         Gp1 = Gp1 + (At-Gft) !Local empty states
         Aall(:,:,c1) = At; Gnall(:,:,c1) = Gft
         Tr1(:,c1) = Trt1; !Tr2(:,c1) = Trt2 ! Transmission for each subband
         Id1(:,c1) = Io*I1t*dw; Id2(:,c1) = Io*I2t*dw; Id3(:,c1) = Io*I3t*dw; !Total current/subband
         Rh1(:,c1) = NN
       end if

       Tprtt(:,:) = Tprtt(:,:) + Tprt1(:,:)
       Isb = Io*Itr*dw !Subband current
       Itot = Itot + Isb !Total current
       Qtot = Qtot + sum(At(:,c1))*E(c1)*dw
     
       if(debug_l>2) write(6,'(2X,A,ES12.5,1X,ES12.5,1X,ES12.5,1X,ES12.5)') 'w, J, Jtot, Qtot = ',sqrt(Ds),Isb,Itot,Qtot

     if(c0 .eq. 1) call Loop_Plot
     call flush(6)

     !if(Isb/Itot .lt. sc_dv) csi = csi + 1 !stop once contribution is below threshold

   end do ! End of subband loop

     Tprt1 = Tprtt !reuse array
     call Negf_Temperature !Calculate Average Temperature
     call Output_Tprt(Tprtt,Rh1,c0)
     call Negf_Temperature_Norm !Calulate Norm of Temperature
     Tprt_old(:) = Tprtt(:,1) !save old temperature

     write(*,*) 'Computed Thermal Conductivity = ',abs(Itot)*Lt/dT,' W/K-m'
     write(*,*) 'Computed Thermal Conductance = ',abs(Itot)/dT,' W/K-m^2'

     if(c0 .eq. 1) then
       call Mpinterface_Write_Eig(sqrt(DH*hbar_eV**2),Nb,debug_l)
       call Mpinterface_Write_Fermi(f_e(:,3),NE,debug_l)
     end if
     call Mpinterface_Write_Gn(Gnall,NE,Eo,Ef,Np,Nb,debug_l,chng_tprt)
     call Mpinterface_Write_Gp((Aall-Gnall),NE,Np,Nb,debug_l,chng_tprt)
     call Mpinterface_Read_Scph(Gn_ee,Np,Nb,sc_ee,debug_l)

   end do !phonon scatter couple

   i = 0 ! End client processes
   call MPI_Bcast(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_Barrier(MPI_COMM_WORLD,ierr)
   if(debug_l>0) write(6,*) 'END OF MAIN CURRENT-VOLTAGE CALCULATION' !Write to screen
   write(op,'(A)') 'END OF MAIN CURRENT-VOLTAGE CALCULATION' !Write to file
   call debug_output(op)
   if(debug_l>2) write(6,*) 'Number of used(calculated) phonon modes:',c1
   if(debug_l>2) write(6,*) 'Total number of phonon modes:',Nb !Write to file
   call flush(6)

   !Tprt1 = Tprtt !reuse array
   !call Negf_Temperature !Calculate Average Temperature

end subroutine Negf_Root

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Temperature_Norm - Calculate the norm of average temperature along device
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Temperature_Norm

   double precision dlange
   external dlange

   real(kind=8) :: norm1, norm2

   norm1 = dlange('i',Np,1,(Tprtt(:,1)-Tprt_old(:)),LDA,WORK)
   norm2 = dlange('i',Np,1,Tprt_old(:),LDA,WORK) 
write(*,*) 'norm1, norm2 =',norm1,norm2
   if(norm2 .ne. 0.0) then
     chng_tprt = norm1/norm2
   else
     chng_tprt = 0.0
   end if

end subroutine Negf_Temperature_Norm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Temperature - Calculate average temperature along device
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Temperature

  integer(kind=4) :: k, n
  real(kind=8) :: f1,m,m2,sum_DOS,sum_E,sum_f1
  real(kind=8) :: Tp_t


   Tprtt(:,:) = 0.0 !zero array
   do n = 1, Np !Tprt1 is a density
     sum_E=0;sum_DOS=0;sum_f1=0
     do k = 1, NE
       call Startup_DOS(vsi,E(k),kT1,No1)
       call Startup_DOS(vsi,E(k),kT2,No2)
       !Derived the fermi function for Bose-Einstein particles, not the same as electronic NEGF
       f1 = -No1*dlog(1-dexp((-E(k))/kT1)) ! 2D Source fermi function. Units of 1/m^2 
       f2 = -No2*dlog(1-dexp((-E(k))/kT2)) ! 2D drain fermi function. Units of 1/m^2
       m = (f2-f1)/Np
       m2 = (No2-No1)/Np
       sum_DOS = sum_DOS + No1
       sum_E = sum_E + (f1 + Tprt1(n,k))*E(k)
       sum_f1 = sum_f1 + f1

       !this routine calculates the average temperature at every energy level
       Tp_t = -E(k)/(kBm*dlog(1-dexp(-((m*(n-1)+f1)+sum(Tprt1(n,:)))/(m2*(n-1)+No1))))
       Tprt_NE(n,k) = Tp_t
       Tprtt(n,1) = Tprtt(n,1) - E(k)/(kBm*dlog(1-dexp(-((m*(n-1)+f1)+sum(Tprt1(n,:)))/(m2*(n-1)+No1))))/NE
       if(Tp_t .gt. Tprtt(n,2)) Tprtt(n,2) = Tp_t !max temperature
       !Tprtt(n,2) = Tprtt(n,2) - E(k)/(kBm*dlog(1-dexp(-(sum(Tprt1(n,:)))/(m2*(n-1)+No1))))/NE
     end do
!write(*,*) 'sum_E,sum_DOS =',sum_E,sum_DOS,sum_f1

     call Negf_Temperature_Search(sum_E,Tprtt(n,3),n)
   end do

   !do n = 1, Np
   !   Tprtt(n,2) = sum_E/(kBm*dlog(1-dexp(-((m*(n-1)+f1)+sum(Tprt1(n,:)))/(sum_DOS))))
   !end do 

end subroutine Negf_Temperature
  

subroutine Negf_Debye(E_db,T)

  real(kind=8), intent(out) :: E_db
  real(kind=8), intent(in) :: T
  integer(kind=4) :: n,ndo
  real(kind=8) :: D_db, dx, x1, vs

  E_db = 0; D_db = 0
  dx = dE/(hbar)

  do n = 1, NE
     !2d Dos
     call Startup_DOS(vsi,E(n),kBm*T,No1)
     !Derived the fermi function for Bose-Einstein particles, not the same as electronic NEGF
     f1 = -No1*dlog(1-dexp((-E(n))/(kBm*T))) ! 2D Source fermi function. Units of 1/m^2 
     E_db = E_db + f1*E(n)

     !3d Dos
     !x1 = E(n)/(hbar)
     !E_db = E_db + 3*kbm**4*T**4/(2*pi**2*vs**2*hbar**3)*x1**3/(exp(x1)-1)*dx
     !E_db = E_db + 3*hbar/(2*pi**2*vs**3)*x1**3/(exp((x1*hbar)/(kbm*T))-1)*dx
     !E_db = E_db + hbar*x1/(exp((x1*hbar)/(kbm*T))-1)*dx
     !D_db = x1**2/(2*pi**2*vs)
  end do

  !E_db = E_db !/D_db

end subroutine Negf_Debye

subroutine Negf_Temperature_Search(E_pp,T,locn)

  integer(kind=4), intent(in) :: locn
  real(kind=8), intent(in) :: E_pp
  real(kind=8), intent(out) :: T
  integer(kind=4) :: nlp
  real(kind=8) :: Eq, Eq_old, Tq, Tq_old, cntrl, dTq
  real(kind=8) :: Tq_1, Eq_1, Tq_2, Eq_2
  !quess at temperature
  Tq = 100; dTq = 1000; Tq_old = Tq
  ! calculate energy at temperature
  call Negf_Debye(Eq,Tq)
  !call Negf_Debye(Eq_old,Tq_old)

  cntrl = 1; nlp = 1
  Eq_1 = Eq; Eq_2 = Eq+dTq
  Tq_1 = Tq; Tq_2 = Tq+dTq

  do while (cntrl > 1e-10)
    write(*,*) 'Trying to guess Tq,E_pp,Eq',Tq,E_pp,Eq
    !if(nlp .ne. 1) Tq_1 = Tq
    !if(nlp .ne. 1) Eq_1 = Eq
    if(Eq_1 .lt. E_pp .and. Eq_2 .gt. E_pp) then
        Tq = Tq_1+(Tq_2-Tq_1)/(Eq_2-Eq_1)*(Eq_2-Eq_1)/2
    else if(Eq_1 .gt. E_pp .and. Eq_2 .lt. E_pp) then
        Tq = Tq_1+(Tq_2-Tq_1)/(Eq_2-Eq_1)*(Eq_2-Eq_1)/2
    else if(Eq_2 .gt. E_pp) then
        Tq = Tq_1+(Tq_2-Tq_1)/(Eq_2-Eq_1)*(Eq_2-Eq_1)*3/2
    else if(Eq_2 .lt. E_pp) then
        Tq = Tq_1+(Tq_2-Tq_1)/(Eq_2-Eq_1)*(Eq_2-Eq_1)*3/2
    end if
    Tq_1 = Tq_2; Eq_1 = Eq_2
    Tq_2 = Tq
    call Negf_Debye(Eq,Tq)
    Eq_2 = Eq
    !call Negf_Debye(Eq,Tq)
    if(nlp .gt. 1) cntrl = abs((Eq-E_pp)/E_pp)
    if(abs(Tq_old-Tq) .lt. 1e-12) cntrl = 0
write(*,*) 'cntrl = ',cntrl,Eq,E_pp
    nlp = nlp + 1;
  end do

  T = Tq

end subroutine Negf_Temperature_Search

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Scatter - Calculate scattering matricies - Inelastic Scattering
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Scatter

    double precision dlange
    external dlange

    integer(4) :: off, n, j, fp, fe
    real(8) :: NW, NE_off, w_ph

    SigInpNew(:,:) = 0.0; SigOutpNew(:,:) = 0.0

    fp = floor(sqrt(DH(1))/dw)-Ns !starting eigenvalue, note first bin not at zero so subtract Ns
    fe = ceiling(sqrt(DH(Np))/dw)-Ns

    do n=1, Np
      w_ph = DH(n) !omega^2
      if(mwa .eq. 1) then !use Maxwellian Approx for Dist
        Nw = exp(-sqrt(DH(n))*hbar/kT) ! Simplified Bose-Einstein distribution
      else
        Nw = 1/(exp(sqrt(DH(n))*hbar/kT)-1) !true BE dist
      end if
      off = ceiling(sqrt(DH(n))/dw)-Ns

      !! Reassigning density of states to accomodate scattering
      !! Need to take the view point of the final electron state 
      !nse = eoshift(Gnt, shift = -off, dim = 2) !E+hw emit phonon,  !1/omega^2-m^2
      !nsa = eoshift(Gnt, shift = off, dim = 2) !E-hw absorb phonon, !1/omega^2-m^2
      !pse = eoshift(Gpt, shift = -off, dim = 2) !E+hw emit phonon,  !1/omega^2-m^2
      !psa = eoshift(Gpt, shift = off, dim = 2) !E-hw absorb phonon, !1/omega^2-m^2
      call Matxopt_eoshift(Gnt,NE,off,nse) !E-hw final state is after emission
      call Matxopt_eoshift(Gnt,NE,-off,nsa)  !E+hw final state is after absorption
      call Matxopt_eoshift(Gpt,NE,off,pse) !E-hw
      call Matxopt_eoshift(Gpt,NE,-off,psa)  !E+hw

      ! need to zero matrix below lowest eigenvalue because +hw can't result
      ! in phonon scattered below allowed modes or above
      call Matxopt_eozero(nsa,NE,fe)
      call Matxopt_eozero(psa,NE,fe)
      call Matxopt_eozero(nse,NE,fe)
      call Matxopt_eozero(pse,NE,fe)

      call Matxopt_eozero(nsa,NE,-fp)
      call Matxopt_eozero(psa,NE,-fp) 
      call Matxopt_eozero(nse,NE,-fp)
      call Matxopt_eozero(pse,NE,-fp)

      !unroll multiplications
      !use mattherson rule to add scattering rates
      !1/omega^2-m^2*omega^4 = omega^2/m^2
      !SigInpNewt and SigOutpNewt Calculated at beginning at top
      SigInpNew  = SigInpNew + So/(2*pi)*w_ph*SigInpNewt*((1+Nw)*nse + Nw*nsa)/Np
      SigOutpNew =  SigOutpNew + So/(2*pi)*w_ph*SigInpNewt*((1+Nw)*(nse+psa) + Nw*(nsa+pse))/Np

      if(debug_l >5) Si(:,:,n) = sqrt(So/(2*pi)*w_ph*SigInpNewt*((1+Nw)*(psa+nse) + Nw*(nsa+pse))/Np) !scattering rate 1/s-m^2

    end do

    ! Poisson convergence method
    select case (coc_s)
      case (1) ! Simple mixing
        call Scatter_lin_smix(SigInpNew, SigInp)
        call Scatter_lin_smix(SigOutpNew, SigOutp)
      case (2) ! Anderson mixing
        call Scatter_anderson(SigInpNew,SigOutpNew,SigInp,SigOutp)
    end select

    norm(1) = dlange('i',Np,NE,(SigInpNew-SigInp),LDA,WORK)
    norm(2) = dlange('i',Np,NE,SigInp,LDA,WORK)    
    norm(3) = dlange('i',Np,NE,(SigOutpNew-SigOutp),LDA,WORK)
    norm(4) = dlange('i',Np,NE,SigOutp,LDA,WORK)

    if(norm(2) .ne. 0.0 .and. norm(4) .ne. 0.0) then
      chng = norm(1)/norm(2) + norm(3)/norm(4)
    else
      if(norm(1) .ne. 0.0 .and. norm(3) .ne. 0.0) then
        chng = 1;
      else
        chng = 0;
      end if 
    end if

    if(debug_l>2) write(6,'(/,A,ES11.4)') 'Scatter chg = ',chng ! Write change of charge to screen

end subroutine Negf_Scatter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Scatter_Single - Calculate scattering matricies - Inelastic Scattering
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Scatter_Single

    double precision dlange
    external dlange

    integer(4) :: off, n, j, ps, fp, fe, ur, lr
    real(8) :: NW, NE_off, w_ph, w_inv

    SigInpNew(:,:) = 0.0; SigOutpNew(:,:) = 0.0

      fp = floor(sqrt(DH(1))/dw)-Ns
      fe = ceiling(sqrt(DH(Np))/dw)-Ns

      ps = Np/2 !Select phonon in middle of BZ
      w_ph = DH(ps) !omega^2
      w_inv = DH(ps)
      if(mwa .eq. 1) then !use Maxwellian Approx for Dist
        Nw = exp(-sqrt(DH(ps))*hbar/kT) ! Simplified Bose-Einstein distribution
      else
        Nw = 1/(exp(sqrt(DH(ps))*hbar/kT)-1) !true BE dist
      end if
      off = ceiling(sqrt(DH(ps))/dw)-Ns

      !! Reassigning density of states to accomodate scattering
      !! Need to take the view point of the final electron state 
      !nse = eoshift(Gnt, shift = -off, dim = 2) !E+hw emit phonon,  !1/omega^2-m^2
      !nsa = eoshift(Gnt, shift = off, dim = 2) !E-hw absorb phonon, !1/omega^2-m^2
      !pse = eoshift(Gpt, shift = -off, dim = 2) !E+hw emit phonon,  !1/omega^2-m^2
      !psa = eoshift(Gpt, shift = off, dim = 2) !E-hw absorb phonon, !1/omega^2-m^2
      call Matxopt_eoshift(Gnt,NE,off,nse) !E-hw final state is after emission
      call Matxopt_eoshift(Gnt,NE,-off,nsa)  !E+hw final state is after absorption
      call Matxopt_eoshift(Gpt,NE,off,pse) !E-hw
      call Matxopt_eoshift(Gpt,NE,-off,psa)  !E+hw

      ! need to zero matrix below lowest eigenvalue because +hw can't result
      ! in phonon scattered below allowed modes
      call Matxopt_eozero(nsa,NE,fe)
      call Matxopt_eozero(psa,NE,fe)
      call Matxopt_eozero(nse,NE,fe)
      call Matxopt_eozero(pse,NE,fe)

      call Matxopt_eozero(nsa,NE,-fp)
      call Matxopt_eozero(psa,NE,-fp) 
      call Matxopt_eozero(nse,NE,-fp)
      call Matxopt_eozero(pse,NE,-fp)

      !unroll multiplications
      !use mattherson rule to add scattering rates
      !1/omega^2-m^2*omega^4 = omega^2/m^2
      !SigInpNewt and SigOutpNewt Calculated at beginning at top
      SigInpNew  = SigInpNew + So/(2*pi)*w_ph*SigInpNewt*((1+Nw)*nse + Nw*nsa)
      SigOutpNew =  SigOutpNew + So/(2*pi)*w_ph*SigInpNewt*((1+Nw)*(nse+psa) + Nw*(nsa+pse))

      if(debug_l >5) Si(:,:,1) = sqrt(So/(2*pi)*w_ph*SigInpNewt*((1+Nw)*(psa+nse) + Nw*(nsa+pse))) !scattering rate 1/s-m

    ! Poisson convergence method
    select case (coc_s)
      case (1) ! Simple mixing
        call Scatter_lin_smix(SigInpNew, SigInp)
        call Scatter_lin_smix(SigOutpNew, SigOutp)
      case (2) ! Anderson mixing
        call Scatter_anderson(SigInpNew,SigOutpNew,SigInp,SigOutp)
    end select

    norm(1) = dlange('i',Np,NE,(SigInpNew-SigInp),LDA,WORK)
    norm(2) = dlange('i',Np,NE,SigInp,LDA,WORK)    
    norm(3) = dlange('i',Np,NE,(SigOutpNew-SigOutp),LDA,WORK)
    norm(4) = dlange('i',Np,NE,SigOutp,LDA,WORK)

    if(norm(2) .ne. 0.0 .and. norm(4) .ne. 0.0) then
      chng = norm(1)/norm(2) + norm(3)/norm(4)
    else
      if(norm(1) .ne. 0.0 .and. norm(3) .ne. 0.0) then
        chng = 1;
      else
        chng = 0;
      end if 
    end if

    if(debug_l>2) write(6,'(/,A,ES11.4)') 'Scatter chg = ',chng ! Write change of charge to screen

end subroutine Negf_Scatter_Single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Scatter - Calculate scattering matricies - Inelastic Scattering
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Scatter_Elast_Single

    double precision dlange
    external dlange

    integer(4) :: off, n, j, ps
    real(8) :: NW, NE_off, w_ph, w_inv

    SigInpNew(:,:) = 0.0; SigOutpNew(:,:) = 0.0

      ps = Np/2 !Select phonon in middle of BZ
      w_ph = DH(ps) !omega^2
      w_inv = DH(ps)
      if(mwa .eq. 1) then !use Maxwellian Approx for Dist
        Nw = exp(-sqrt(DH(ps))*hbar/kT) ! Simplified Bose-Einstein distribution
      else
        Nw = 1/(exp(sqrt(DH(ps))*hbar/kT)-1) !true BE dist
      end if
      off = ceiling(sqrt(DH(ps))/dw)-Ns

      !! Reassigning density of states to accomodate scattering
      !! Need to take the view point of the final electron state 
      nse = Gnt !eoshift(Gnt, shift = -off, dim = 2) !E+hw emit phonon,  !1/omega^2-m^2
      nsa = Gnt !eoshift(Gnt, shift = off, dim = 2) !E-hw absorb phonon, !1/omega^2-m^2
      pse = Gpt !eoshift(Gpt, shift = -off, dim = 2) !E+hw emit phonon,  !1/omega^2-m^2
      psa = Gpt !eoshift(Gpt, shift = off, dim = 2) !E-hw absorb phonon, !1/omega^2-m^2

      !unroll multiplications
      !use mattherson rule to add scattering rates
      !1/omega^2-m^2*omega^4 = omega^2/m^2
      !SigInpNewt and SigOutpNewt Calculated at beginning at top
      SigInpNew  = SigInpNew + So/(2*pi)*w_ph*SigInpNewt*((1+Nw)*nse + Nw*nsa)
      SigOutpNew =  SigOutpNew + So/(2*pi)*w_ph*SigInpNewt*((1+Nw)*(nse+psa) + Nw*(nsa+pse))

      if(debug_l >5) Si(:,:,1) = sqrt(So/(2*pi)*w_ph*SigInpNewt*((1+Nw)*(psa+nse) + Nw*(nsa+pse))) !scattering rate 1/s-m^2

    ! Poisson convergence method
    select case (coc_s)
      case (1) ! Simple mixing
        call Scatter_lin_smix(SigInpNew, SigInp)
        call Scatter_lin_smix(SigOutpNew, SigOutp)
      case (2) ! Anderson mixing
        call Scatter_anderson(SigInpNew,SigOutpNew,SigInp,SigOutp)
    end select

    norm(1) = dlange('i',Np,NE,(SigInpNew-SigInp),LDA,WORK)
    norm(2) = dlange('i',Np,NE,SigInp,LDA,WORK)    
    norm(3) = dlange('i',Np,NE,(SigOutpNew-SigOutp),LDA,WORK)
    norm(4) = dlange('i',Np,NE,SigOutp,LDA,WORK)

    if(norm(2) .ne. 0.0 .and. norm(4) .ne. 0.0) then
      chng = norm(1)/norm(2) + norm(3)/norm(4)
    else
      if(norm(1) .ne. 0.0 .and. norm(3) .ne. 0.0) then
        chng = 1;
      else
        chng = 0;
      end if 
    end if

    if(debug_l>2) write(6,'(/,A,ES11.4)') 'Scatter chg = ',chng ! Write change of charge to screen

end subroutine Negf_Scatter_Elast_Single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Scatter_parallel - Calculate scattering matricies in parallel - Inelastic Scattering
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Scatter_parallel

    double precision dlange
    external dlange

    integer(4) :: off, n, j, ldc, lds, ierr, fp, fe
    real(8) :: NW, NE_off, w_ph

    SigInpNew(:,:) = 0.0; SigOutpNew(:,:) = 0.0;

    fp = floor(sqrt(DH(1))/dw)-Ns
    fe = ceiling(sqrt(DH(Np))/dw)-Ns

    ldc = disp_np(rank+1)
    lds = ldc + rcount_np(rank+1)

    !Need to broadcast matricies to all processors
    call MPI_Bcast(Gnt,Np*NE,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(Gpt,Np*NE,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    ! this loop will split the calculation arcross multiple proc
    if(ldc .ge. 0) then
    do n=ldc+1, lds
      w_ph = DH(n) !omega^2
      if(mwa .eq. 1) then !use Maxwellian Approx for Dist
        Nw = exp(-sqrt(DH(n))*hbar/kT) ! Simplified Bose-Einstein distribution
      else
        Nw = 1/(exp(sqrt(DH(n))*hbar/kT)-1) !true BE dist
      end if
      off = ceiling(sqrt(DH(n))/dw)-Ns

      !! Reassigning density of states to accomodate scattering
      !! Need to take the view point of the final electron state 
      !nse = eoshift(Gnt, shift = -off, dim = 2) !E+hw emit phonon,  !1/omega^2-m^2
      !nsa = eoshift(Gnt, shift = off, dim = 2) !E-hw absorb phonon, !1/omega^2-m^2
      !pse = eoshift(Gpt, shift = -off, dim = 2) !E+hw emit phonon,  !1/omega^2-m^2
      !psa = eoshift(Gpt, shift = off, dim = 2) !E-hw absorb phonon, !1/omega^2-m^2
      !These are custom eoshift routines that are more efficient
      call Matxopt_eoshift(Gnt,NE,off,nse) !E-hw final state is after emission
      call Matxopt_eoshift(Gnt,NE,-off,nsa)  !E+hw final state is after absorption
      call Matxopt_eoshift(Gpt,NE,off,pse) !E-hw
      call Matxopt_eoshift(Gpt,NE,-off,psa)  !E+hw

      ! need to zero matrix below lowest eigenvalue because +hw can't result
      ! in phonon scattered below allowed modes
      call Matxopt_eozero(nsa,NE,fe)
      call Matxopt_eozero(psa,NE,fe)
      call Matxopt_eozero(nse,NE,fe)
      call Matxopt_eozero(pse,NE,fe)

      call Matxopt_eozero(nsa,NE,-fp)
      call Matxopt_eozero(psa,NE,-fp) 
      call Matxopt_eozero(nse,NE,-fp)
      call Matxopt_eozero(pse,NE,-fp)

      !moved outside loop, this a constant.
      ! This loop scatters with all the other modes present in equilibrium dist
      ! there are Np phonon modes
      !SigInpNewt = 0.0; SigOutpNewt = 0.0
      !do j=1, Np
      !  !unroll inv
      !  DH_inv = 1/DH(j)
      !  !note: there are Np number of eigenvalues: DH(j), this is just a summation over all eigenvalues
      !  SigInpNewt  = SigInpNewt  + DH_inv !1/omega^2-m^2*omega^4 = omega/m^2
      !  SigOutpNewt = SigOutpNewt + DH_inv !omega^2/m^2
      !end do

      !unroll multiplications
      !use mattherson rule to add scattering rates
      !1/omega^2-m^2*omega^4 = omega/m^2
      !SigInpNewt and SigOutpNewt Calculated at beginning at top 1/omega^2
      SigInpNew  = SigInpNew + So/(2*pi)*w_ph*SigInpNewt*((1+Nw)*nse + Nw*nsa)/Np
      SigOutpNew =  SigOutpNew + So/(2*pi)*w_ph*SigInpNewt*((1+Nw)*(nse+psa) + Nw*(nsa+pse))/Np

      if(debug_l >5) Si(:,:,n) = sqrt(So/(2*pi)*w_ph*SigInpNewt*((1+Nw)*(psa+nse) + Nw*(nsa+pse))/Np) !scattering rate 1/s-m^2
    end do
    end if

    Gnt(:,:) = 0.0; Gpt(:,:) = 0.0
    !This will add up all SigInp and SigOutp across all processors
    call MPI_Reduce(SigInpNew,Gnt,Np*NE,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_Reduce(SigOutpNew,Gpt,Np*NE,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    if(rank .eq. 0) then
      !Gather scattering matrix for plotting scattering rate. This contains scatter per phonon freq.
      if(debug_l >5) then
        Sit(:,:,1:rcount_np(1)) = Si(:,:,1:rcount_np(1)) !gatherv source and destination can't be same
                                                         !so have to use temp variable. This can get big!
        !expensive routine, sending 3d array to root
        call MPI_Gatherv(Sit,Np*NE*rcount_np(rank+1),MPI_DOUBLE_PRECISION,Si,rcount_scat,disp_scat, &
                             MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      end if

      !Copy back to SigInpNew, MPI won't send and receive on the same variables so have to use temp var
      SigInpNew = Gnt
      SigOutpNew = Gpt

      ! Poisson convergence method
      select case (coc_s)
        case (1) ! Simple mixing
          call Scatter_lin_smix(SigInpNew, SigInp)
          call Scatter_lin_smix(SigOutpNew, SigOutp)
        case (2) ! Anderson mixing
          call Scatter_anderson(SigInpNew,SigOutpNew,SigInp,SigOutp)
      end select

      !compute norm of scattering matricies, this routine is faster then the mix
      norm(1) = dlange('i',Np,NE,(SigInpNew-SigInp),LDA,WORK)
      norm(2) = dlange('i',Np,NE,SigInp,LDA,WORK)    
      norm(3) = dlange('i',Np,NE,(SigOutpNew-SigOutp),LDA,WORK)
      norm(4) = dlange('i',Np,NE,SigOutp,LDA,WORK)

      if(norm(2) .ne. 0.0 .and. norm(4) .ne. 0.0) then
        chng = norm(1)/norm(2) + norm(3)/norm(4)
      else
        if(norm(1) .ne. 0.0 .and. norm(3) .ne. 0.0) then
          chng = 1;
        else
          chng = 0;
        end if 
      end if

      if(debug_l>2) write(6,'(/,A,ES11.4)') 'Scatter chg = ',chng ! Write change of charge to screen
    else
      !send scattering rates to root for output
      if(debug_l > 5) then
        call MPI_Gatherv(Si,Np*NE*rcount_np(rank+1),MPI_DOUBLE_PRECISION,Si,rcount_scat,disp_scat, &
                             MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      end if
    end if

end subroutine Negf_Scatter_parallel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Scatter_Elast - Calculate scattering matricies - Elastic Scattering
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Scatter_Elast

    double precision dlange
    external dlange

    integer(4) :: off, n, j
    real(8) :: NW, NE_off, w_ph

    SigInpNew(:,:) = 0.0; SigOutpNew(:,:) = 0.0

    do n=1, Np
      w_ph = DH(n) !omega^2
      if(mwa .eq. 1) then !use Maxwellian Approx for Dist
        Nw = exp(-sqrt(DH(n))*hbar/kT) ! Simplified Bose-Einstein distribution
      else
        Nw = 1/(exp(sqrt(DH(n))*hbar/kT)-1) !true BE dist
      end if
      off = 0 !Nint(sqrt(DH(n))/dw)

      !! Reassigning density of states to accomodate scattering
      !! Need to take the view point of the final electron state
      nse = Gnt !eoshift(Gnt, shift = -off, dim = 2) !E+hw emit phonon, !1/omega^2-m^2
      nsa = Gnt !eoshift(Gnt, shift = off, dim = 2) !E-hw absorb phonon,!1/omega^2-m^2
      pse = Gpt !eoshift(Gpt, shift = -off, dim = 2) !E+hw emit phonon, !1/omega^2-m^2
      psa = Gpt !eoshift(Gpt, shift = off, dim = 2) !E-hw absorb phonon,!1/omega^2-m^2
        
      !unroll multiplications
      !use mattherson rule to add scattering rates
      !1/omega^2-m^2*omega^4 = omega/m^2
      SigInpNew  = SigInpNew + So/(2*pi)*w_ph*SigInpNewt*((1+Nw)*nse + Nw*nsa)/Np
      SigOutpNew =  SigOutpNew + So/(2*pi)*w_ph*SigInpNewt*((1+Nw)*(nse+psa) + Nw*(nsa+pse))/Np

      if(debug_l >5) Si(:,:,n) = sqrt(So/(2*pi)*w_ph*SigInpNewt*((1+Nw)*(psa+nse) + Nw*(nsa+pse))/Np) !scattering rate 1/s-m^2
    end do

    ! Poisson convergence method
    select case (coc_s)
      case (1) ! Simple mixing
        call Scatter_lin_smix(SigInpNew, SigInp)
        call Scatter_lin_smix(SigOutpNew, SigOutp)
      case (2) ! Anderson mixing
        call Scatter_anderson(SigInpNew,SigOutpNew,SigInp,SigOutp)
    end select

    norm(1) = dlange('i',Np,NE,(SigInpNew-SigInp),LDA,WORK)
    norm(2) = dlange('i',Np,NE,SigInp,LDA,WORK)    
    norm(3) = dlange('i',Np,NE,(SigOutpNew-SigOutp),LDA,WORK)
    norm(4) = dlange('i',Np,NE,SigOutp,LDA,WORK)

    if(norm(2) .ne. 0.0 .and. norm(4) .ne. 0.0) then
      chng = norm(1)/norm(2) + norm(3)/norm(4)
    else
      if(norm(1) .ne. 0.0 .and. norm(3) .ne. 0.0) then
        chng = 1;
      else
        chng = 0;
      end if 
    end if

    if(debug_l>2) write(6,'(/,A,ES11.4)') 'Scatter chg = ',chng ! Write change of charge to screen

end subroutine Negf_Scatter_Elast

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Matmul - Matrix multiplication with matrix algebra simplification
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Matmul
   
   integer(kind=4) :: n, r
   real(kind=8) :: G1_i, G2_i, G3_i, G4_i, Rt

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

     I1(k) = I1(k) + G3_i*(f1-f2) !omega^2/m^2
     I2(k) = I2(k) + G4_i*(f2-f1) !omega^2/m^2
     Tr(k,1) = Tr(k,1) + G3_i
     Tr(k,2) = Tr(k,2) + G4_i
     G2(:,:) = 0.0
     do n = 1, Np
        if(f5 .ne. 0) G2(n,n) = SigInp(n,k)/f5 !determine gamma_in_phonon !omega^2
     end do

     !Butticker Probe 
     G1(:,:) = (0,0);
     G1(1,1) =  Gam1(1,1)
     !G3 = matmul(matmul(G2,G),matmul(G1,Gtc)) !FIXME: so faster multiply T31
     call Matxopt_zdff(G2,G,Np,G3)
     call Matxopt_zsff(G1,1,Gtc,Np,G4)
     call Matxopt_zffd2(G3,G4,Np,G5)
     call Matxopt_tracesum(G5,Np,G3_i)

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
     if(G3_i+G4_i .ne. 0) then
       f5 = (f1*G3_i+f2*G4_i-(I1(k)+I2(k)))/(G3_i+G4_i)
     end if

     I3(k) = 0.0
     !back calculate current at scattering contact
     !conservation of energy
     if(G3_i+G4_i .ne. 0.0) then
        I3(k) = (-G3_i*f1 - G4_i*f2 + f5*(G3_i+G4_i))
     end if

     I1(k) = I1(k)*E(k)
     I2(k) = I2(k)*E(k) 
     I3(k) = I3(k)*E(k) 
     IT(k) = I1(k) !FIXME: only need to do at end of k loop

     !Gn
     G1(:,:) = (0,0);
     G1(1,1) =  Gam1(1,1)*f1 !omega^2/m^2
     G1(Np,Np) = Gam2(Np,Np)*f2 !omega^2/m^2
     !Sigoutp from phonon scattering already contains fermi function
     do n = 1, Np
       G1(n,n) = G1(n,n) + SigInp(n,k)*f5/f5_o !(omega^2/m^2) simplify multiplication full matrix times diag 
     end do

     !G2 = matmul(matmul(G,G1),Gtc) !FIXME: so faster multiply
     call Matxopt_zfdf(G,G1,Np,G3) ! 1/omega^2 * omega^2/m^2
     call Matxopt_zffd(G3,Gtc,Np,G2) !1/m^2 * 1/omega^2 = 1/omega^2-m^2

     !Gn = correlation function
     do n = 1, Np
       Gn(n,k) = DBLE(G2(n,n)) !1/omega^2-m^2
     end do

     do n = 1, Np
       A(n,k) = DBLE(i*(G(n,n) - Gtc(n,n))) !*(f1+f2+f5_o) ! Channel avaliable local density of states matrix (1/omega^2)
       Gf(n,k) = Gn(n,k)/(f1+f2+f5_o) ! Channel density of states matrix (1/omega^2)
     end do

     do n = 1, Np !
       !Gp(n,k) = abs(A(n,k) - Gn(n,k)) ! Unfilled states, take view point of contact 1 !1/omega^2-m^2
       Gp(n,k) = A(n,k)
       Rt = ((dE/hbar)**2/(2.0*pi))*Gn(n,k)
       Rh(n) = Rh(n) + Rt ! Calculating electron density in the channel (1/m^2)
       Tprt(n,k) = Rt
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
   character*100 :: var

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
     !if(rank .eq. 0 .and. k .eq. 1) call dtime(t_start,t_end)
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
     !if(rank .eq. 0 .and. k .eq. En) call dtime(t_start,t_end)
     !if(rank .eq. 0 .and. k .eq. En) write(*,*) 'Lapack(Gaussian Gen) Time (time) = ',t_end,' sec' !Write to file
   else if(csp .eq. 2) then
     !if(rank .eq. 0 .and. k .eq. 1) call dtime(t_start,t_end)
     Np_t = Np
     ! Not working, not open source. Get seg fault in closed library
     ! Suppose to be faster then superLU
     !call Sparse_zgsmp_init(Np_t,G)
     !call Sparse_zgsmp(Np_t,ca,G)
     stop
     !if(rank .eq. 0 .and. k .eq. En) call dtime(t_start,t_end)
     !if(rank .eq. 0 .and. k .eq. En) write(*,*) 'WASP(Sparse) Time (time) = ',t_end,' sec' !Write to file
   else if(csp .eq. 3) then
     !if(rank .eq. 0 .and. k .eq. 1) call dtime(t_start,t_end)
     Np_t = Np
     !call dtime(t_start,t_end)
     if(cfr .eq. 0) call Sparse_superlu_init(Np_t,G) ! only initalize once
     call Sparse_superlu(Np_t,ca,G)
     !if(rank .eq. 0 .and. k .eq. En) call dtime(t_start,t_end)
     !if(rank .eq. 0 .and. k .eq. En) write(*,*) 'SuperLU(Sparse) Time (time) = ',t_end,' sec' !Write to file
   else if(csp .eq. 4) then
     !if(rank .eq. 0 .and. k .eq. 1) call dtime(t_start,t_end)
     !lapack tridiagonal gaussian solver
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
     if (INFO <  0) then   
       write (6,*) 'the i-th argument had an illegal value'
       stop 'ABORT'
     else if(INFO > 0) then
       write (6,*) 'U(i,i) is exactly zero; the matrix is singular and its inverse could not be computed'
       stop 'ABORT'
     end if
     !if(rank .eq. 0 .and. k .eq. En) call dtime(t_start,t_end)
     !if(rank .eq. 0 .and. k .eq. En) write(*,*) 'Lapack(Gaussian Tri) Time (time) = ',t_end,' sec' !Write to file
   else if(csp .eq. 5) then
     !if(rank .eq. 0 .and. k .eq. 1) call dtime(t_start,t_end)
     !lapack tridiagonal direct solver
     call Matxopt_ztrace_sub(G,1,Np,DL)
     call Matxopt_ztrace_sup(G,1,Np,DU)
     call Matxopt_ztrace(G,Np,D)
     call ZDTTRF(Np,DL,D,DU,INFO)
     if (INFO <  0) then   
       write (6,*) 'the i-th argument had an illegal value'
       stop 'ABORT'
     else if(INFO > 0) then
       write (6,*) 'U(i,i) is exactly zero; the matrix is singular'
       stop 'ABORT'
     end if
     !call Matxopt_zput_sub(DL,1,Np,G)
     !call Matxopt_zput_sup(DU,1,Np,G)
     !call Matxopt_zput(D,Np,G)

     G1 = eye
     call ZGTSV(Np,Np,DL,D,DU,G1,Np,INFO)
     G = G1
     if (INFO <  0) then   
       write (6,*) 'the i-th argument had an illegal value'
       stop 'ABORT'
     else if(INFO > 0) then
       write (6,*) 'U(i,i) is exactly zero; the matrix is singular and its inverse could not be computed'
       stop 'ABORT'
     end if
     !if(rank .eq. 0 .and. k .eq. En) call dtime(t_start,t_end)
     !if(rank .eq. 0 .and. k .eq. En) write(*,*) 'Lapack(Direct Tri) Time (time) = ',t_end,' sec' !Write to file
   end if

end subroutine Negf_Invers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Bound - Boundary conditions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Bound

   integer(kind=4) :: n
   real(kind=8) :: tm_ss, tm_se, wcut, Np1, Np2, dN

   tm_ss = abs(H(1,1))
   tm_se = abs(H(Np,Np))
   wcut = 2*sqrt(abs(tm_ss))

   !write(*,*) 'wcut, ecut = ',wcut,hbar*wcut
   !if(hbar*wcut .lt. E(k)) then
     !E(k) = 0
   !  write(*,*) 'Warning: Above cut off frequency',E(k),hbar*wcut
   !end if

   !! Calculating Fermi functions of electrons entering through Silicon Source and Drain !! 
   call Startup_DOS(vsi,E(k),kT1,No1)
   call Startup_DOS(vsi,E(k),kT2,No2)
   !Derived the fermi function for Bose-Einstein particles, not the same as electronic NEGF
   f1 = -No1*dlog(1-dexp((-E(k))/kT1)) ! 2D Source fermi function. Units of 1/m^2 
   f2 = -No2*dlog(1-dexp((-E(k))/kT2)) ! 2D drain fermi function. Units of 1/m^2
   f3 = No1*(1-dlog(1-dexp((-(E(k)))/kT1))) ! 2D Source fermi function. Units of 1/m^2
   f4 = No2*(1-dlog(1-dexp((-(E(k)))/kT2))) ! 2D drain fermi function. Units of 1/m^2

   f5 = f_e(k,3) !fermi function of scatters - don't know so we determine such that the current is conserved

   !f1 = 1/(dexp(E(k)/kT1)-1)
   !f2 = 1/(dexp(E(k)/kT2)-1)
   !write(*,*) 'f1, f2 = ',f1,f2
   !dN = E(k)/kT1*exp(E(k)/kT1)/(exp(E(k)/kT1)-1)**2*dT
   !write(*,*) 'df, dN = ',f2-f1,dN

   Sig1(:,:) = 0.0; Gam1(:,:) = 0.0; SigIn1(:,:) = 0.0; SigOut1(:,:) = 0.0
   ka1 =  i*(1-(1 - ((dcmplx((E(k)/hbar)**2) + zplus - dcmplx(Ds))/(2.0*dcmplx(tm_ss))))**2)**0.5
   ka1 = -i*cdlog((1 - ((dcmplx((E(k)/hbar)**2) + zplus - dcmplx(Ds))/(2.0*dcmplx(tm_ss)))) + ka1) ! Real Wave vector of broadened wavGp at drain (ACOS can be written in LOG Format See Matlab Manual) 		
   !ka1 = acos(1-real(((E(k)/hbar)**2 + zplus - Ds)/(2.0*dcmplx(tm_ss))))
   Sig1(1,1) = dcmplx(-tm_ss)*zexp(i*ka1) ! (omega^2) Self energy of source , CMPLX=Converts to Complex Data Type
   Gam1(1,1) = i*(Sig1(1,1) - conjg(Sig1(1,1))) ! (omega^2) Source broadening matrix  
   SigIn1(1,1) = Gam1(1,1)*dcmplx(f1) ! (omega^2/m^2) Inscattering term for source 
   SigOut1(1,1) = Gam1(1,1)*dcmplx(f3)

   Sig2(:,:) = 0.0; Gam2(:,:) = 0.0; SigIn2(:,:) = 0.0; SigOut2(:,:) = 0.0
   ka2 =  i*(1-(1 - ((dcmplx((E(k)/hbar)**2) + zplus - dcmplx(Ds))/(2.0*dcmplx(tm_se))))**2)**0.5 
   ka2 = -i*cdlog((1 - ((dcmplx((E(k)/hbar)**2) + zplus - dcmplx(Ds))/(2.0*dcmplx(tm_se)))) + ka2) ! Real Wave vector of broadened wavGp at drain (ACOS can be written in LOG Format See Matlab Manual)  
   !ka2 = acos(1-real(((E(k)/hbar)**2 + zplus - Ds)/(2.0*dcmplx(tm_se))))
   Sig2(Np,Np) = dcmplx(-tm_se)*zexp(i*ka2) ! (omega^2) Self energy of drain
   Gam2(Np,Np) = i*(Sig2(Np,Np) - conjg(Sig2(Np,Np))) ! (omega^2) Drain broadening matrix  
   SigIn2(Np,Np) = Gam2(Np,Np)*dcmplx(f2) ! (omega^2/m^2) Inscattering term for drain
   SigOut2(Np,Np) = Gam2(Np,Np)*dcmplx(f4) 

   Gamp(:) = 0.0
   !if(f5 .lt. 1e-50) f5 = 1e-12
   !if(maxval(SigOutp(:,k)) .gt. 0) call Matxopt_hilbert(0.5*SigOutp(:,k),(E/hbar)**2,dw,NE,SigInp(:,k))
   if(f5 .ne. 0.0 .and. tst .lt. 5) Gamp = -SigInp(:,k)/f5-0.5*i*(SigOutp(:,k)+SigInp(:,k))/f5 !omega^2

   G(:,:) = (0,0); 
   do n = 1, Np
     G1(n,n) = (0,0); G2(n,n) = (0,0)
   end do

   G1(1,1) = Sig1(1,1); G1(Np,Np) = Sig2(Np,Np) !(omega^2)
!if(rank .eq. 0 .and. c0 .eq. 2) write(*,*) 'G1(1,1),G1(Np,Np),Gn_ee(Np/2,c1)',G1(1,1),G1(Np,Np),Gn_ee(Np/2,c1)**2
   do n = 1, Np
     G2(n,n) = ((dcmplx((E(k)/hbar)**2) + zplus)*eye(n,n) + i*dcmplx(Gn_ee(n,c1)**2)) !(omega^2) Simplify Addition Because Only Diagonal Terms 
   end do

   call Matxopt_ztcpy(dcmplx(H),G,Np)
   call Matxopt_sub_zdtf(G2,G,Np) 
   call Matxopt_sub_zdtf2(G1,G,Np)
   call Matxopt_sub_zvtf2(Gamp,G,Np)
   !G = (G2 - dcmplx(H) - G1) ! Green's function (omega^2) 

end subroutine Negf_Bound

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! State_Plot - Plot Density of State
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine State_Plot

    integer(kind=4) :: r, n, j, err
    real(kind=4) :: l, mx
    character*100 :: var

    open (unit=19, file='A_p.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'A_p.dat'
    open (unit=20, file='Tr_p.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'Tr_p.dat'
    open (unit=21, file='I_p.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'I_p.dat'
    open (unit=22, file='Gf_p.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'Gf_p.dat'
    open (unit=23, file='Gp_p.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'Gp_p.dat'
    open (unit=24, file='Aall_p.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'Aall_p.dat'
    open (unit=25, file='Gfall_p.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'Gfall_p.dat'
    open (unit=26, file='Gpall_p.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'Gfall_p.dat'
    open (unit=27, file='Rh_p.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'Rh_p.dat'
    open (unit=28, file='Sc_p.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'Sc_p.dat'
    open (unit=29, file='Tau_p.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'Tau_p.dat'
    open (unit=30, file='Fermi_p.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'Fermi_p.dat'
    open (unit=31, file='Asum_p.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'Asum_p.dat'
    open (unit=32, file='Diff_p.dat', status='NEW', action='write', iostat=err)
    write(var,*) 'Asum_p.dat'
    open (unit=33, file='Tprt_p.dat', status='NEW', action='write', iostat=err)
    call Errors_Fileopen(err,var)
 
    100 FORMAT (5000(ES12.5,1X))

    101 format (3(ES12.5,1X))
    !102 format (4(ES12.5,1X))
    103 format (7(ES12.5,1X))
    104 format (1000(ES12.5,1X))

    ! Plot Transmission and Current
    mx = maxval(Tr1(:,:))

    do r=1, NE
       l = 0.0
       do n=1,Np
                               !A1(Np,Ne)
         write(19,101) E(r)/q,l,A1(n,r)/(2.0*pi) !LDOS Contour Plot
         write(22,101) E(r)/q,l,Gf1(n,r)/(2.0*pi) !LDOS Contour Plot
         write(23,101) E(r)/q,l,Gp1(n,r)/(2.0*pi) !LDOS Contour Plot
         l = Gdx(n) + l
       end do
       write(19,*);write(22,*);write(23,*)
    end do

    do r=1, NE
      write(20,100) E(r)/q,(Tr1(r,n), n=1,4),sum(Tr1(r,:))/(c1-1) !Transmission Plot
      write(21,103) E(r)/q,sum(Id1(r,:))*q,sum(Id2(r,:))*q,sum(Id3(r,:))*q,-sum(Id1(r,:))*q-sum(Id2(r,:))*q
    end do

    do r=1, NE
       l = 0.0
       do n=1,Np
         write(24,104) E(r)/q,l,(Aall(n,r,j)/(2.0*pi), j=1,c1) !LDOS Contour Plot
         write(25,104) E(r)/q,l,(Gnall(n,r,j)/(2.0*pi), j=1,c1) !LDOS Contour Plot
         write(26,104) E(r)/q,l,((Aall(n,r,j)-Gnall(n,r,j))/(2.0*pi), j=1,c1) !LDOS Contour Plot
         l = Gdx(n) + l
       end do
       write(24,*);write(25,*);write(26,*)
    end do

    do r=1, NE
      write(31,104) E(r)/q,sum(Aall(:,r,:)/(2.0*pi)),sum(Gnall(:,r,:)/(2.0*pi)),sum((Aall(:,r,:)-Gnall(:,n,:))/(2.0*pi)) !LDOS Contour Plot
    end do

    !First Derivative
    Tprtt(:,3) = 0.0
    do n=2, Np-1
      !Tprtt(n,3) = (Tprtt(n+1,1)-Tprtt(n,1))/(Gdx(n)) !forward difference
      Tprtt(n,4) = (Tprtt(n+1,1)-Tprtt(n-1,1))/(Gdx(n)+Gdx(n-1)) !central difference
    end do
    ! Plot Rho
    l = 0.0
    do n=1, Np
      write(27,104) l,((Rh1(n,j))/an*1E-6, j=1,4),sum(Rh1(n,:))/an*1E-6 !Electron plot
      if(n .gt. 1 .and. n .lt. Np) then
         write(33,*) l,Tprtt(n,1),Tprtt(n,2),Tprtt(n,3),abs(Itot/Tprtt(n,4))
      else if(n .eq. 1) then
         write(33,*) l,Tprtt(n,1),Tprtt(n,2),Tprtt(n,3),abs(Itot/Tprtt(2,4))
      else if(n .eq. Np) then
         write(33,*) l,Tprtt(n,1),Tprtt(n,2),Tprtt(n,3),abs(Itot/Tprtt(Np-1,4))
      end if 
      l = Gdx(n) + l
    end do 

    do n=1, Np
    G1(n,n) = sum(Rh1(n,:)/an)
    end do 
    call Matxopt_zdff(G1,dcmplx(dz),Np,G2)
    do n=1, Np
      write(32,104) l,real(G2(n,n)),Qtot/real(G2(n,n)),Lt/dT/real(G2(n,n)) !Electron plot
      l = Gdx(n) + l
    end do

    Cv = 0.0; Df=0.0
    do n=1, Np
      Df = Df + Qtot/real(G2(n,n))/Np
      Cv = Cv + Lt/dT/real(G2(n,n))/Np
    end do

    !LDOS Contour Plot
    do r=1, NE
       l = 0
       do n=1,Np
         write(28,103) E(r)/q,l,(Si1(n,r,j),j=1,4),sum(Si1(n,r,:)) !LDOS Contour Plot
         l = Gdx(n) + l
       end do
       write(28,*)
    end do

    do r=1, NE
     write(29,104) E(r)/q,(sum(Si1(:,r,j)),j=1,4),sum(Si1(:,r,:)) !LDOS Contour Plot
    end do

    do r=1, NE
     write(30,104) E(r)/q,f_e(r,1),f_e(r,2),f_e(r,3)
    end do

    close(unit=19); close(unit=20); close(unit=21); close(unit=22); close(unit=23)
    close(unit=24); close(unit=25); close(unit=26); close(unit=27); close(unit=28); 
    close(unit=29); close(unit=30); close(unit=31); close(unit=32); close(unit=33)

end subroutine State_Plot

subroutine Output_Tprt(Tprtt,Rh1,c)

    integer(kind=4) :: r, n, j, err
    real(kind=4) :: l, mx
    character*100 :: var
    character*80 :: Fname

    integer(kind=4), intent(in) :: c
    real(kind=8), dimension(:,:), intent(in) :: Rh1
    real(kind=8), dimension(:,:), intent(inout) :: Tprtt

    Fname = Utility_MakeFileName('Tprt_p',1,c,1,'dat')
    open (unit=33, file=Fname, status='NEW', action='write', iostat=err)
    Fname = Utility_MakeFileName('Tprt_NE',1,c,1,'dat')
    open (unit=34, file=Fname, status='NEW', action='write', iostat=err)
    call Errors_Fileopen(err,var)

    104 format (1000(ES12.5,1X))

    !First Derivative
    Tprtt(:,4) = 0.0
    do n=2, Np-1
      !Tprtt(n,3) = (Tprtt(n+1,1)-Tprtt(n,1))/(Gdx(n)) !forward difference
      Tprtt(n,4) = (Tprtt(n+1,1)-Tprtt(n-1,1))/(Gdx(n)+Gdx(n-1)) !central difference
    end do
    ! Plot Rho
    l = 0.0
    do n=1, Np
      if(n .gt. 1 .and. n .lt. Np) then
         write(33,*) l,Tprtt(n,1),Tprtt(n,2),Tprtt(n,3),abs(Itot/Tprtt(n,4))
      else if(n .eq. 1) then
         write(33,*) l,Tprtt(n,1),Tprtt(n,2),Tprtt(n,3),abs(Itot/Tprtt(2,4))
      else if(n .eq. Np) then
         write(33,*) l,Tprtt(n,1),Tprtt(n,2),Tprtt(n,3),abs(Itot/Tprtt(Np-1,4))
      end if 
      l = Gdx(n) + l
    end do

    do r = 1, NE
       l = 0
       do  n = 1, Np
         write(34,*) E(r),l,Tprt_NE(n,r)
         l = Gdx(n) + l
       end do
       write(34,*) 
    end do

    close(33);close(34)

end subroutine Output_Tprt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Loop_Plot - Appends to file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Loop_Plot

    character*100 :: op

    write(op,'(A,ES12.5,1X,ES12.5,1X,ES12.5)') 'w, J, Jtot = ',sqrt(Ds),Isb,Itot
    call debug_output(op)
    call flush(6)

end subroutine Loop_Plot

end module
