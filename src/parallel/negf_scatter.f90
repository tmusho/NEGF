!  Negf_Scatter.f90 
!
!  Authors: G. Walker - Vanderbilt University
!           T. Musho - Vanderbilt University         
!
!  Created: 1/25/11
!  Last Modified: 1/25/11
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
!  PROGRAM: Scatter_Routine(Module)
!
!  PURPOSE: Contains subroutine for Scatter Routine
!
!  SUBROUTINES: 
!               
!               
!                
!
!****************************************************************************
!
module Negf_Scatter

use Startup
use Matxopt
use Scatter

!use Negf

implicit none

  !include 'mpif.h'
  real(kind=8), allocatable, dimension(:,:) :: nse, nsa, pse, psa, SigInpNew, SigOutpNew
  real(kind=8), allocatable, dimension(:,:) :: phe, pha
  real(kind=8), private, dimension(4) :: norm

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Scatter - Calculate scattering matricies - Inelastic
!                This uses multiple phonon freq calc in the phonon code and the 
!                spatially varying phonon population. However, the
!                density of states are used in place of the BE dist.
!                Because we use the phonon DOS we account for transition that
!                are restricted due to the phonon band gaps. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Scatter_Couple(DH_p,Gnt,Gpt,Gn_pi,Gp_pi,Si,SigInp,SigOutp,chng,LDA,WORK)

    double precision dlange
    external dlange

    integer(kind=4) :: n,j, NE_off, NE_start
    real(kind=8) :: Esb_ph, Do_p, Wsb_ph

    real(kind=8), intent(out) :: chng
    real(kind=8), dimension(:), intent(in) :: DH_p
    real(kind=8), dimension(:,:), intent(in) :: Gnt,Gpt
    real(kind=8), dimension(:,:), intent(inout) :: SigInp,SigOutp
    real(kind=8), dimension(:,:,:), intent(inout) :: Gn_pi,Gp_pi,Si

    integer(kind=8), intent(in) :: LDA
    complex(kind=8), dimension(:), intent(inout) :: WORK

    SigInpNew(:,:) = 0.0; SigOutpNew(:,:) = 0.0;

    !Only take diagonal elements of scattering matrix
    !loop over all phonon eigenvalues
    do j=1,Nb_p

      Esb_ph = DH_p(j)
      Wsb_ph = Esb_ph/hbar_eV
      NE_off = ceiling(Esb_ph/dE)
      NE_start = ceiling(Ds/dE)
      !Nw = exp(-dE*NE_off/kT) ! Simplified Bose-Einstein distribution

      !! Reassigning density of states to accomodate scattering
      call Matxopt_eoshift(Gnt,NE,NE_off,nse) !Gn(E+hw) (1/eV-m^2)
      call Matxopt_eoshift(Gnt,NE,-NE_off,nsa) !Gn(E+hw) (1/eV-m^2)
      call Matxopt_eoshift(Gpt,NE,NE_off,pse) !Gn(E+hw) (1/eV-m^2)
      call Matxopt_eoshift(Gpt,NE,-NE_off,psa) !Gn(E+hw) (1/eV-m^2)
      phe = eoshift(Gn_pi(:,:,j), shift = -NE_start, dim = 2)
      pha = eoshift(Gn_pi(:,:,j)-Gp_pi(:,:,j), shift = -NE_start, dim = 2)

      do n=1,Np
        if(mtyp(n) .eq. 1) then
          op_cut = op_cut_si
        else
          op_cut = op_cut_ge
        end if
        if(Esb_ph .gt. op_cut) then
          call Defpot_Optical_ph(Esb_ph,Do_p,mtyp(n)) !|K|^2|A|^2 [eV^2]
        else
          call Defpot_Acoustic_ph(Esb_ph,Do_p,mtyp(n)) !|K|^2|A|^2 [eV^2]
        end if
        SigInpNew(n,:) = SigInpNew(n,:) + Do_p*(pha(n,j)*nsa(n,:) + phe(n,j)*nse(n,:))/Nb_p !units are now eV/m^2
        SigOutpNew(n,:) = SigOutpNew(n,:) + Do_p*(pha(n,j)*nsa(n,:)+ phe(n,j)*nse(n,:) + &
                            & pha(n,j)*pse(n,:) + phe(n,j)*psa(n,:))/Nb_p !eV/m^2     
      end do
      if(debug_l>5)Si(:,:,j) =  SigOutpNew(:,:)/hbar_eV !scattering rate 1/s-m^2
      !write(*,*) 'sum(phe),sum(pha) = ',sum(phe),sum(pha)
      !write(*,*) 'sum(Nw+1),sum(Nw) = ',Nw+1,Nw
      !write(*,*) 'sum(nse),sum(nsa) = ',sum(nse),sum(nsa)
      !write(*,*) 'sum(SigInpNew),sum(SigOutpNew) = ',sum(SigInpNew),sum(SigOutpNew)

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

end subroutine Negf_Scatter_Couple

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Scatter_Couple_Indv - Calculate scattering matricies
!                            This scattering uses phonon freq calc in phonon code
!                            and sums the phonon distribution spatial to use
!                            instead of the BE dist
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Scatter_Couple_Indv(DH,DH_p,Gnt,Gpt,Gn_pi,Gp_pi,Si,SigInp,SigOutp,chng,LDA,WORK)

    double precision dlange
    external dlange

    integer(kind=4) :: n,j, NE_off, fp, fe
    real(kind=8) :: Esb_ph, Do_p, Wsb_ph

    real(kind=8), intent(out) :: chng
    real(kind=8), dimension(:), intent(in) :: DH, DH_p
    real(kind=8), dimension(:,:), intent(in) :: Gnt, Gpt
    real(kind=8), dimension(:,:), intent(inout) :: SigInp,SigOutp
    real(kind=8), dimension(:,:,:), intent(inout) :: Gn_pi,Gp_pi,Si

    integer(kind=8), intent(in) :: LDA
    complex(kind=8), dimension(:), intent(inout) :: WORK
    SigInpNew(:,:) = 0.0; SigOutpNew(:,:) = 0.0; 
    fp = floor(sqrt(DH(1))/dE)
    fe = ceiling(sqrt(DH(Np))/dE)

    !Only take diagonal elements of scattering matrix
    !loop over all phonon eigenvalues
    do j=1,Nb_p

      Esb_ph = DH_p(j)
      Wsb_ph = Esb_ph/hbar_eV
      NE_off = ceiling(Esb_ph/dE)
      if(mwa .eq. 1) then !use Maxwellian Approx for Dist
        Nw = exp(-dE*NE_off/kT) ! Simplified Bose-Einstein distribution
      else
        Nw = 1/(exp((dE*NE_off)/kT)-1) !true BE dist
      end if

      !! Reassigning density of states to accomodate scattering
      !custom eoshift routines
      call Matxopt_eoshift(Gnt,NE,NE_off,nse) !Gn(E+hw) (1/eV-m^2)
      call Matxopt_eoshift(Gnt,NE,-NE_off,nsa) !Gn(E+hw) (1/eV-m^2)
      call Matxopt_eoshift(Gpt,NE,NE_off,pse) !Gn(E+hw) (1/eV-m^2)
      call Matxopt_eoshift(Gpt,NE,-NE_off,psa) !Gn(E+hw) (1/eV-m^2)

      !accumulate phonon density of states
      !make the phonon distribution a scalar value similar to what you would get from BE dist
      do n=1,Np
        !The phonon density of states has units of 1/omega^2 we need to convert to a number
        !without an energy association
        !Because the energy is interpolated use dE of the Energy spectrum
        phe(n,j) = sum(Gn_pi(n,:,j))*(dE_p/hbar_eV)**2/(2.0*pi) !Convert to a number FIXME: Noninteger values? N [1/1]
        pha(n,j) = sum(Gn_pi(n,:,j)-Gp_pi(n,:,j))*(dE_p/hbar_eV)**2/(2.0*pi) !Convert to a number N+1 [1/1]
      end do

      !calculate scatter matricies
      do n=1,Np
        if(mtyp(n) .eq. 1) then
          op_cut = op_cut_si
        else
          op_cut = op_cut_ge
        end if
        if(Esb_ph .gt. op_cut) then
          call Defpot_Optical_ph(Esb_ph,Do_p,mtyp(n)) !|K|^2|A|^2 [eV^2]
        else
          call Defpot_Acoustic_ph(Esb_ph,Do_p,mtyp(n)) !|K|^2|A|^2 [eV^2]
        end if
        SigInpNew(n,:) = SigInpNew(n,:) + Do_p*(pha(n,j)*nse(n,:) + phe(n,j)*nsa(n,:))/Nb_p !units are now eV/m^2
        SigOutpNew(n,:) = SigInpNew(n,:) + Do_p*(pha(n,j)*psa(n,:) + phe(n,j)*pse(n,:))/Nb_p !units are now eV/m^2
      end do
      if(debug_l>5)Si(:,:,j) =  (SigInpNew(:,:)+SigOutpNew(:,:))/hbar_eV !scattering rate 1/s-m^2
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

end subroutine Negf_Scatter_Couple_Indv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Scatter_Couple_Indv_Parallel - Calculate scattering matricies - Inelastic
!                            This scattering uses phonon freq calc in phonon code
!                            and sums the phonon distribution spatial to use
!                            instead of the BE dist
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Scatter_Couple_Indv_Parallel(DH,DH_p,Gnt,Gpt,Gn_pi,Gp_pi,Si,Sit,SigInp,SigOutp,chng,LDA,WORK, &
                                             disp_scat,rcount_scat,rcount_np,disp_np)

    double precision dlange
    external dlange

    integer(kind=4) :: n,j, NE_off, ldc, lds, ierr
    real(kind=8) :: Esb_ph, Do_p, Wsb_ph

    real(kind=8), intent(out) :: chng
    real(kind=8), dimension(:), intent(in) :: DH,DH_p
    real(kind=8), dimension(:,:), intent(in) :: Gnt,Gpt
    real(kind=8), dimension(:,:), intent(inout) :: SigInp,SigOutp
    real(kind=8), dimension(:,:,:), intent(inout) :: Gn_pi,Gp_pi,Si,Sit

    integer(kind=8), intent(in) :: LDA
    complex(kind=8), dimension(:), intent(inout) :: WORK
    integer(kind=4), dimension(:), intent(in) :: disp_scat,rcount_scat,rcount_np,disp_np

    SigInpNew(:,:) = 0.0; SigOutpNew(:,:) = 0.0;

    ldc = disp_np(rank+1)
    lds = ldc + rcount_np(rank+1)

    !Need to broadcast matricies to all processors
    call MPI_Bcast(Gnt,Np*NE,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(Gpt,Np*NE,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    !Only take diagonal elements of scattering matrix
    !loop over all phonon eigenvalues
    if(ldc .ge. 0) then
    do j=ldc+1, lds
      Esb_ph = DH_p(j) !eV phonon eigenvalues
      Wsb_ph = Esb_ph/hbar_eV !omega phonon frequency
      NE_off = ceiling(Esb_ph/dE)
      if(mwa .eq. 1) then !use Maxwellian Approx for Dist
        Nw = exp(-dE*NE_off/kT) ! Simplified Bose-Einstein distribution
      else
        Nw = 1/(exp((dE*NE_off)/kT)-1) !true BE dist
      end if

      !! Reassigning density of states to accomodate scattering
      !custom eoshift routines
      call Matxopt_eoshift(Gnt,NE,NE_off,nse) !Gn(E+hw) (1/eV-m^2)
      call Matxopt_eoshift(Gnt,NE,-NE_off,nsa) !Gn(E+hw) (1/eV-m^2)
      call Matxopt_eoshift(Gpt,NE,NE_off,pse) !Gn(E+hw) (1/eV-m^2)
      call Matxopt_eoshift(Gpt,NE,-NE_off,psa) !Gn(E+hw) (1/eV-m^2)

      !accumulate phonon density of states
      !make the phonon distribution a scalar value similar to what you would get from BE dist
      !done up front

      !calculate scatter matricies
      do n=1,Np !this iterates spatially across device
        if(mtyp(n) .eq. 1) then
          op_cut = op_cut_si
        else
          op_cut = op_cut_ge
        end if
        if(Esb_ph .gt. op_cut) then
          call Defpot_Optical_ph(Esb_ph,Do_p,mtyp(n)) !|K|^2|A|^2 [eV^2]
        else
          call Defpot_Acoustic_ph(Esb_ph,Do_p,mtyp(n)) !|K|^2|A|^2 [eV^2]
        end if
        SigInpNew(n,:) = SigInpNew(n,:) + Do_p*(pha(n,j)*nse(n,:) + phe(n,j)*nsa(n,:))/Nb_p !units are now eV/m^2
        SigOutpNew(n,:) = SigOutpNew(n,:) + Do_p*(pha(n,j)*psa(n,:) + phe(n,j)*pse(n,:))/Nb_p !units are now eV/m^2

      end do
      !for outputting scattering rate per frequency
      if(debug_l>5)Si(:,:,j) =  (SigInpNew(:,:)+SigOutpNew(:,:))/hbar_eV !scattering rate 1/s-m^2
    end do
      call MPI_Barrier(MPI_COMM_WORLD,ierr)  !Sync
    else
      call MPI_Barrier(MPI_COMM_WORLD,ierr)  !Sync
    end if

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

end subroutine Negf_Scatter_Couple_Indv_Parallel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Scatter_Couple_Indv_Elast - Calculate scattering matricies - Elastic
!                            This scattering uses phonon freq calc in phonon code
!                            and sums the phonon distribution spatial to use
!                            instead of the BE dist.
!                            Similar to above except the electron interaction is
!                            elastic in nature.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Scatter_Couple_Indv_Elast(DH_p,Gnt,Gpt,Gn_pi,Gp_pi,Si,SigInp,SigOutp,chng,LDA,WORK)

    double precision dlange
    external dlange

    integer(kind=4) :: n,j, NE_off
    real(kind=8) :: Esb_ph, Do_p, Wsb_ph

    real(kind=8), intent(out) :: chng
    real(kind=8), dimension(:), intent(in) :: DH_p
    real(kind=8), dimension(:,:), intent(in) :: Gnt,Gpt
    real(kind=8), dimension(:,:), intent(inout) :: SigInp,SigOutp
    real(kind=8), dimension(:,:,:), intent(inout) :: Gn_pi,Gp_pi,Si

    integer(kind=8), intent(in) :: LDA
    complex(kind=8), dimension(:), intent(inout) :: WORK
    SigInpNew(:,:) = 0.0; SigOutpNew(:,:) = 0.0;

    !Only take diagonal elements of scattering matrix
    !loop over all phonon eigenvalues
    do j=1,Nb_p

      Esb_ph = DH_p(j)
      Wsb_ph = Esb_ph/hbar_eV
      NE_off = 0
      if(mwa .eq. 1) then !use Maxwellian Approx for Dist
        Nw = exp(-dE*NE_off/kT) ! Simplified Bose-Einstein distribution
      else
        Nw = 1/(exp((dE*NE_off)/kT)-1) !true BE dist
      end if

      !! Reassigning density of states to accomodate scattering
      call Matxopt_eoshift(Gnt,NE,NE_off,nse) !Gn(E+hw) (1/eV-m^2)
      call Matxopt_eoshift(Gnt,NE,-NE_off,nsa) !Gn(E+hw) (1/eV-m^2)
      call Matxopt_eoshift(Gpt,NE,NE_off,pse) !Gn(E+hw) (1/eV-m^2)
      call Matxopt_eoshift(Gpt,NE,-NE_off,psa) !Gn(E+hw) (1/eV-m^2)

      !accumulate phonon density of states
      !make the phonon distribution a scalar value similar to what you would get from BE dist
      do n=1,Np
        !The phonon density of states has units of 1/omega^2 we need to convert to a number
        !without an energy association
        !Because the energy is interpolated use dE of the Energy spectrum
        phe(n,1) = sum(Gn_pi(n,:,j))*(dE_p/hbar_eV)**2/(2.0*pi) !Convert to a number FIXME: Noninteger values? [1/1]
        pha(n,1) = sum(Gn_pi(n,:,j)-Gp_pi(n,:,j))*(dE_p/hbar_eV)**2/(2.0*pi) !Convert to a number N+1 [1/1]
      end do

      !write(*,*) 'phe(n,1) = ',sum(phe(:,1)),f_p(j)
      !write(*,*) 'pha(n,1) = ',sum(pha(:,1)),f_p(j)

      !calculate scatter matricies
      do n=1,Np
        if(mtyp(n) .eq. 1) then
          op_cut = op_cut_si
        else
          op_cut = op_cut_ge
        end if
        if(Esb_ph .gt. op_cut) then
          call Defpot_Optical_ph(Esb_ph,Do_p,mtyp(n)) !|K|^2|A|^2 [eV^2]
        else
          call Defpot_Acoustic_ph(Esb_ph,Do_p,mtyp(n)) !|K|^2|A|^2 [eV^2]
        end if
        SigInpNew(n,:) = SigInpNew(n,:) + Do_p*(pha(n,1)*nse(n,:) + phe(n,:)*nsa(n,:))/Nb_p !units are now eV/m^2
        SigOutpNew(n,:) = SigOutpNew(n,:) + Do_p*(pha(n,1)*psa(n,:) + phe(n,1)*pse(n,:))/Nb_p !units are now eV/m^2
      end do
      if(debug_l>5)Si(:,:,j) =  (SigInpNew(:,:)+SigOutpNew(:,:))/hbar_eV !scattering rate 1/s-m^2
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

end subroutine Negf_Scatter_Couple_Indv_Elast

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Defpot_Optical_ph - This is |K|^2|A|^2, derived from the coupling of classic
!                     lattice movement and quantum. See Lundstrom txt pg78
!                     For eletron-optical phonon interaction
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Defpot_Optical_ph(Eph,D,typ)

    real(kind=8) :: D_o, rho, m_a, m_e, m_b, omega, omega_cut, beta, vs
    integer(kind=4), intent(in) :: typ
    real(kind=8), intent(in) :: Eph
    real(kind=8), intent(out) :: D

    !Ref: C. Jacoboni and L. Reggiani, 1983.
    !Ref: C. Jacoboni and P. Lugli, 1989. 
    if(typ .eq. 1) then !Si
      !D_o = 2.2e10 !(eV/m)
      D_o = 0.7e10 !(eV/m)
      rho = 2.338e3 !(kg/m^3)
      m_a = 1.6913e-26 ! kg mass of atom
      m_b = ms_b ! bulk effective mass
      vs = 9033 ! sound speed
      m_e = ms ! effective mass
    else !Ge
      !D_o = 5.5e10 !(eV/m) deformation term
      rho = 5.32e3 !(kg/m^3) density
      D_o = 0.79e10 !(eV/m)
      m_a = 4.3744e-26 ! kg mass of atom
      m_b =  mg_b ! bulk effective mass
      vs = 5310 !sound speed
      m_e = mg ! effective mass
    end if
    omega = Eph/hbar_eV
    omega_cut = op_cut/hbar_eV
    beta = omega/vs
!write(*,*) 'omega,beta =',omega,beta
    !D = pi*m_e*D_o**2*q**2/(hbar*omega*beta*m_a*vs*m_a)*hbar_eV !(eV)
    !D = 0.5*pi*D_o**2*q**2/(hbar*vs*beta*m_a*omega)*hbar_eV !(eV)
    ! D = pi*D_o**2/(2*rho*omega)
    D = D_o**2*beta**2*hbar_eV*q/(32*rho*omega*an)
    !D = D_o**2*beta**2*2*hbar/(rho*omega*an)
!write(*,*) 'Dodp =',D
end subroutine Defpot_Optical_ph

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Defpot_Acoustic_ph - This is |K|^2|A|^2, derived from the coupling of classic
!                      lattice movement and quantum. See Lundstrom txt pg78
!                      For eletron-acoustic phonon interaction
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Defpot_Acoustic_ph(Eph,D,typ)

    real(kind=8) :: D_o, rho, m_a, m_e, m_b, omega, beta, vs
    integer(kind=4), intent(in) :: typ
    real(kind=8), intent(in) :: Eph
    real(kind=8), intent(out) :: D

    !Ref: C. Jacoboni and L. Reggiani, 1983.
    !Ref: C. Jacoboni and P. Lugli, 1989. 
    if(typ .eq. 1) then !Si
      !D_o = 7.2 !(eV) Si
      D_o = 9.0 
      rho = 2.338e3 !(kg/m^3)
      m_a = 1.6913e-26 ! kg mass of atom
      m_b = ms_b ! bulk effective mass
      vs = 9033
      m_e = ms ! effective mass
    else !Ge
      !D_o = 9.58 !(eV) Ge 
      D_o = 9.0
      rho = 5.32e3 !(kg/m^3)
      m_a = 4.3744e-26 ! kg mass of atom
      m_b =  mg_b ! bulk effective mass
      vs = 5310
      m_e = mg ! effective mass
    end if
    omega = Eph/hbar_eV
    beta = omega/vs
    !D = pi*m_e*D_o**2*q**2/(hbar*m_a*vs*beta*hbar*beta)*hbar_eV !(eV)
    !D = 0.5*pi*beta*D_o**2*q**2/(hbar*vs*m_a*omega)*hbar_eV !(1/eV)
    D = D_o**2*beta**3*hbar_eV*q/(24*rho*vs*an)
!write(*,*) 'Dadp =',D
end subroutine Defpot_Acoustic_ph

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Scatter_Multiple - Calculate scattering matricies - Inelastic
!                Use the phonon freq from phonon code but use simple BE Dist
!                to determine population. Phonons uniform throughout structure
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Scatter_Multiple(DH_p,Gnt,Gpt,Si,SigInp,SigOutp,chng,LDA,WORK)

    double precision dlange
    external dlange

    integer(kind=4) :: n, j, NE_off
    real(kind=8) :: Esb_ph, Do_p

    real(kind=8), intent(out) :: chng
    real(kind=8), dimension(:), intent(in) :: DH_p
    real(kind=8), dimension(:,:), intent(in) :: Gnt,Gpt
    real(kind=8), dimension(:,:), intent(inout) :: SigInp,SigOutp
    real(kind=8), dimension(:,:,:), intent(inout) :: Si

    integer(kind=8), intent(in) :: LDA
    complex(kind=8), dimension(:), intent(inout) :: WORK
    SigInpNew = 0; SigOutpNew = 0

    do j=1,Nb_p

      Esb_ph = DH_p(j)
      NE_off = ceiling(Esb_ph/dE)
      if(mwa .eq. 1) then !use Maxwellian Approx for Dist
        Nw = exp(-dE*NE_off/kT) ! Simplified Bose-Einstein distribution
      else
        Nw = 1/(exp((dE*NE_off)/kT)-1) !true BE dist
      end if

      !! Reassigning density of states to accomodate scattering
      call Matxopt_eoshift(Gnt,NE,NE_off,nse) !Gn(E+hw)  (1/eV-m^2)
      call Matxopt_eoshift(Gnt,NE,-NE_off,nsa) !Gn(E+hw) (1/eV-m^2)
      call Matxopt_eoshift(Gpt,NE,NE_off,pse) !Gn(E+hw) (1/eV-m^2)
      call Matxopt_eoshift(Gpt,NE,-NE_off,psa) !Gn(E+hw) (1/eV-m^2)
        
     ! Nw = phonon fermi function, So = deformation potent
     ! Dem=(1+Nw)*So     Dab=Nw*So
     ! There is a mistake in Datta's book when he defines SigIn
     ! Gam expression is correct
      do n=1,Np
        if(mtyp(n) .eq. 1) then
          op_cut = op_cut_si
        else
          op_cut = op_cut_ge
        end if
        if(Esb_ph .gt. op_cut) then
          call Defpot_Optical_ph(Esb_ph,Do_p,mtyp(n)) !|K|^2|A|^2 [eV^2]
        else
          call Defpot_Acoustic_ph(Esb_ph,Do_p,mtyp(n)) !|K|^2|A|^2 [eV^2]
        end if
      SigInpNew(n,:) = SigInpNew(n,:) + Do_p*((1+Nw)*nse(n,:) + Nw*nsa(n,:))/Nb_p !units are now eV/m^2
      !SigOutpNew(n,:) = SigOutpNew(n,:) + Do_p*((1+Nw)*nse(n,:)+ Nw*nsa(n,:) + &
      !                      & (1+Nw)*psa(n,:) + Nw*pse(n,:))*Esb_ph !eV/m^2
      SigOutpNew(n,:) = SigOutpNew(n,:) + Do_p*((1+Nw)*psa(n,:) + Nw*pse(n,:))/Nb_p !eV/m^2
      end do
      if(debug_l>5)Si(:,:,j) =  (SigInpNew(:,:)-SigOutpNew(:,:))/hbar_eV !scattering rate 1/s-m^2    
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

end subroutine Negf_Scatter_Multiple

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Scatter_Multiple_Elast - Calculate scattering matricies - Elastic Scattering
!                Use the phonon freq from phonon code but use simple BE Dist
!                to determine population. Phonons uniform throughout structure
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Scatter_Multiple_Elast(DH_p,Gnt,Gpt,Si,SigInp,SigOutp,chng,LDA,WORK)

    double precision dlange
    external dlange

    integer(kind=4) :: n, j, NE_off
    real(kind=8) :: Esb_ph, Do_p

    real(kind=8), intent(out) :: chng
    real(kind=8), dimension(:), intent(in) :: DH_p
    real(kind=8), dimension(:,:), intent(in) :: Gnt,Gpt
    real(kind=8), dimension(:,:), intent(inout) :: SigInp,SigOutp
    real(kind=8), dimension(:,:,:), intent(inout) :: Si

    integer(kind=8), intent(in) :: LDA
    complex(kind=8), dimension(:), intent(inout) :: WORK
    SigInpNew = 0; SigOutpNew = 0

    do j=1,Nb_p

      Esb_ph = DH_p(j)
      NE_off = 0
      if(mwa .eq. 1) then !use Maxwellian Approx for Dist
        Nw = exp(-dE*NE_off/kT) ! Simplified Bose-Einstein distribution
      else
        Nw = 1/(exp((dE*NE_off)/kT)-1) !true BE dist
      end if

      !! Reassigning density of states to accomodate scattering
      call Matxopt_eoshift(Gnt,NE,NE_off,nse) !Gn(E+hw) (1/eV-m^2)
      call Matxopt_eoshift(Gnt,NE,-NE_off,nsa) !Gn(E+hw) (1/eV-m^2)
      call Matxopt_eoshift(Gpt,NE,NE_off,pse) !Gn(E+hw) (1/eV-m^2)
      call Matxopt_eoshift(Gpt,NE,-NE_off,psa) !Gn(E+hw) (1/eV-m^2)
        
     ! Nw = phonon fermi function, So = deformation potent
     ! Dem=(1+Nw)*So     Dab=Nw*So
     ! There is a mistake in Datta's book when he defines SigIn
     ! Gam expression is correct
      do n=1,Np
        if(mtyp(n) .eq. 1) then
          op_cut = op_cut_si
        else
          op_cut = op_cut_ge
        end if
        if(Esb_ph .gt. op_cut) then
          call Defpot_Optical_ph(Esb_ph,Do_p,mtyp(n)) !|K|^2|A|^2 [eV^2]
        else
          call Defpot_Acoustic_ph(Esb_ph,Do_p,mtyp(n)) !|K|^2|A|^2 [eV^2]
        end if
      SigInpNew(n,:) = SigInpNew(n,:) + Do_p*((1+Nw)*nse(n,:) + Nw*nsa(n,:))/Nb_p !units are now eV/m^2
      SigOutpNew(n,:) = SigOutpNew(n,:) + Do_p*((1+Nw)*psa(n,:) + Nw*pse(n,:))/Nb_p !eV/m^2
      end do
      if(debug_l>5)Si(:,:,j) =  (SigInpNew(:,:)+SigOutpNew(:,:))/hbar_eV !scattering rate 1/s-m^2    
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

end subroutine Negf_Scatter_Multiple_Elast

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Scatter_Single - Calculate scattering matricies - Inelastic Scattering
!                Use a user specified single phonon freq and BE dist to scatter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Scatter_Single(Gnt,Gpt,Si,SigInp,SigOutp,chng,LDA,WORK)

    double precision dlange
    external dlange

    integer(kind=4) :: n, NE_off, NE_start, sh_flag, j
    real(kind=8) :: Esb_ph, Do_p

    real(kind=8), intent(out) :: chng
    real(kind=8), dimension(:,:), intent(in) :: Gnt,Gpt
    real(kind=8), dimension(:,:), intent(inout) :: SigInp,SigOutp
    real(kind=8), dimension(:,:,:), intent(out) :: Si

    integer(kind=8), intent(in) :: LDA
    complex(kind=8), dimension(:), intent(inout) :: WORK
    SigInpNew(:,:) = 0.0; SigOutpNew(:,:) = 0

    sh_flag = 0

    do n=1,Np
      if(mtyp(n) .eq. 1) then
        Esb_ph = Eph_si
        op_cut = op_cut_si
      else
        Esb_ph = Eph_ge
        op_cut = op_cut_ge
      end if
      NE_off = ceiling(Esb_ph/dE)
      if(mwa .eq. 1) then !use Maxwellian Approx for Dist
        Nw = exp(-dE*NE_off/kT) ! Simplified Bose-Einstein distribution
      else
        Nw = 1/(exp((dE*NE_off)/kT)-1) !true BE dist
      end if

      !! Reassigning density of states to accomodate scattering
      !! Need to take the view point of the final electron state
      if(sh_flag .ne.  mtyp(n)) then !negative = shift up, positive = shift down
        call Matxopt_eoshift(Gnt,NE,NE_off,nse)  !Gn(E+hw) (1/eV-m^2)
        call Matxopt_eoshift(Gnt,NE,-NE_off,nsa) !Gn(E+hw) (1/eV-m^2)
        call Matxopt_eoshift(Gpt,NE,NE_off,pse)  !Gp(E+hw) (1/eV-m^2)
        call Matxopt_eoshift(Gpt,NE,-NE_off,psa) !Gp(E+hw) (1/eV-m^2)
        sh_flag = mtyp(n)
      end if

      ! Nw = phonon fermi function, So = deformation potent
      ! Dem=(1+Nw)*So     Dab=Nw*So
      ! There is a mistake in Datta's book when he defines SigIn
      ! Gam expression is correct
      if(Esb_ph .gt. op_cut) then
        call Defpot_Optical_ph(Esb_ph,Do_p,mtyp(n)) !|K|^2|A|^2 [eV^2]
      else
        call Defpot_Acoustic_ph(Esb_ph,Do_p,mtyp(n)) !|K|^2|A|^2 [eV^2]
      end if
      SigInpNew(n,:) = SigInpNew(n,:) + Do_p*((1+Nw)*nse(n,:) + Nw*nsa(n,:)) !units are now eV/m^2
      SigOutpNew(n,:) = SigOutpNew(n,:) + Do_p*((1+Nw)*psa(n,:) + Nw*pse(n,:)) !eV/m^2
    end do
    if(debug_l>5)Si(:,:,1) =  (SigInpNew(:,:)+SigOutpNew(:,:))/hbar_eV !scattering rate 1/s-m^2    
    
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
          chng = 1.0;
        else
          chng = 0.0;
        end if 
      end if

    if(debug_l>2) write(6,'(/,A,ES11.4)') 'Scatter chg = ',chng ! Write change of charge to screen

end subroutine Negf_Scatter_Single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Negf_Scatter_Single_Elast - Calculate scattering matricies - Elastic Scatter
!                Use a user specified single phonon freq and BE dist to scatter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Negf_Scatter_Single_Elast(Gnt,Gpt,Si,SigInp,SigOutp,chng,LDA,WORK)

    double precision dlange
    external dlange

    integer(kind=4) :: n, NE_off, NE_start, sh_flag
    real(kind=8) :: Esb_ph, Do_p

    real(kind=8), intent(out) :: chng
    real(kind=8), dimension(:,:), intent(in) :: Gnt,Gpt
    real(kind=8), dimension(:,:), intent(inout) :: SigInp,SigOutp
    real(kind=8), dimension(:,:,:), intent(inout) :: Si

    integer(kind=8), intent(in) :: LDA
    complex(kind=8), dimension(:), intent(inout) :: WORK
    SigInpNew(:,:) = 0.0; SigOutpNew(:,:) = 0

    sh_flag = 0

    do n=1,Np
      if(mtyp(n) .eq. 1) then
        Esb_ph = Eph_si
        op_cut = op_cut_si
      else
        Esb_ph = Eph_ge
        op_cut = op_cut_ge
      end if
      NE_off = 0 !elastic
      if(mwa .eq. 1) then !use Maxwellian Approx for Dist
        Nw = exp(-dE*NE_off/kT) ! Simplified Bose-Einstein distribution
      else
        Nw = 1/(exp((dE*NE_off)/kT)-1) !true BE dist
      end if

      if(sh_flag .ne.  1) then
        call Matxopt_eoshift(Gnt,NE,NE_off,nse) !Gn(E+hw)
        call Matxopt_eoshift(Gnt,NE,-NE_off,nsa) !Gn(E+hw)
        call Matxopt_eoshift(Gpt,NE,NE_off,pse) !Gn(E+hw)
        call Matxopt_eoshift(Gpt,NE,-NE_off,psa) !Gn(E+hw)
        sh_flag = 1
      end if

      ! Nw = phonon fermi function, So = deformation potent
      ! Dem=(1+Nw)*So     Dab=Nw*So
      ! There is a mistake in Datta's book when he defines SigIn
      ! Gam expression is correct
      if(Esb_ph .gt. op_cut) then
        call Defpot_Optical_ph(Esb_ph,Do_p,mtyp(n)) !|K|^2|A|^2 [eV^2]
      else
        call Defpot_Acoustic_ph(Esb_ph,Do_p,mtyp(n)) !|K|^2|A|^2 [eV^2]
      end if
      SigInpNew(n,:) = SigInpNew(n,:) + Do_p*((1+Nw)*nse(n,:) + Nw*nsa(n,:))*Esb_ph !units are now eV/m^2
      SigOutpNew(n,:) = SigOutpNew(n,:) + Do_p*((1+Nw)*psa(n,:) + Nw*pse(n,:))*Esb_ph !eV/m^2
      end do
      if(debug_l>5)Si(:,:,1) =  (SigInpNew(:,:)+SigOutpNew(:,:))/hbar_eV !scattering rate 1/s-m^2    

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

end subroutine Negf_Scatter_Single_Elast
end module
