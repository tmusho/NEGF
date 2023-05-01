!  init.f90 
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
!  PROGRAM: Init(Module)
!
!  PURPOSE: Contains subroutines de/allocating variables 
!
!  SUBROUTINES: Init_Main - Main Startup, subroutine calls
!               Init_Allocate - Allocates arrays for NEGF
!               Init_Deallocate - Deallocates arrays 
!
!****************************************************************************
!
module Init

!Use global variables from module
use Negf

implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Init_Main - Main Startup, subroutine calls
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Init_Main

   call Init_Allocate
   call Init_Var

   !! LAPACK Variable !!
   LDA = Np ! The leading dimension of the array LDA >= max(1,M)

end subroutine Init_Main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Init_Allocate - Allocates arrays for NEGF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Init_Allocate

   integer(kind=4) :: err
   character*120 :: var

   !! Allocation Variables !!
   if(debug_l>5 .and. rank .eq. 0) then !Memory conservation - not plotting don't allocate
     allocate(A1(Np,NE), Gn1(Np,NE), stat=err)
     allocate(Aall(Np,NE,Nb), Ua1(Np,NV+2), stat=err)
     allocate(Id1(NE,Nb), Id2(NE,Nb), Id3(NE,Nb), Tr1(NE,Nb), Tr2(NE,Nb), U1(Np,Nb), stat=err)
   end if

   allocate(Si(Np,NE,Nb_p), stat=err)

   if(rank .eq. 0) then
     allocate(Sit(Np,NE,Nb_p), Si1(Np,NE,Nb_p), stat=err)
   end if

   if(psc .eq. 1) then 
     allocate(nse(Np,NE), nsa(Np,NE), pse(Np,NE), psa(Np,NE), SigInpNew(Np,NE), SigOutpNew(Np,NE), stat=err)
     allocate(phe(Np,NE), pha(Np,NE), stat=err)
   else if(rank .eq. 0) then !only rank0 needs scattering matrices when not parallel scattering
     allocate(nse(Np,NE), nsa(Np,NE), pse(Np,NE), psa(Np,NE), SigInpNew(Np,NE), SigOutpNew(Np,NE), stat=err)
     allocate(phe(Np,NE), pha(Np,NE), stat=err)
   end if

   allocate(Sig1(Np,Np), Gam1(Np,Np), SigIn1(Np,Np), Sig2(Np,Np), Gam2(Np,Np), stat=err) ! Allocate Matrix Variables 
   allocate(SigIn2(Np,Np), G(Np,Np), Gtc(Np,Np), Gn(Np,NE), Gp(Np,NE), Rh(Np), T(Np,Np), stat=err)
   allocate(G1(Np,Np), G2(Np,Np), G3(Np,Np), G4(Np,Np), G5(Np,Np), Tr(NE,2), stat=err)
   allocate(NN(Np), Ud(Np), No(Np,Nb), IT(NE), I1(NE), I2(NE), I3(NE), stat=err) 
   allocate(A(Np,NE), Gf(Np,NE), f_e(NE), Rh1(Np), Ua(Np), stat=err)
   allocate(SigOut1(Np,Np), SigOut2(Np,Np), SigInp(Np,NE), SigOutp(Np,NE), Gamp(Np), stat=err)
   allocate(IV(NV+1,2),MV(NV+1,2), stat=err)

   write(var,*) 'NEGF Matricies Allocate'
   call Errors_Allocate(err,var)

end subroutine Init_Allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Init_Var - Allocates arrays for NEGF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Init_Var

  if(debug_l>5 .and. rank .eq. 0) then !Memory conservation - not plotting don't allocate
    Aall(:,:,:) = 0.0; A1(:,:) = 0.0; Gn1(:,:)=0.0
    Tr1(:,:) = 0.0; Tr2(:,:) = 0.0; Ua1(:,:) =0.0
    U1(:,:) = 0.0; Id1(:,:)= 0.0; Id2(:,:)= 0.0; Id3(:,:)= 0.0
  end if

  Si(:,:,:)=0.0;
   if(rank .eq. 0) then
     Sit(:,:,:) = (0.0,0.0); Si1(:,:,:) = 0.0; 
   end if

  if(psc .eq. 1) then 
    nse(:,:) = (0.0,0.0); nsa(:,:) = (0.0,0.0); pse(:,:) = (0.0,0.0);
    psa(:,:) = (0.0,0.0); SigInpNew(:,:) = (0.0,0.0); SigOutpNew(:,:) = (0.0,0.0)
    phe(:,:) = (0.0,0.0); pha(:,:) = (0.0,0.0)
  else if(rank .eq. 0) then
    nse(:,:) = (0.0,0.0); nsa(:,:) = (0.0,0.0); pse(:,:) = (0.0,0.0);
    psa(:,:) = (0.0,0.0); SigInpNew(:,:) = (0.0,0.0); SigOutpNew(:,:) = (0.0,0.0)
    phe(:,:) = (0.0,0.0); pha(:,:) = (0.0,0.0)
  end if

  No(:,:) = 0.0 ! Set Array to Zero
  NN(:) = 0.0; No(:,:) = 0.0; Ud(:) = 0.0
  G(:,:) = (0.0,0.0); Gtc(:,:) = (0.0,0.0)
  G1(:,:) = (0.0,0.0); G2(:,:) = (0.0,0.0); G4(:,:) = (0.0,0.0); G5(:,:) = (0.0,0.0)
  Gn(:,:) = (0.0,0.0); Rh(:) = (0.0,0.0); Gp(:,:) = (0.0,0.0)
  Sig1(:,:) = (0.0,0.0); Sig2(:,:) = (0.0,0.0)
  SigIn1(:,:) = (0.0,0.0); SigIn2(:,:) = (0.0,0.0)
  Gam1(:,:) = (0.0,0.0); Gam2(:,:) = (0.0,0.0)
  Rh1(:) = 0.0; Tr(:,:) = 0.0; A(:,:) = 0.0;Gamp(:) = (0.0,0.0);
  SigOut1(:,:) = (0.0,0.0); SigOut2(:,:) = (0.0,0.0); SigInp(:,:) = (0.0,0.0); SigOutp(:,:) = (0.0,0.0)
  I1(:)=0.0; I2(:)=0.0; I3(:)=0.0; Ua(:)=0.0
  IV(:,:) = 0.0; MV(:,:) = 0.0; Gf(:,:) = 0.0


end subroutine Init_Var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Init_Deallocate - Deallocates arrays
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Init_Deallocate

   integer(kind=4) :: err
   character*120 :: var

   !! Deallocate Dynamic Memory !!
   if(debug_l>5 .and. rank .eq. 0) then !Memory conservation - not plotting don't allocate
     deallocate(A1, Gn1, stat=err)
     deallocate(Aall, stat=err)
     deallocate(Id1, Id2, Id3, Tr1, Tr2, U1, Ua1, stat=err)
   end if

   deallocate(Si, stat=err)
   if(rank .eq. 0) then
     deallocate(Sit, Si1, stat=err)
   end if

   if(psc .eq. 1) then 
    deallocate(nse, nsa, pse, psa, SigInpNew, SigOutpNew, stat=err)
    deallocate(phe, pha, stat=err)
   else if(rank .eq. 0) then
    deallocate(nse, nsa, pse, psa, SigInpNew, SigOutpNew, stat=err)
    deallocate(phe, pha, stat=err)
   end if

   deallocate(Sig1, Gam1, SigIn1, Sig2, Gam2, stat=err) ! Allocate Matrix Variables 
   deallocate(SigIn2, G, Gtc, Gn, Gp, Rh, T, stat=err)
   deallocate(G1, G2, G3, G4, G5, Tr, stat=err)
   deallocate(NN, Ud, No, IT, I1, I2, I3, stat=err) 
   deallocate(A, Gf, f_e, Rh1, Ua, stat=err)
   deallocate(SigOut1, SigOut2, SigInp, SigOutp,Gamp, stat=err)
   deallocate(IV, MV, stat=err)

   write(var,*) 'NEGF Matricies Deallocate'
   call Errors_Allocate(err,var)

end subroutine Init_Deallocate

end module
