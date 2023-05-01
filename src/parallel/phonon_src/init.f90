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
   character*100 :: var

   !! Allocation Variables !!
   if(debug_l > 5 .and. rank .eq. 0) then
     allocate(A1(Np,NE), Gp1(Np,NE), Gf1(Np,NE), Tr1(NE,Np), Id1(NE,Np), Id2(NE,Np), Id3(NE,Np), stat=err)
     allocate(Si1(Np,NE,Np),Sit(Np,NE,Np), stat=err)
   end if

   if(rank .eq. 0) then
     allocate(Aall(Np,NE,Np),Gnall(Np,NE,Np), stat=err)
   end if

   if(debug_l > 5) then
     allocate(Si(Np,NE,Np), stat=err)
   end if

   if(psc .eq. 1) then !parallel scatter every rank needs variables
     allocate(nse(Np,NE), nsa(Np,NE), pse(Np,NE), psa(Np,NE), SigInpNew(Np,NE), SigOutpNew(Np,NE), stat=err)
     allocate(SigInpNewt(Np,NE), SigOutpNewt(Np,NE), stat=err)
   else if(rank .eq. 0) then
     allocate(nse(Np,NE), nsa(Np,NE), pse(Np,NE), psa(Np,NE), SigInpNew(Np,NE), SigOutpNew(Np,NE), stat=err)
     allocate(SigInpNewt(Np,NE), SigOutpNewt(Np,NE), stat=err)
   end if

   allocate(Sig1(Np,Np), Gam1(Np,Np), SigIn1(Np,Np), Sig2(Np,Np), Gam2(Np,Np), stat=err) ! Allocate Matrix Variables 
   allocate(SigIn2(Np,Np), G(Np,Np), Gtc(Np,Np), Gf(Np,NE), Gp(Np,NE), Gn(Np,NE), T(Np,Np), stat=err)
   allocate(G1(Np,Np), G2(Np,Np), G3(Np,Np), G4(Np,Np), G5(Np,Np), stat=err)
   allocate(IT(NE), I1(NE), I2(NE), I3(NE), Rh(Np), Tr(NE,2), Rh1(Np,Np), stat=err)
   allocate(A(Np,NE), Gp(Np,NE), Gf(Np,NE), stat=err)
   allocate(SigOut1(Np,Np), SigOut2(Np,Np), SigInp(Np,NE), SigOutp(Np,NE), Gamp(Np), stat=err)
   allocate(f_e(NE,3),Tprt(Np,NE),Tprtt(Np,NE), Gn_ee(Np,Np), Tprt_NE(Np,NE), stat=err)
   write(var,*) 'NEGF Matricies Allocate'
   call Errors_Allocate(err,var)

end subroutine Init_Allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Init_Var - Allocates arrays for NEGF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Init_Var

  !! Allocation Variables !!
  if(debug_l > 5 .and. rank .eq. 0) then
    Gp1(:,:) = 0.0; Gf1(:,:) = 0.0
    Id1(:,:)= 0.0; Id2(:,:)= 0.0; Id3(:,:) = 0.0; Tr1(:,:) = 0.0; A1(:,:) = 0.0
    Si1(:,:,:)=0.0; Sit(:,:,:)=0.0
  end if

   if(rank .eq. 0) then
     Aall(:,:,:) = 0.0; Gnall(:,:,:) = 0.0;
   end if

   if(debug_l > 5) then
     Si(:,:,:)=0.0
   end if

  if(psc .eq. 1) then 
    nse(:,:)=0.0; nsa(:,:)=0.0; pse(:,:)=0.0; psa(:,:)=0.0; SigInpNew(:,:)=0.0; SigOutpNew(:,:) = 0.0
    SigInpNewt(:,:)=0.0; SigOutpNewt(:,:)=0.0
  else if(rank .eq. 0) then
    nse(:,:)=0.0; nsa(:,:)=0.0; pse(:,:)=0.0; psa(:,:)=0.0; SigInpNew(:,:)=0.0; SigOutpNew(:,:) = 0.0
    SigInpNewt(:,:)=0.0; SigOutpNewt(:,:)=0.0
  end if

  G(:,:) = (0.0,0.0); Gtc(:,:) = (0.0,0.0)
  G1(:,:) = (0.0,0.0); G2(:,:) = (0.0,0.0); G3(:,:) = (0.0,0.0); G4(:,:) = (0.0,0.0); G5(:,:) = (0.0,0.0)
  Gn(:,:) = (0.0,0.0);
  Sig1(:,:) = (0.0,0.0); Sig2(:,:) = (0.0,0.0)
  SigIn1(:,:) = (0.0,0.0); SigIn2(:,:) = (0.0,0.0)
  Gam1(:,:) = (0.0,0.0); Gam2(:,:) = (0.0,0.0)
  Tr(:,:) = 0.0; A(:,:) = 0.0; Rh1(:,:) = 0.0
  I1(:) = 0.0; I2(:) = 0.0; I3(:) = 0.0
  Gp(:,:) = 0.0; Gf(:,:) = 0.0
  Rh(:) = 0.0; Tprt(:,:) = 0.0; Tprtt(:,:) = 0.0; Gn_ee(:,:) = 0.0
  SigOut1(:,:)=0.0; SigOut2(:,:)=0.0; SigInp(:,:)=0.0; SigOutp(:,:)=0.0; Gamp(:)=0.0; Tprt_NE(:,:)=0.0 

end subroutine Init_Var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Init_Deallocate - Deallocates arrays
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Init_Deallocate

    !! Deallocate Dynamic Memory !!
   integer(kind=4) :: err
   character*100 :: var

   !! Allocation Variables !!
   if(debug_l > 5 .and. rank .eq. 0) then
     deallocate(A1, Gp1, Gf1, Tr1, Id1, Id2, Id3, stat=err)
     deallocate(Si1, Sit, stat=err)
   end if

   if(rank .eq. 0) then
     deallocate(Aall,Gnall, stat=err)
   end if

   if(debug_l > 5) then
     deallocate(Si, stat=err)
   end if

   if(psc .eq. 1) then 
     deallocate(nse, nsa, pse, psa, SigInpNew, SigOutpNew, stat=err)
     deallocate(SigInpNewt, SigOutpNewt, stat=err)
   else if(rank .eq. 0) then
     deallocate(nse, nsa, pse, psa, SigInpNew, SigOutpNew, stat=err)
     deallocate(SigInpNewt, SigOutpNewt, stat=err)
   end if

   deallocate(Sig1, Gam1, SigIn1, Sig2, Gam2, stat=err) ! Allocate Matrix Variables 
   deallocate(SigIn2, G, Gtc, Gf, Gp, Gn, T, stat=err)
   deallocate(G1, G2, G3, G4, G5, stat=err)
   deallocate(IT, I1, I2, I3, Rh, Rh1, Tr, stat=err)
   deallocate(A, Gp, Gf, stat=err)
   deallocate(SigOut1, SigOut2, SigInp, SigOutp, Gamp, stat=err)
   deallocate(f_e, Tprt, Tprtt, Gn_ee, Tprt_NE, stat=err)

   write(var,*) 'NEGF Matricies Dellocate'
   call Errors_Allocate(err,var)

end subroutine Init_Deallocate

end module
