!  interpolate.f90 
!
!  Authors: G. Walker - Vanderbilt University
!           T. Musho - Vanderbilt University
!
!  Created: 9/20/10
!  Last Modified: 9/20/10
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
!  PROGRAM: interpolate(Module)
!
!  PURPOSE: Contains interpolatation routines to interpolate between phonon
!           and electron energy and spatial coords
!
!  SUBROUTINES: 
!               
!                
!
!****************************************************************************
!
module Interpolate

use Startup, ONLY : pi,q,tst
use Debug
use Errors

implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Interpolate - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Interpolate_ph2el(GG,Gph,NE,NE_p,Np,Np_p,Nb,Nb_p,Ef,Ef_p,Eo,Eo_p,dE,dE_p,Lt,debug_l)

   integer(kind=4) :: r, j, q, icp, err
   integer(kind=4), intent(in) :: NE, Np, Nb
   integer(kind=4), intent(in) :: NE_p, Np_p, Nb_p
   integer(kind=4), intent(in) :: debug_l
   real(kind=8) :: Ef, Ef_p, Eo, Eo_p, dE, dE_p,Lt
   real(kind=8), allocatable, intent(in), dimension(:,:,:) :: Gph
   real(kind=8), allocatable, intent(out), dimension(:,:,:) :: GG

   integer(kind=4) :: NE_n, NE_s, Np_s
   real(kind=8) :: NE_bin, Np_bin, NE_bin_r, Np_bin_r
   character*120 :: var

   if(debug_l>5)write(6,*)'** Interpolating Phonon Spectra with Electron Spectra'

   if(debug_l>5)write(6,*)' Phonon Energy Bins, Electron Energy Bins = ',NE_p, NE
   if(debug_l>5)write(6,'(1X,A,ES12.5,ES12.5)')' Phonon Energy Range (Eo,Ef) = ',Eo_p,Ef_p
   if(debug_l>5)write(6,'(1X,A,ES12.5,ES12.5)')' Electron Energy Range (Eo,Ef) = ',Eo,Ef
   if(debug_l>5)write(6,*)' Phonon Spatial Bins, Electron Spatial Bins = ',Np_p, Np

   NE_bin = real(NE_p)/real(NE)*(Ef-Eo)/(Ef_p-Eo_p) ! Number of phonon energy bins per electron energy bins
   Np_bin = real(Np_p)/real(Np)

   !calc remander
   if(NE_bin .gt. 1.0) then
     NE_bin_r = NE_bin - floor(NE_bin)
   else
     NE_bin_r = 1/NE_bin - floor(1/NE_bin)
   end if
   if(Np_bin .gt. 1.0) then
     Np_bin_r = Np_bin - floor(Np_bin)
   else
     Np_bin_r = 1/Np_bin - floor(1/Np_bin)
   end if

   if(debug_l>5)write(6,'(1X,A,ES12.5,ES12.5)')' Phonon Energy Bins per Electron Energy Bins = ',NE_bin, NE_bin_r
   if(debug_l>5)write(6,'(1X,A,ES12.5,ES12.5)')' Phonon Spatial Bins per Electron Spatial Bins = ',Np_bin, Np_bin_r
   if(debug_l>5)write(6,'(1X,A,ES12.5,ES12.5)')' Phonon dE, Electron dE = ',dE_p, dE

   allocate(GG(NE,Np,Nb_p), stat=err) ! temp matrix
   write(var,*) 'GG(NE,Np,Nb_p)'
   call Errors_Allocate(err,var) !Check Errors
   GG(:,:,:) = 0.0

   NE_n = ceiling((Ef_p-Eo_p)/dE)
   NE_s = floor(((Eo_p+NE_bin*dE_p)-Eo)/dE)
   Np_s = floor(Np_bin/2)
   if(debug_l>5)write(6,*)' Phonon Energy Range in Electron Bins = ',NE_n
   if(debug_l>5)write(6,*)' Phonon Energy Starting Electron Bin = ',NE_s
   if(debug_l>5)write(6,*)' Phonon Energy Starting Position Bin = ',Np_s

   !Energy conservation remesh
   if(NE_bin .gt. 1.0 .and. Np_bin .gt. 1.0) icp = 0 !NE_p>NE Np_p>Np 
   if(NE_bin .lt. 1.0 .and. Np_bin .gt. 1.0) icp = 1 !NE_p<NE Np_p>Np
   if(NE_bin .gt. 1.0 .and. Np_bin .lt. 1.0) icp = 2 !NE_p>NE Np_p<Np
   if(NE_bin .lt. 1.0 .and. Np_bin .lt. 1.0) icp = 3 !NE_p<NE Np_p<Np

   !Spline remesh
   if(NE_bin .gt. 1.0 .and. Np_bin .lt. 1.0) icp = 4 !Cubic hermite spline

   select case (icp)
     case (0)
       if(debug_l>5)write(6,*) ' Interpolate Choice 1 - Conservation of Energy'
       if(tst .eq. 2 .or. tst .eq. 3) then !these two scattering methods don't care about how the energy is distributed
         NE_bin = 1                        !so by making 1bin in phonon energy = 1bin electron energy better interpolation
         call Interpolate_e1p1(GG,Gph,NE_n,NE_s,Np,Np_s,Nb_p,NE_bin,Np_bin,NE_bin_r,Np_bin_r) !NE_p>NE Np_p>Np
       else
         call Interpolate_e1p1(GG,Gph,NE_n,NE_s,Np,Np_s,Nb_p,NE_bin,Np_bin,NE_bin_r,Np_bin_r) !NE_p>NE Np_p>Np
       end if
     case (1)
       if(debug_l>5)write(6,*) ' Interpolate Choice 2 - Conservation of Energy'
       if(tst .eq. 2 .or. tst .eq. 3) then
         NE_bin = ceiling(real(NE)/real(NE_p))
         call Interpolate_e2p1(GG,Gph,NE_n,NE_s,Np,Np_s,Nb_p,NE_bin,Np_bin,NE_bin_r,Np_bin_r) !NE_p<NE Np_p>Np
       else
         call Interpolate_e2p1(GG,Gph,NE_n,NE_s,Np,Np_s,Nb_p,NE_bin,Np_bin,NE_bin_r,Np_bin_r) !NE_p<NE Np_p>Np
       end if
     case (2)
       if(debug_l>5)write(6,*) ' Interpolate Choice 3 - Conservation of Energy'
       if(tst .eq. 2 .or. tst .eq. 3) then
         NE_bin = 1
         call Interpolate_e1p2(GG,Gph,NE_n,NE_s,Np,Np_s,Nb_p,NE_bin,Np_bin,NE_bin_r,Np_bin_r) !NE_p>NE Np_p<Np
       else
         call Interpolate_e1p2(GG,Gph,NE_n,NE_s,Np,Np_s,Nb_p,NE_bin,Np_bin,NE_bin_r,Np_bin_r) !NE_p>NE Np_p<Np
       end if
     case (3)
       if(debug_l>5)write(6,*) ' Interpolate Choice 4 - Conservation of Energy'
       if(tst .eq. 2 .or. tst .eq. 3) then
         NE_bin = ceiling(real(NE)/real(NE_p))
         call Interpolate_e2p2(GG,Gph,NE_n,NE_s,Np,Np_s,Nb_p,NE_bin,Np_bin,NE_bin_r,Np_bin_r) !NE_p<NE Np_p<Np
       else
         call Interpolate_e2p2(GG,Gph,NE_n,NE_s,Np,Np_s,Nb_p,NE_bin,Np_bin,NE_bin_r,Np_bin_r) !NE_p<NE Np_p<Np
       end if
     case (4)
       if(debug_l>5)write(6,*) ' Interpolate Choice 5 - Cubic Hermit Spline Interpolation'
       if(tst .eq. 2 .or. tst .eq. 3) then
         if(NE_p .gt. NE) then
         NE_n = NE; NE_s =0
         else
         NE_n = NE_p; NE_s =0
         end if
         call Interpolate_chs2(GG,Gph,NE,NE_n,NE_p,NE_s,Np,Np_p,Np_s,Nb_p,NE_bin,Np_bin,Ef,Ef_p,Eo,Eo_p,Lt) !NE_p<NE Np_p<Np
       else
         call Interpolate_chs(GG,Gph,NE,NE_n,NE_p,NE_s,Np,Np_p,Np_s,Nb_p,NE_bin,Np_bin,Ef,Ef_p,Eo,Eo_p,Lt) !NE_p<NE Np_p<Np
       end if
   end select

   call Interpolate_check(GG,Gph,dE,dE_p,Np,Np_p,NE,NE_p)

   if(debug_l>5)write(6,*)' Interpolation Done'

end subroutine Interpolate_ph2el

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Interpolate_check - compare mesh energies
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Interpolate_check(GG,Gph,dE,dE_p,Np,Np_p,NE,NE_p)

   integer(kind=4), intent(in) :: Np, NE, Np_p, NE_p
   real(kind=8), intent(in) :: dE,dE_p
   real(kind=8), allocatable, intent(in), dimension(:,:,:) :: Gph
   real(kind=8), allocatable, intent(in), dimension(:,:,:) :: GG

   integer(kind=4) :: n
   real(kind=8) :: m_old, m_new, mr_old, mr_new

   if(debug_l>5)write(6,*)'** Cheking Interpolated Spectra'

   m_new = 0.0; m_old = 0.0
   do n = 1, Np
     m_new = m_new + sum(GG(:,n,1))*dE !integrate
   end do
   do n = 1, Np_p
     m_old = m_old + sum(Gph(:,n,1))*dE_p !integrate
   end do

   !mr_new = 0.0; mr_old = 0.0
   !do n = 1, NE
   !  mr_new = mr_new + GG(n,:,1)
   !end do
   !do n = 1, NE_p
   !  mr_old = mr_old + Gph(n,:,1)
   !end do

   !write(6,*)' Spectral Energy - Original Mesh (col) = ',m_old
   !write(6,*)' Spectral Energy - Interpolated Mesh (col) = ',m_new
   !write(6,*)' Spectral Energy Percent Difference (col) =',abs(m_old-m_new)/m_old*100,'%'

   !write(6,*)' Spectral Energy - Original Mesh (row) = ',mr_old
   !write(6,*)' Spectral Energy - Interpolated Mesh (row) = ',mr_new
   !write(6,*)' Spectral Energy Percent Difference (row) =',abs(mr_old-mr_new)/mr_old*100,'%'

   !m_old = m_old + mr_old; m_new = m_new + mr_new;
   if(debug_l>5)write(6,*)' Spectral Energy - Original Mesh (total) = ',m_old
   if(debug_l>5)write(6,*)' Spectral Energy - Interpolated Mesh (total) = ',m_new
   if(debug_l>5)write(6,*)' Spectral Energy Percent Difference (total) =',abs(m_old-m_new)/m_old*100,'%'

end subroutine Interpolate_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Interpolate_e1p1 - NE_p>NE Np_p>Np - Choice 1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Interpolate_e1p1(GG,Gph,NE,NE_s,Np,Np_s,Nb_p,NE_bin,Np_bin,NE_bin_r,Np_bin_r)

   integer(kind=4), intent(in) :: NE, Np, Nb_p, NE_s, Np_s
   real(kind=8), intent(in) :: NE_bin, Np_bin, NE_bin_r, Np_bin_r
   real(kind=8), allocatable, intent(in), dimension(:,:,:) :: Gph
   real(kind=8), allocatable, intent(inout), dimension(:,:,:) :: GG

   integer(kind=4) :: r, j, q
   integer(kind=4) :: qq, q_p, jj, j_p

   do r = 1, Nb_p
     do j = 1+NE_s, NE+NE_s
       j_p = ((j-NE_s-1)-1)*ceiling(NE_bin)+1
         !write(*,*) 'j_p = ',j_p
         do jj = 0, ceiling(NE_bin)-1
           !write(*,*) 'jj = ',jj
           if(jj .lt. ceiling(NE_bin)-1 .or. NE_bin_r .eq. 0.0) then
             do q = 1, Np
               q_p = (q-1)*ceiling(Np_bin)+1
               !write(*,*) 'q_p = ',q_p
               do qq = 0, ceiling(Np_bin)-1
                 if(qq .lt. ceiling(Np_bin)-2 .or. Np_bin_r .eq. 0.0) then
                   GG(j,q,r) = GG(j,q,r)+Gph(j_p+jj,q_p+qq,r)
                 else
                   GG(j,q,r) = GG(j,q,r)+Gph(j_p+jj,q_p+qq,r)*Np_bin_r
                   if(q .gt. 1) then
                     GG(j,q,r) = GG(j,q,r)+Gph(j_p,q_p-1,r)*(1-Np_bin_r)
                   end if
                 end if
               end do
             end do
           else
             do q = 1, Np
               q_p = (q-1)*ceiling(Np_bin)+1
               !write(*,*) 'q_p2 = ',q_p
               do qq = 0, ceiling(Np_bin)-1
                 !write(*,*) 'qq2 = ',qq
                 if(qq .lt. ceiling(Np_bin)-2 .or. Np_bin_r .eq. 0.0) then
                   GG(j,q,r) = GG(j,q,r)+Gph(j_p+jj,q_p+qq,r)*NE_bin_r
                 else
                   GG(j,q,r) = GG(j,q,r)+Gph(j_p+jj,q_p+qq,r)*Np_bin_r*NE_bin_r
                   if(q .gt. 1) then
                     GG(j,q,r) = GG(j,q,r)+Gph(j_p,q_p-1,r)*(1-Np_bin_r)*NE_bin_r
                   end if
                 end if
               end do
             end do
             if(j .gt. 1) then
               GG(j,q,r) = GG(j,q,r)+Gph(j_p-1,q_p-1,r)*(1-NE_bin_r)*NE_bin_r
             end if
           end if
         end do
     end do
   end do

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Interpolate_e2p1 - NE_p<NE Np_p>Np - Choice 2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Interpolate_e2p1(GG,Gph,NE,NE_s,Np,Np_s,Nb_p,NE_bin,Np_bin,NE_bin_r,Np_bin_r)

   integer(kind=4), intent(in) :: NE, Np, Nb_p, NE_s, Np_s
   real(kind=8), intent(in) :: NE_bin, Np_bin, NE_bin_r, Np_bin_r
   real(kind=8), allocatable, intent(in), dimension(:,:,:) :: Gph
   real(kind=8), allocatable, intent(inout), dimension(:,:,:) :: GG

   integer(kind=4) :: r, j, q
   integer(kind=4) :: qq, q_p, jj, j_p

   do r = 1, Nb_p
     do j = 1+NE_s, NE+NE_s
       j_p = floor((j-NE_s-1)*NE_bin)+1
         !write(*,*) 'j_p = ',j_p
         do jj = 0, ceiling(1/NE_bin)
           !write(*,*) 'jj = ',jj
             do q = 1, Np
               q_p = (q-1)*ceiling(Np_bin)+1
               !write(*,*) 'q_p = ',q_p
               do qq = 0, ceiling(Np_bin)-1
                 if(qq .lt. ceiling(Np_bin)-2 .or. Np_bin_r .eq. 0.0) then
                 !write(*,*) 'G(j,q,r) = ',j+jj,q,r
                 !write(*,*) 'Gph(j,q,r) = ',j_p,q_p+qq,r
                   GG(j+jj,q,r) = GG(j+jj,q,r)+Gph(j_p,q_p+qq,r)
                 else
                   GG(j+jj,q,r) = GG(j+jj,q,r)+Gph(j_p,q_p+qq,r)*Np_bin_r
                   if(q .gt. 1) then
                     GG(j+jj,q,r) = GG(j+jj,q,r)+Gph(j_p,q_p-1,r)*(1-Np_bin_r)
                   end if
                 end if
               end do
             end do
         end do
     end do
   end do

end subroutine Interpolate_e2p1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Interpolate_e1p2 - NE_p>NE Np_p<Np - Choice 3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Interpolate_e1p2(GG,Gph,NE,NE_s,Np,Np_s,Nb_p,NE_bin,Np_bin,NE_bin_r,Np_bin_r)

   integer(kind=4), intent(in) :: NE, Np, Nb_p, NE_s, Np_s
   real(kind=8), intent(in) :: NE_bin, Np_bin, NE_bin_r, Np_bin_r
   real(kind=8), allocatable, intent(in), dimension(:,:,:) :: Gph
   real(kind=8), allocatable, intent(inout), dimension(:,:,:) :: GG

   integer(kind=4) :: r, j, q
   integer(kind=4) :: qq, q_p, jj, j_p
   do r = 1, Nb_p
     do j = 1+NE_s, NE+NE_s
       j_p = floor((j-NE_s-1)*(NE_bin+1))+1
         !write(*,*) 'NE_bin =', NE_bin
         !write(*,*) 'j_p = ',j_p
         do jj = 0, floor(NE_bin)
           !write(*,*) 'jj = ',jj
           if(j .gt. 1 .and. jj .eq. 0) then
             do q = 1, Np-Np_s
               q_p = floor(q*Np_bin)+1
               if(q_p .lt. floor((q+1)*Np_bin)+1 .and. q .ne. Np) then
                 GG(j,q+Np_s,r) = GG(j,q+Np_s,r) + Gph(j_p,q_p,r)*Np_bin*(1-NE_bin_r)
                 GG(j,q+Np_s+1,r) = GG(j,q+Np_s+1,r) + Gph(j_p,q_p,r)*Np_bin*Np_bin_r*(1-NE_bin_r)
               else if(q_p .gt. floor((q-1)*Np_bin)+1 .and. q .ne. 1) then
                 GG(j,q+Np_s,r) = GG(j,q+Np_s,r) + Gph(j_p,q_p,r)*Np_bin*(1-Np_bin_r)*(1-NE_bin_r)
               else
                 GG(j,q+Np_s,r) = GG(j,q+Np_s,r) + Gph(j_p,q_p,r)*Np_bin*(1-NE_bin_r)
               end if
               !end do
             end do
           else if(jj .eq. floor(NE_bin)) then
             do q = 1, Np-Np_s
               q_p = floor(q*Np_bin)+1
               if(q_p .lt. floor((q+1)*Np_bin)+1 .and. q .ne. Np) then
                 GG(j,q+Np_s,r) = GG(j,q+Np_s,r) + Gph(j_p+jj,q_p,r)*Np_bin*NE_bin_r
                 GG(j,q+Np_s+1,r) = GG(j,q+Np_s+1,r) + Gph(j_p+jj,q_p,r)*Np_bin*Np_bin_r*NE_bin_r
               else if(q_p .gt. floor((q-1)*Np_bin)+1 .and. q .ne. 1) then
                 GG(j,q+Np_s,r) = GG(j,q+Np_s,r) + Gph(j_p+jj,q_p,r)*Np_bin*(1-Np_bin_r)*NE_bin_r
               else
                 GG(j,q+Np_s,r) = GG(j,q+Np_s,r) + Gph(j_p+jj,q_p,r)*Np_bin*NE_bin_r
               end if
             end do
           else
             do q = 1, Np-Np_s
               q_p = floor(q*Np_bin)+1 !Np_p<Np
               !write(*,*) 'Np_bin =', Np_bin
               !write(*,*) 'q,q_p,Np_bin_r = ',q,q_p,Np_bin_r
               !do qq = 0, ceiling(1/Np_bin)
               !write(*,*) 'G(j,q,r) = ',j,q,r
               !write(*,*) 'Gph(j,q,r) = ',j_p+jj,q_p,r
               if(q_p .lt. floor((q+1)*Np_bin)+1 .and. q .ne. Np) then
                 GG(j,q+Np_s,r) = GG(j,q+Np_s,r) + Gph(j_p+jj,q_p,r)*Np_bin
                 GG(j,q+Np_s+1,r) = GG(j,q+Np_s+1,r) + Gph(j_p+jj,q_p,r)*Np_bin*Np_bin_r
               else if(q_p .gt. floor((q-1)*Np_bin)+1 .and. q .ne. 1) then
                 GG(j,q+Np_s,r) = GG(j,q+Np_s,r) + Gph(j_p+jj,q_p,r)*Np_bin*(1-Np_bin_r)
               else
                 GG(j,q+Np_s,r) = GG(j,q+Np_s,r) + Gph(j_p+jj,q_p,r)*Np_bin
               end if
               !end do
             end do
           end if
         end do
     end do
   end do

end subroutine Interpolate_e1p2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Interpolate_e1p1 - NE_p<NE Np_p<Np - Choice 4
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Interpolate_e2p2(GG,Gph,NE,NE_s,Np,Np_s,Nb_p,NE_bin,Np_bin,NE_bin_r,Np_bin_r)

   integer(kind=4), intent(in) :: NE, Np, Nb_p, NE_s, Np_s
   real(kind=8), intent(in) :: NE_bin, Np_bin, NE_bin_r, Np_bin_r
   real(kind=8), allocatable, intent(in), dimension(:,:,:) :: Gph
   real(kind=8), allocatable, intent(inout), dimension(:,:,:) :: GG

   integer(kind=4) :: r, j, q
   integer(kind=4) :: qq, q_p, jj, j_p

   do r = 1, Nb_p
     do j = 1+NE_s, NE+NE_s
       j_p = floor((j-NE_s-1)*NE_bin)+1
       !write(*,*) 'j_p = ',j_p
       !do jj = 0, ceiling(1/NE_bin)
             do q = 1, Np
               q_p = floor(q*Np_bin)+1
               !write(*,*) 'q_p = ',q_p
               !do qq = 0, ceiling(1/Np_bin)
                 !write(*,*) 'G(j,q,r) = ',j+jj,q+qq,r
                 !write(*,*) 'Gph(j,q,r) = ',j_p,q_p,r
                 GG(j,q,r) = GG(j,q,r) + Gph(j_p,q_p,r)*Np_bin*NE_bin
               !end do
             end do
       !end do
     end do
   end do

end subroutine Interpolate_e2p2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Interpolate_chs - 1D cubic hermite spline - Choice 5
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Interpolate_chs(GG,Gph,NE,NE_n,NE_p,NE_s,Np,Np_p,Np_s,Nb_p,NE_bin,Np_bin,Ef,Ef_p,Eo,Eo_p,Lt)

   integer(kind=4), intent(in) :: NE, NE_n, NE_p, Np, Np_p, Nb_p, NE_s, Np_s
   real(kind=8), intent(in) :: NE_bin, Np_bin,Ef,Ef_p,Eo,Eo_p,Lt
   real(kind=8), allocatable, intent(in), dimension(:,:,:) :: Gph
   real(kind=8), allocatable, intent(inout), dimension(:,:,:) :: GG

   integer(kind=4) :: r, j, q, incr=1, err=0
   integer(kind=4) :: qq, q_p, jj, j_p
   logical :: skip = .FALSE.
   real(kind=8) :: Ls = 0.0
   real(kind=8), allocatable, dimension(:) :: Ee, Ep, Np_ee, Np_pp, Ee_p
   real(kind=8), allocatable, dimension(:) :: D, GG_i_t
   real(kind=8), allocatable, dimension(:,:,:) :: GG_i
   character*100 :: var

   allocate(Ee(NE),Ep(NE_p),Np_ee(Np),Np_pp(Np_p),GG_i(NE,Np_p,Nb_p), stat=err) ! Need intermediate array
   write(var,*) 'Ee(NE),Ep(NE_p),Np_ee(Np),Np_pp(Np_p),GG_i(NE,Np_p,Nb_p)'
   call Errors_Allocate(err,var) !Check Errors
   allocate(D(NE_p), stat=err) ! Need intermediate array
   write(var,*) 'D(NE_p)'
   call Errors_Allocate(err,var) !Check Errors
   allocate(GG_i_t(NE_n),Ee_p(NE_n), stat=err) ! Need intermediate array
   write(var,*) 'GG_i_t(NE_n),Ee_p(NE_n)'
   call Errors_Allocate(err,var) !Check Errors
   GG_i(:,:,:) = 0.0; GG_i_t(:) = 0.0

   call linspace(Ee,NE,Eo,Ef) ! Specifying electron range
   do q=1+NE_s, NE_n+NE_s
     Ee_p(q-NE_s) = Ee(q) !use exact node energy values
   end do
   !call linspace(Ee_p,NE_n,Ee(NE_s),Ee(NE_s+NE_n)) ! Specifying electron grid phonon range
   call linspace(Ep,NE_p,Eo_p,Ef_p) ! Specifying phonon range

   !SLACTEC Library calls
   if(debug_l>5)write(6,*) '  Starting Energy Sweep Interpolation'
   if(debug_l>5)write(6,*) '    Interpolating from ',NE_p,'->',NE_n
   !Interpolate Energy Sweep
   do r = 1, Nb_p
     do j = 1, Np_p
       !call DPCHIM(N,X,F,D,INCFD,IERR) !Calculate hermite spline derivatives
       call DPCHIM(NE_p,Ep,Gph(:,j,r),D,incr,err)
       call Debug_DPCHIM(err)
       !call DPCHFE(N,X,F,D,INCFD,SKIP,NE,XE,FE,IERR) !Evaluate
       call DPCHFE(NE_p,Ep,Gph(:,j,r),D,incr,skip,NE_n,Ee_p,GG_i_t,err)
       call Debug_DPCHFE(err)
       !Wasn't able to extrapolate outside of spline so had to create temp array
       do q = 1+NE_s, NE_n+NE_s
         GG_i(q,j,r) = GG_i_t(q-NE_s) !puts array in matrix, matrix init to zero
       end do
     end do
   end do

   deallocate(D)
   allocate(D(Np_p), stat=err) ! Need intermediate array
   write(var,*) 'D(NE_p)'
   call Errors_Allocate(err,var) !Check Errors

   call linspace(Np_ee,Np,Ls,Lt)
   call linspace(Np_pp,Np_p,Ls,Lt)

   if(debug_l>5)write(6,*) '  Starting Spatial Sweep Interpolation'
   if(debug_l>5)write(6,*) '    Interpolating from ',Np_p,'->',Np
   !Interpolate Spatial Sweep
   do r = 1, Nb_p
     do j = 1+NE_s, NE_n+NE_s
       !call DPCHIM(N,X,F,D,INCFD,IERR) !Calculate hermite spline derivatives
       call DPCHIM(Np_p,Np_pp,GG_i(j,:,r),D,incr,err)
       call Debug_DPCHIM(err)
       !call DPCHFE(N,X,F,D,INCFD,SKIP,NE,XE,FE,IERR) !Evaluate
       call DPCHFE(Np_p,Np_pp,GG_i(j,:,r),D,incr,skip,Np,Np_ee,GG(j,:,r),err) 
       call Debug_DPCHFE(err)
     end do
   end do

   deallocate(GG_i,Ee,Ep,D)

end subroutine Interpolate_chs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Interpolate_chs2 - 1D cubic hermite spline - Choice 5a - Spatial Interp only
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Interpolate_chs2(GG,Gph,NE,NE_n,NE_p,NE_s,Np,Np_p,Np_s,Nb_p,NE_bin,Np_bin,Ef,Ef_p,Eo,Eo_p,Lt)

   integer(kind=4), intent(in) :: NE, NE_n, NE_p, Np, Np_p, Nb_p, NE_s, Np_s
   real(kind=8), intent(in) :: NE_bin, Np_bin,Ef,Ef_p,Eo,Eo_p,Lt
   real(kind=8), allocatable, intent(in), dimension(:,:,:) :: Gph
   real(kind=8), allocatable, intent(inout), dimension(:,:,:) :: GG

   integer(kind=4) :: r, j, q, incr=1, err=0
   integer(kind=4) :: qq, q_p, jj, j_p
   logical :: skip = .FALSE.
   real(kind=8) :: Ls = 0.0
   real(kind=8), allocatable, dimension(:) :: Ee, Ep, Np_ee, Np_pp, Ee_p
   real(kind=8), allocatable, dimension(:) :: D, GG_i_t
   real(kind=8), allocatable, dimension(:,:,:) :: GG_i
   character*100 :: var

   allocate(Ee(NE),Ep(NE_p),Np_ee(Np),Np_pp(Np_p),GG_i(NE,Np_p,Nb_p), stat=err) ! Need intermediate array
   write(var,*) 'Ee(NE),Ep(NE_p),Np_ee(Np),Np_pp(Np_p),GG_i(NE,Np_p,Nb_p)'
   call Errors_Allocate(err,var) !Check Errors
   allocate(D(NE_p), stat=err) ! Need intermediate array
   write(var,*) 'D(NE_p)'
   call Errors_Allocate(err,var) !Check Errors
   allocate(GG_i_t(NE_n),Ee_p(NE_n), stat=err) ! Need intermediate array
   write(var,*) 'GG_i_t(NE_n),Ee_p(NE_n)'
   call Errors_Allocate(err,var) !Check Errors
   GG_i(:,:,:) = 0.0; GG_i_t(:) = 0.0

   call linspace(Ee,NE,Eo,Ef) ! Specifying electron range
   do q=1+NE_s, NE_n+NE_s
     Ee_p(q-NE_s) = Ee(q) !use exact node energy values
   end do
   !call linspace(Ee_p,NE_n,Ee(NE_s),Ee(NE_s+NE_n)) ! Specifying electron grid phonon range
   call linspace(Ep,NE_p,Eo_p,Ef_p) ! Specifying phonon range

   !SLACTEC Library calls
   if(debug_l>5)write(6,*) '  Starting Energy Sweep Interpolation'
   if(debug_l>5)write(6,*) '    Interpolating from ',NE_p,'->',NE_n
   !Interpolate Energy Sweep
   do r = 1, Nb_p
     do j = 1, Np_p
       GG_i(1:NE_n,j,r)=Gph(1:NE_n,j,r)
     end do
   end do

   deallocate(D)
   allocate(D(Np_p), stat=err) ! Need intermediate array
   write(var,*) 'D(NE_p)'
   call Errors_Allocate(err,var) !Check Errors

   call linspace(Np_ee,Np,Ls,Lt)
   call linspace(Np_pp,Np_p,Ls,Lt)

   if(debug_l>5)write(6,*) '  Starting Spatial Sweep Interpolation'
   if(debug_l>5)write(6,*) '    Interpolating from ',Np_p,'->',Np
   !Interpolate Spatial Sweep
   do r = 1, Nb_p
     do j = 1, NE_n
       !call DPCHIM(N,X,F,D,INCFD,IERR) !Calculate hermite spline derivatives
       call DPCHIM(Np_p,Np_pp,GG_i(j,:,r),D,incr,err)
       call Debug_DPCHIM(err)
       !call DPCHFE(N,X,F,D,INCFD,SKIP,NE,XE,FE,IERR) !Evaluate
       call DPCHFE(Np_p,Np_pp,GG_i(j,:,r),D,incr,skip,Np,Np_ee,GG(j,:,r),err) 
       call Debug_DPCHFE(err)
     end do
   end do

   deallocate(GG_i,Ee,Ep,D)

end subroutine Interpolate_chs2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Interpolate_chs2d - 2D cubic hermite spline - Choice 6
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Interpolate_chs2d(GG,Gph,NE,NE_p,NE_s,Np,Np_p,Np_s,Nb_p,NE_bin,Np_bin,Ef,Ef_p,Eo,Eo_p,Lt)

   integer(kind=4), intent(in) :: NE, NE_p, Np, Np_p, Nb_p, NE_s, Np_s
   real(kind=8), intent(in) :: NE_bin, Np_bin,Ef,Ef_p,Eo,Eo_p,Lt
   real(kind=8), allocatable, intent(in), dimension(:,:,:) :: Gph
   real(kind=8), allocatable, intent(inout), dimension(:,:,:) :: GG

   integer(kind=4) :: r, j, q, incr=1, err=0
   integer(kind=4) :: qq, q_p, jj, j_p
   logical :: skip = .false.
   real(kind=8) :: Ls = 0.0
   real(kind=8), allocatable, dimension(:) :: Ee, Ep, Np_ee, Np_pp
   real(kind=8), allocatable, dimension(:) :: D
   real(kind=8), allocatable, dimension(:,:,:) :: GG_i
   character*100 :: var

   allocate(Ee(NE),Ep(NE_p),Np_ee(Np),Np_pp(Np_p),GG_i(NE,Np_p,Nb_p), stat=err) ! Need intermediate array
   write(var,*) 'Ee(NE),Ep(NE_p),Np_ee(Np),Np_pp(Np_p),GG_i(NE,Np_p,Nb_p)'
   call Errors_Allocate(err,var) !Check Errors
   allocate(D(NE_p*Np_p), stat=err) ! Need intermediate array
   write(var,*) 'D(NE_p)'
   call Errors_Allocate(err,var) !Check Errors
   GG_i(:,:,:) = 0.0

   call linspace(Ee,NE,Eo,Ef) ! Specifying incoming electron range
   call linspace(Ep,NE_p,Eo_p,Ef_p) ! Specifying incoming phonon range

   !FIXME
   !SLACTEC Library calls
   !Interpolate Energy Sweep
   do r = 1, Nb_p
     !call DPCHIM(N,X,F,D,INCFD,IERR) !Calculate hermite spline derivatives
     call DPCHIM(NE_p,Ep,Gph(:,:,r),D,Np_p,err)
     call Debug_DPCHIM(err)
     !call DPCHFE(N,X,F,D,INCFD,SKIP,NE,XE,FE,IERR) !Evaluate
     call DPCHFE(NE_p,Ep,Gph(:,:,r),D,Np_p,skip,NE,Ee,GG(:,:,r),err)
     call Debug_DPCHFE(err)
   end do

   deallocate(GG_i,Ee,Ep,D)

end subroutine Interpolate_chs2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Interpolate_chs2 - 1D cubic hermite spline - Choice 5a - Spatial Interp only
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Interpolate_ph1d(GG,Gph,Np,Np_p,Nb_p,Lt)

   integer(kind=4), intent(in) :: Np, Np_p, Nb_p
   real(kind=8), intent(in) :: Lt
   real(kind=8), allocatable, intent(inout), dimension(:,:,:) :: Gph
   real(kind=8), allocatable, intent(in), dimension(:,:,:) :: GG

   integer(kind=4) :: r, j, q, incr=1, err=0
   integer(kind=4) :: qq, q_p, jj, j_p
   logical :: skip = .FALSE.
   real(kind=8) :: Ls = 0.0
   real(kind=8), allocatable, dimension(:) :: Np_ee, Np_pp
   real(kind=8), allocatable, dimension(:) :: D
   character*120 :: var

   allocate(Np_ee(Np), Np_pp(Np_p), stat=err) ! Need intermediate array
   write(var,*) 'Np_pp(Np_p)'
   call Errors_Allocate(err,var) !Check Errors

   allocate(D(Np), stat=err) ! Need intermediate array
   write(var,*) 'D(NE_p)'
   call Errors_Allocate(err,var) !Check Errors

   call linspace(Np_ee,Np,Ls,Lt)
   call linspace(Np_pp,Np_p,Ls,Lt)

   if(debug_l>5)write(6,*) '  Starting Spatial Sweep Interpolation'
   if(debug_l>5)write(6,*) '    Interpolating from ',Np,'->',Np_p
   !Interpolate Spatial Sweep
   do r = 1, Nb_p
       !call DPCHIM(N,X,F,D,INCFD,IERR) !Calculate hermite spline derivatives
       call DPCHIM(Np,Np_ee,GG(:,1,r),D,incr,err)
       call Debug_DPCHIM(err)
       !call DPCHFE(N,X,F,D,INCFD,SKIP,NE,XE,FE,IERR) !Evaluate
       call DPCHFE(Np,Np_ee,GG(:,1,r),D,incr,skip,Np_p,Np_pp,Gph(:,1,r),err) 
       call Debug_DPCHFE(err)
   end do

   deallocate(Np_ee,Np_pp,D)

end subroutine Interpolate_ph1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Interpolate_Plot - Plot Interpolatation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Interpolate_plot(GG,E,NE,Nb_p,Np,Gdx,fname)

   integer(kind=4) :: r, n, j, err
   real(kind=4) :: l
   character*120 :: var
   integer(kind=4), intent(in) :: Nb_p,NE,Np
   real(kind=8), allocatable, intent(in), dimension(:) :: E,Gdx
   real(kind=8), allocatable, intent(in), dimension(:,:,:) :: GG
   character*10, intent(in) :: fname

   open (unit=18, file=fname, status='NEW', action='write', iostat=err)
   write(var,*) fname
   call Errors_Fileopen(err,var)

   104 format (10000(ES12.5,1X))

   do r=1, NE
      l = 0.0
      do n=1,Np
        write(18,104) E(r),l,(GG(r,n,j)/(2.0*pi), j=1,Nb_p) !LDOS Contour Plot
        l = Gdx(n) + l
      end do
      write(18,*)
   end do

   close(18);

end subroutine Interpolate_plot

end module


