!  matxopt.f90 
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
!  PROGRAM: Potent(Module)
!
!  PURPOSE: Contains subroutine to operate matrix
!
!  SUBROUTINES: Matxopt_trace - Take trace of matrix
!               Matxopt_diag - Puts vector along matrix diag
!               
!                
!
!****************************************************************************
!
module Matxopt

!Use global variables from module
use Errors
use Debug

implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_ztcpy - Copy tridiagonal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_ztcpy(MM,TT,ns)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns
   complex(kind=8), dimension(:,:), intent(in) :: MM
   complex(kind=8), dimension(:,:), intent(inout) :: TT !in and out

   do n = 1, ns
      TT(n,n) = MM(n,n) ! Main diagonal
   end do

   do n = 1, ns-1
      TT(n,n+1) = MM(n,n+1) ! Main diagonal
   end do

   do n = 2, ns
      TT(n,n-1) = MM(n,n-1) ! Main diagonal
   end do

end subroutine Matxopt_ztcpy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_sub_zdtf2 - Subtract diagonal to tridiagonal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_sub_zdtf2(DD,TT,ns)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns
   complex(kind=8), dimension(:,:), intent(in) :: DD
   complex(kind=8), dimension(:,:), intent(inout) :: TT !in and out

   do n = 1, ns
      TT(n,n) = -DD(n,n)+TT(n,n) ! Main diagonal
   end do

end subroutine Matxopt_sub_zdtf2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_sub_zvtf2 - Subtract vector to tridiagonal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_sub_zvtf2(VV,TT,ns)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns
   complex(kind=8), dimension(:), intent(in) :: VV

   complex(kind=8), dimension(:,:), intent(inout) :: TT !in and out

   do n = 1, ns
      TT(n,n) = -VV(n)+TT(n,n) ! Main diagonal
   end do

end subroutine Matxopt_sub_zvtf2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_sub_zvtf - Subtract vector to tridiagonal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_sub_zvtf(VV,TT,ns)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns
   complex(kind=8), dimension(:), intent(in) :: VV
   complex(kind=8), dimension(:,:), intent(inout) :: TT !in and out

   do n = 1, ns
      TT(n,n) = VV(n)-TT(n,n) ! Main diagonal
   end do

end subroutine Matxopt_sub_zvtf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_sub_zstf2 - Subtract diagonal edges(1,ns) to tridiagonal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_sub_zstf2(DD,TT,ns)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns
   complex(kind=8), dimension(:,:), intent(in) :: DD
   complex(kind=8), dimension(:,:), intent(inout) :: TT !

   TT(1,1) = -DD(1,1)+TT(1,1) ! Main diagonal
   TT(ns,ns) = -DD(ns,ns)+TT(ns,ns) ! Main diagonal

end subroutine Matxopt_sub_zstf2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_sub_zdtf - Subtract diagonal to tridiagonal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_sub_zdtf(DD,TT,ns)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns
   complex(kind=8), dimension(:,:), intent(in) :: DD
   complex(kind=8), dimension(:,:), intent(inout) :: TT !in and out

   do n = 1, ns
      TT(n,n) = DD(n,n)-TT(n,n) ! Main diagonal
   end do

end subroutine Matxopt_sub_zdtf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_sub_zstf - Subtract diagonal edges(1,ns) to tridiagonal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_sub_zstf(DD,TT,ns)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns
   complex(kind=8), dimension(:,:), intent(in) :: DD
   complex(kind=8), dimension(:,:), intent(inout) :: TT !

   TT(1,1) = DD(1,1)-TT(1,1) ! Main diagonal
   TT(ns,ns) = DD(ns,ns)-TT(ns,ns) ! Main diagonal

end subroutine Matxopt_sub_zstf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_add_zdtf - Add diagonal to tridiagonal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_add_zdtf(DD,TT,ns)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns
   complex(kind=8), dimension(:,:), intent(in) :: DD
   complex(kind=8), dimension(:,:), intent(inout) :: TT !in and out

   do n = 1, ns
      TT(n,n) = DD(n,n)+TT(n,n) ! Main diagonal
   end do

end subroutine Matxopt_add_zdtf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_add_zstf - Add diagonal edges(1,ns) to tridiagonal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_add_zstf(DD,TT,ns)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns
   complex(kind=8), dimension(:,:), intent(in) :: DD
   complex(kind=8), dimension(:,:), intent(inout) :: TT !

   TT(1,1) = DD(1,1)+TT(1,1) ! Main diagonal
   TT(ns,ns) = DD(ns,ns)+TT(ns,ns) ! Main diagonal

end subroutine Matxopt_add_zstf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_trace - Take trace of matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_trace(MM,ns,V)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns
   real(kind=8), dimension(:,:), intent(in) :: MM
   real(kind=8), dimension(:), intent(out) :: V

   do n = 1, ns
      V(n) = MM(n,n) ! Stores diagonal elements
   end do

end subroutine Matxopt_trace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_trace_sub - Take trace off subdiagonal of matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_trace_sub(MM,df,ns,V)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns, df
   real(kind=8), dimension(:,:), intent(in) :: MM
   real(kind=8), dimension(:), intent(out) :: V

   do n = 1+df, ns
      V(n-df) = MM(n,n-df) ! Stores diagonal elements
   end do

end subroutine Matxopt_trace_sub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_trace_sup - Take trace off superdiagonal of matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_trace_sup(MM,df,ns,V)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns, df
   real(kind=8), dimension(:,:), intent(in) :: MM
   real(kind=8), dimension(:), intent(out) :: V

   do n = 1, (ns-df)
      V(n) = MM(n,n+df) ! Stores diagonal elements
   end do

end subroutine Matxopt_trace_sup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_trace - Take trace of matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_ztrace(MM,ns,V)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns
   complex(kind=8), dimension(:,:), intent(in) :: MM
   complex(kind=8), dimension(:), intent(out) :: V

   do n = 1, ns
      V(n) = MM(n,n) ! Stores diagonal elements
   end do

end subroutine Matxopt_ztrace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_ztrace_sub - Take trace off subdiagonal of matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_ztrace_sub(MM,df,ns,V)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns, df
   complex(kind=8), dimension(:,:), intent(in) :: MM
   complex(kind=8), dimension(:), intent(out) :: V

   do n = 1+df, ns
      V(n-df) = MM(n,n-df) ! Stores diagonal elements
   end do

end subroutine Matxopt_ztrace_sub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_ztrace_sup - Take trace off superdiagonal of matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_ztrace_sup(MM,df,ns,V)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns, df
   complex(kind=8), dimension(:,:), intent(in) :: MM
   complex(kind=8), dimension(:), intent(out) :: V

   do n = 1, (ns-df)
      V(n) = MM(n,n+df) ! Stores diagonal elements
   end do

end subroutine Matxopt_ztrace_sup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_tracesum - Take summation trace of matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_tracesum(MM,ns,R)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns
   real(kind=8), dimension(:,:), intent(in) :: MM
   real(kind=8), intent(out) :: R

   R = 0.0
   do n = 1, ns
      R = R + MM(n,n) ! Stores diagonal elements
   end do

end subroutine Matxopt_tracesum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_zput - Puts vector along matrix diag
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_zput(V,ns,MM)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns
   complex(kind=8), dimension(:), intent(in) :: V
   complex(kind=8), dimension(:,:), intent(out) :: MM

   do n = 1, ns
      MM(n,n) = V(n)
   end do

end subroutine Matxopt_zput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_zput_sub - Puts vector along matrix diag
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_zput_sub(V,df,ns,MM)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns, df
   complex(kind=8), dimension(:), intent(in) :: V
   complex(kind=8), dimension(:,:), intent(out) :: MM

   do n = 1+df, ns
     MM(n,n-df) = V(n-df)! Stores diagonal elements
   end do

end subroutine Matxopt_zput_sub

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_zput_sup - Puts vector along matrix diag
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_zput_sup(V,df,ns,MM)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns,df
   complex(kind=8), dimension(:), intent(in) :: V
   complex(kind=8), dimension(:,:), intent(out) :: MM

   do n = 1, ns-df
    MM(n,n+df) = V(n)! Stores diagonal elements
   end do

end subroutine Matxopt_zput_sup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_diag - Puts vector along matrix diag
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_diag(V,ns,MM)

   integer(kind=4) :: n
   integer(kind=4), intent(in) :: ns
   real(kind=8), dimension(:), intent(in) :: V
   real(kind=8), dimension(:,:), intent(out) :: MM

   do n = 1, ns
      MM(n,n) = V(n)
   end do

end subroutine Matxopt_diag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_ddff - (double) single value times full matrix, full matrix out
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_dsff(D,l,F,ns,MM)

   integer(kind=4) :: n, m
   integer(kind=4), intent(in) :: l, ns !l=equals location of single value
   real(kind=8), dimension(:,:), intent(in) :: D, F
   real(kind=8), dimension(:,:), intent(out) :: MM

   MM(:,:) = 0.0
   do n = 1, ns
     MM(n,l) = D(l,l) * F(n,l)
   end do

end subroutine Matxopt_dsff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_ddff - (double) diagonal times full matrix, full matrix out
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_ddff(D,F,ns,MM)

   integer(kind=4) :: n, m
   integer(kind=4), intent(in) :: ns
   real(kind=8), dimension(:,:), intent(in) :: D,F
   real(kind=8), dimension(:,:), intent(out) :: MM

   MM(:,:) = 0.0
   do n = 1, ns
      do m = 1, ns
        MM(n,m) = D(n,n) * F(n,m)
      end do
   end do

end subroutine Matxopt_ddff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_dfdf - (double) full matrix times diagonal, full matrix out
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_dfdf(F,D,ns,MM)

   integer(kind=4) :: n, m
   integer(kind=4), intent(in) :: ns
   real(kind=8), dimension(:,:), intent(in) :: D,F
   real(kind=8), dimension(:,:), intent(out) :: MM

   MM(:,:) = 0.0
   do n = 1, ns
      do m = 1, ns
        MM(n,m) = D(n,n) * F(n,m)
      end do
   end do

end subroutine Matxopt_dfdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_dffd - (double) full matrix time full matrix, diagonal out
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_dffd(F1,F2,ns,MM)

   integer(kind=4) :: n, m
   integer(kind=4), intent(in) :: ns
   real(kind=8), dimension(:,:), intent(in) :: F1,F2
   real(kind=8), dimension(:,:), intent(out) :: MM

   MM(:,:) = 0.0
   do n = 1, ns
      do m = 1, ns
        MM(n,n) = MM(n,n) + F1(m,n) * F2(n,m)
      end do
   end do
   !M = matmul(D,F)

end subroutine Matxopt_dffd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_zdff - (complex) diagonal times full matrix, full matrix out
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_zdff(D,F,ns,MM)

   integer(kind=4) :: n, m
   integer(kind=4), intent(in) :: ns
   complex(kind=8), dimension(:,:), intent(in) :: D,F
   complex(kind=8), dimension(:,:), intent(out) :: MM

   MM(:,:) = 0.0
   do n = 1, ns
      do m = 1, ns
        MM(n,m) = D(n,n) * F(n,m)
      end do
   end do

end subroutine Matxopt_zdff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_zfdf - (complex) full matrix times diagonal, full matrix out
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_zfdf(F,D,ns,MM)

   integer(kind=4) :: n, m
   integer(kind=4), intent(in) :: ns
   complex(kind=8), dimension(:,:), intent(in) :: D,F
   complex(kind=8), dimension(:,:), intent(out) :: MM

   MM(:,:) = 0.0
   do n = 1, ns
      do m = 1, ns
        MM(n,m) = D(n,n) * F(n,m)
      end do
   end do

end subroutine Matxopt_zfdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_zffd - (complex) full matrix time full matrix, diagonal out
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_zffd(F1,F2,ns,MM)

   integer(kind=4) :: n, m
   integer(kind=4), intent(in) :: ns
   complex(kind=8), dimension(:,:), intent(in) :: F1,F2
   complex(kind=8), dimension(:,:), intent(out) :: MM

   MM(:,:) = 0.0
   do n = 1, ns
      do m = 1, ns
        MM(n,n) = MM(n,n) + F1(m,n) * F2(n,m)
      end do
   end do

end subroutine Matxopt_zffd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_zffd2 - (complex) full matrix time full matrix, diagonal out (real)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_zffd2(F1,F2,ns,MM)

   integer(kind=4) :: n, m
   integer(kind=4), intent(in) :: ns
   complex(kind=8), dimension(:,:), intent(in) :: F1,F2
   real(kind=8), dimension(:,:), intent(out) :: MM

   MM(:,:) = 0.0
   do n = 1, ns
      do m = 1, ns
        MM(n,n) = MM(n,n) + real(F1(m,n) * F2(n,m))
      end do
   end do

end subroutine Matxopt_zffd2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_zdff - (complex) single value times full matrix, full matrix out
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_zsff(D,l,F,ns,MM)

   integer(kind=4) :: n, m
   integer(kind=4), intent(in) :: l, ns !l=equals location of single value
   complex(kind=8), dimension(:,:), intent(in) :: D, F
   complex(kind=8), dimension(:,:), intent(out) :: MM

   MM(:,:) = 0.0
   do n = 1, ns
     MM(n,l) = D(l,l) * F(n,l)
   end do

end subroutine Matxopt_zsff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_eoshift - (double) faster eoshift
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_eoshift(MM,ns,off,TT)

   integer(kind=4) :: n, m
   integer(kind=4), intent(in) :: ns, off !l=equals location of single value
   real(kind=8), dimension(:,:), intent(in) :: MM
   real(kind=8), dimension(:,:), intent(out) :: TT

   if (abs(off) .eq. 0) then
     TT(:,:) = MM
   else if(abs(off) .ge. ns-1) then
     TT(:,:) = 0.0
   else
     if(off .gt. 0) then !shift down
       TT(:,ns-(off+1):ns) = 0.0 !don't need to zero all, only values shift from outside domain
       !write(*,*) '1:ns-off,off+1:ns =',1-(ns-off),(off+1)-ns
       TT(:,1:ns-off) = MM(:,off+1:ns) !copy only values that will be moved
     else if(off .lt. 0) then !shift up
       TT(:,1:-off) = 0.0 !don't need to zero all, try to speed up
       !write(*,*) '-off:ns,ns+off =',(-off+1)-ns,1-(ns+off),ns
       TT(:,-off+1:ns) = MM(:,1:ns+off)
     end if
   end if

end subroutine Matxopt_eoshift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_eozero - (double) zero select range of matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_eozero(MM,ns,off)

   integer(kind=4) :: n, m
   integer(kind=4), intent(in) :: ns, off !negative off = zero below, positive off = zero above

   real(kind=8), dimension(:,:), intent(inout) :: MM

   if(abs(off) .ge. ns) then
     MM(:,:) = 0.0 !zero everything
   else if(off .lt. 0) then
     if(-off .gt. 0) MM(:,1:-off) = 0.0 !zero below off
   else
     MM(:,off:ns) = 0.0 !zero above off
   end if

end subroutine Matxopt_eozero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Matxopt_hilbert - Calculated the real part from the imag using Hilbert transformation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matxopt_hilbert(MM,E,dE,ns,TT)

   integer(kind=4) :: n, m
   integer(kind=4), intent(in) :: ns ! number of elements
   real(kind=8), intent(in) :: dE ! energy spacing
   real(kind=8), dimension(:), intent(in) :: E !energy range
   real(kind=8), dimension(:), intent(in) :: MM !imag input
   real(kind=8), dimension(:), intent(out) :: TT !real result

   TT(:) = 0.0
   do n = 1, ns
     do m = 1, ns
       if(E(m)-E(n) .ne. 0) TT(n) = TT(n) + 1/3.14159 * MM(m)/(E(m)-E(n))*dE !Hilbert transform
     end do
   end do

end subroutine Matxopt_hilbert

end module
