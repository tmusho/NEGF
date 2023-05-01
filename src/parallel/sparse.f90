!  sparse.f90 
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
!  PROGRAM: Sparse(Module)
!
!  PURPOSE: Contains subroutine to sparse LU decomp matrix
!
!  SUBROUTINES: 
!               
!               
!                
!
!****************************************************************************
!
module Sparse

!Use global variables from module
use Startup, ONLY: eye, rank
use Errors
use Debug

implicit none

   integer :: ldb, nrhs, nds, symbolic, numeric
   integer, dimension(64) :: iparm
   integer, dimension(20) :: control
   integer, dimension(90) :: umf_info
   integer, allocatable, dimension(:) :: ia, ja
   real, dimension(64) :: dparm
   complex(kind=8), allocatable, dimension(:) :: avals, b
   complex(kind=8) :: rmisc

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Sparse_superLU - Sparse factorization and back solver for inv(G), 10x faster
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Sparse_superLU(n,ca,M)

   external c_fortran_zgssv

   integer :: ierr, iopt, info, i
   integer(kind=8) :: factors
   character*120 :: var
   integer, intent(in) :: n
   integer, intent(in) :: ca
   complex(kind=8), dimension(:,:), intent(inout) :: M

   call dnscsr(n, n, nds, M, n, avals, ja, ia, ierr) !put matrix in CRC format

   !Calculate LU factors using sparse solver
   iopt = 1
   call c_fortran_zgssv( iopt, n, nds, nrhs, avals, ja, ia, b, ldb, factors, info )
   !write(*,*) 'zgssv =',info

   !Calculate inverse col by col
   iopt = 2
   do i=1, n
     b = eye(i,:)
     call c_fortran_zgssv( iopt, n, nds, nrhs, avals, ja, ia, b, ldb, factors, info )
     M(i,:) = b
   end do

   !Free LU storage
   ! if you don't deallocate there is a memory gathering and ram runs out.
   iopt = 3
   call c_fortran_zgssv( iopt, n, nds, nrhs, avals, ja, ia, b, ldb, factors, info )
   !call csrdns(n, n, avals, ja, ia, M, n, ierr)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Sparse_superLU_init - Opensource LU solver initallization
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Sparse_superLU_init(n,M)

   integer :: ierr
   integer, intent(in) :: n
   complex(kind=8), dimension(:,:), intent(inout) :: M

   iparm(:) = 0; dparm(:) = 0 !zero entries
   iparm(33) = 1 ! 1 cpu only, stupid closed source

   call dnnze(n,n,M,nds,ierr)
   if(rank .eq. 0) then
     if(debug_l>5)write(*,*) 'Number of dense matrix entries: ',n*n
     if(debug_l>5)write(*,*) 'Number of non-zero entries: ',nds
   end if
   ldb=n; nrhs=1

   allocate(ia(n+1),ja(nds)) !ja .eq. number of non zeros FIXME
   allocate(avals(nds),b(n))

   b(:) = 1.0

   call dnscsr(n, n, nds, M, n, avals, ja, ia, ierr) !put matrix in CRC format

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Sparse_zgsmp - WASP sparse solver (segfault)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Sparse_zgsmp(n,ca,M)

   integer :: ierr
   integer, intent(in) :: n
   integer, intent(in) :: ca
   complex(kind=8), dimension(:,:), intent(inout) :: M

   call dnscsr(n, n, nds, M, n, avals, ja, ia, ierr) !put matrix in CRC format
   write(*,*) 'put in crc format',ierr

   !Don't need to analyze each time if same permutations
   if(ca .eq. 0) then
     ! Analysis.
     iparm(2) = 1; iparm(3) = 1
     !call zgsmp(n, ia, ja, avals, b, ldb, nrhs, rmisc, iparm, dparm)
     if(iparm(64) .ne. 0) stop
     write(*,*) 'analysis done',iparm(64)
   end if

   ! Factorization.
   iparm(2) = 2; iparm(3) = 2
   !call zgsmp(n, ia, ja, avals, b, ldb, nrhs, rmisc, iparm, dparm)
   if(iparm(64) .ne. 0) stop
   write(*,*) 'Factorization done',iparm(64)

   ! Backsubstitution
   iparm(2) = 3; iparm(3) = 3
   !call zgsmp (n, ia, ja, avals, b, ldb, nrhs, rmisc, iparm, dparm)
   if(iparm(64) .ne. 0) stop
   write(*,*) 'backsub done',iparm(64),avals

   call csrdns(n, n, avals, ja, ia, M, n, ierr)

end subroutine Sparse_zgsmp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Sparse_zgsmp
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine Sparse_umfpack(n,ca,M)
!
!   integer :: ierr
!   integer, intent(in) :: n
!   integer, intent(in) :: ca
!   complex(kind=8), dimension(:,:), intent(inout) :: M
!
!   call dnscsr(n, n, nds, M, n, avals, ja, ia, ierr) !put matrix in CRC format
!   write(*,*) 'put in crc format',ierr
!
!   !Don't need to analyze each time if same permutations
!   if(ca .eq. 0) then
!     ! Analysis.
!     iparm(2) = 1; iparm(3) = 1
!     !call zgsmp(n, ia, ja, avals, b, ldb, nrhs, rmisc, iparm, dparm)
!     call umf4zsym(n, n, ia, ja, avals, symbolic, control, umf_info)
!     write(*,*) 'analysis done',iparm(64)
!   end if
!
!   ! Factorization.
!   iparm(2) = 2; iparm(3) = 2
!   !call zgsmp(n, ia, ja, avals, b, ldb, nrhs, rmisc, iparm, dparm)
!   call umf4znum(ia, ja, avals, symbolic, control, umf_info)
!   write(*,*) 'Factorization done',iparm(64)!
!
!   ! Backsubstitution
!   iparm(2) = 3; iparm(3) = 3
!   !call zgsmp (n, ia, ja, avals, b, ldb, nrhs, rmisc, iparm, dparm)
!   call umf4zsol(sys, x, b, numeric, control, umf_info)
 !  write(*,*) 'backsub done',iparm(64),avals
!   call csrdns(n, n, avals, ja, ia, M, n, ierr)
!
!end subroutine Sparse_umfpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Sparse_zgsmp_init - Sparse solver, proprietary code
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Sparse_zgsmp_init(n,M)

   ! This solver is suppose to be faster then superLU but isn't open source
   ! I was getting a segfault in the library so I gave up.
   ! http://www.alphaworks.ibm.com/tech/wsmp
   
   integer :: ierr
   integer, intent(in) :: n
   complex(kind=8), dimension(:,:), intent(inout) :: M

   iparm(:) = 0; dparm(:) = 0 !zero entries
   iparm(33) = 1 ! 1 cpu only, stupid closed source
   !call wsetmaxthrds(1)

   call dnnze(n,n,M,nds,ierr)
   write(*,*) 'Number of dense matrix entries: ',n*n
   write(*,*) 'Number of non-zero entries: ',nds
   ldb=n; nrhs=1

   allocate(ia(n+1),ja(nds)) !ja .eq. number of non zeros FIXME
   allocate(avals(nds),b(n))

   b(:) = 1.0

   call dnscsr(n, n, nds, M, n, avals, ja, ia, ierr) !put matrix in CRC format

   iparm(1) = 0; iparm(2) = 0; iparm(3) = 0
   !call zgsmp (n, ia, ja, avals, b, ldb, nrhs, rmisc, iparm, dparm)
   write(*,*) 'initallize zgsmp',iparm(64),iparm(33),iparm(24)
   if(iparm(64) .ne. 0) stop

end subroutine Sparse_zgsmp_init

subroutine dnnze(nrow, ncol, dns, nds, ierr)
!*****************************************************************************80
!
!! dnnnze Computes number of non-zero entries.
!
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NCOL, the column dimension of the matrix.
!
!    Input, real DNS(NDNS,NCOL), an NROW by NCOL dense matrix.
!
  integer i, j
  integer ncol
  integer ndns
  integer nrow
  integer nds
  integer ierr
  complex(kind = 8), intent(in) :: dns(nrow,ncol)
  nds = 0

  do i = 1, nrow
    do j = 1, ncol
      if ( dns(i,j) /= (0.0D+00,0.0D+00) ) then
        nds = nds + 1
      end if
    end do
  end do

end subroutine dnnze

subroutine dnscsr(nrow, ncol, nzmax, dns, ndns, a, ja, ia, ierr)

!*****************************************************************************80
!
!! DNSCSR converts Dense to Compressed Row Sparse format.
!
!  Discussion:
!
!    This routine converts a densely stored matrix into a row orientied
!    compactly sparse matrix.  It is the reverse of CSRDNS.
!
!    This routine does not check whether an element is small.  It considers 
!    that A(I,J) is zero only if it is exactly equal to zero.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NCOL, the column dimension of the matrix.
!
!    Input, integer NZMAX, the maximum number of nonzero elements 
!    allowed.  This should be set to be the lengths of the arrays A and JA.
!
!    Input, real DNS(NDNS,NCOL), an NROW by NCOL dense matrix.
!
!    Input, integer NDNS, the first dimension of DNS, which must be
!    at least NROW.
!
!    Output, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, integer IERR, error indicator.
!    0 means normal return;
!    I, means that the the code stopped while processing row I, because
!       there was no space left in A and JA, as defined by NZMAX.
!
  implicit none

  integer ncol
  integer ndns
  integer nrow

  complex (kind = 8), intent(out) :: a(*)
  complex (kind = 8), intent(in) :: dns(ndns,ncol)
  integer i
  integer ia(nrow+1)
  integer ierr
  integer j
  integer ja(*)
  integer next
  integer nzmax

  ierr = 0
  next = 1
  ia(1) = 1

  do i = 1, nrow

    do j = 1, ncol

      if ( dns(i,j) /= (0.0D+00,0.0D+00) ) then

        if ( nzmax < next ) then
          ierr = i
          return
        end if

        ja(next) = j
        a(next) = dns(i,j)
        next = next + 1

      end if

    end do

    ia(i+1) = next

  end do

end subroutine dnscsr

subroutine csrdns(nrow, ncol, a, ja, ia, dns, ndns, ierr)

!*****************************************************************************80
!
!! CSRDNS converts Compressed Sparse Row to Dense format.
!
!  Discussion:
!
!    This routine converts a row-stored sparse matrix into a densely stored one.
!
!  Modified:
!
!    07 January 2004
!
!  Author:
!
!    Youcef Saad
!
!  Parameters:
!
!    Input, integer NROW, the row dimension of the matrix.
!
!    Input, integer NCOL, the column dimension of the matrix.
!
!    Input, real A(*), integer JA(*), IA(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
!    Output, real DNS(NDNS,NDNS), the dense array containing a
!    copy of the matrix.
!
!    Input, integer NDNS, the dimension of the DNS array.
!
!    Output, integer IERR, error indicator.
!    0, means normal return
!    i, means that the code has stopped when processing
!       row number i, because it found a column number > ncol.
!
  implicit none

  integer ncol
  integer ndns

  complex(kind = 8),intent(in) :: a(*)
  complex(kind = 8),intent(out) :: dns(ndns,ncol)
  integer i
  integer ia(*)
  integer ierr
  integer j
  integer ja(*)
  integer k
  integer nrow
  
  ierr = 0
  dns(1:nrow,1:ncol) = (0.0D+00,0.0D+00)

  do i = 1, nrow
    do k = ia(i), ia(i+1)-1
      j = ja(k)
      if ( ncol < j ) then
        ierr = i
        return
      end if
      dns(i,j) = a(k)
    end do
  end do

  return
end subroutine

end module
