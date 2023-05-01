!******************************************************************************
!
! File : anderson_realmix.f90  -- (modified version of anderson_mixing.f90)
!   by : Alan Tackett
!   on : 07/19/99
!  for : density mixing (originally written for pwpaw code)
!
!  This module contains routines to implement the extended Anderson
!  mixing method as outlined in V. Eyert, J. Comp. Phys. 124,  271(1996).
!
!  The method is defined by eq's: 2.1, 8.2, 7.7
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

module Anderson_mixing

  implicit none

  type Anderson_Context  !** Anderson Mixing context
     real(8)    :: NewMix    !** Amount of new vectors to mix, ie beta in paper.
     integer :: Nmax      !** max number of vectors to keep 
     integer :: N         !** Current number of vectors in list
     integer :: Slot      !** New fill Slot
     integer :: VecSize   !** Size of each vector
     integer :: Err_Unit  !** Error unit
     real(8), pointer :: Matrix(:,:)
     real(8), pointer :: Gamma(:)  !** Gamma as defined in 7.6
     real(8), pointer :: DF(:,:)   !** Delta F
     real(8), pointer :: Fprev(:)
     real(8), pointer :: DX(:,:)
     real(8), pointer :: Xprev(:)

     ! temporary constants and arrays needed for each call to Anderson_Mix
     integer, pointer :: IPIV(:)
     real(8), pointer :: S(:)
     real(8), pointer :: RWork(:)
     real(8), pointer :: U(:,:)
     real(8), pointer :: VT(:,:)
     real(8), pointer :: Work(:)
     real(8), pointer :: DupMatrix(:,:)
     integer          :: Lwork
     integer          :: LRwork
     real(8)          :: ConditionNo
     real(8)          :: MachAccur
  end type Anderson_Context

contains

  !******************************************************************************
  !
  ! Anderson_Mix - Performs the actual mixing of the input vector with the
  !                history and retuns the result.
  !
  !   AC - Anderson context
  !   X  - Current vector on input and new guess on output
  !   F  - F(X) - X. Nonlinear mixing of input vector
  !
  ! Modified to call SVD routines
  !******************************************************************************

  subroutine Anderson_Mix(AC, X, F)
    type  (Anderson_Context), intent(inout) :: AC
    real(8), intent(inout) :: X(:)
    real(8), intent(in) :: F(:)

    integer :: i, slot, currentdim , n, j
    real(8) :: term
    real(8)  :: tmp

    !** First determine where to store the new correction vectors ***
    AC%slot = AC%slot + 1
    if (AC%Slot>AC%Nmax) AC%Slot = 1

    if ((AC%N < 0) .or. (AC%Nmax == 0)) then  !** Simple mixing for 1st time ***
       AC%Xprev = X
       X = X + AC%NewMix*F
    else
       slot = AC%Slot

       AC%DF(:,slot) = F - AC%Fprev   !** Make new DF vector
       AC%DX(:,slot) = X - AC%Xprev   !** Make new DX vector

       currentdim=min(AC%N+1,AC%Nmax)
       do i=1, currentdim              !*** Add row/col to matrix
          term = dot_product(AC%DF(:,i), AC%DF(:,slot))
          AC%Matrix(i,slot) = term
          if (i /= slot) AC%Matrix(slot,i) = (term)

          AC%Gamma(i) = dot_product(AC%DF(:,i), F)
       end do

       AC%DupMatrix = AC%Matrix

       n = AC%Nmax; j= currentdim
       call DGESDD('A',j,j,AC%DupMatrix(1,1),n,AC%S(1), &
            AC%U(1,1),n,AC%VT(1,1),n,AC%Work(1),AC%Lwork, AC%IPIV(1),i)
       if (i /= 0) then
          write(AC%Err_Unit,*) 'Anderson_Mix: Error in DGESDD. Error=',i
          tmp = 0
          tmp = 1.d0/tmp
          stop
       end if

       !WRITE(6,*) 'in Anderson_Mix -- completed SVD with values'
       !WRITE(6,'(1p5e15.7)') (AC%S(i),i=1,j)

       AC%Work(1:j) = AC%Gamma(1:j)
       AC%Gamma = 0
       tmp=max(abs(AC%S(1))/AC%ConditionNo,AC%Machaccur)
       do i=1,j
          if (abs(AC%S(i)).gt.tmp) then
             AC%Gamma(1:j)=AC%Gamma(1:j)+&
                  (AC%VT(i,1:j))*dot_product(AC%U(1:j,i),AC%Work(1:j))/AC%S(i)
          end if
       enddo

       AC%Xprev = X

       !*** Now calculate the new vector ***
       X = X + AC%NewMix*F

       !do i=1, min(AC%N+1,AC%Nmax)     !*** Add row/col to matrix
       do i=1, currentdim               ! updated vector
          X = X - AC%Gamma(i)*(AC%DX(:,i) + AC%NewMix*AC%DF(:,i))
       end do
    end if

    AC%Fprev = F

    AC%N = AC%N + 1
    if (AC%N > AC%Nmax) AC%N = AC%Nmax

    return
  end subroutine Anderson_Mix

  !******************************************************************************
  !
  !  Anderson_ResetMix - Resets the mixing history to None
  !
  !     AC - Anderson context to reset
  !
  !******************************************************************************

  subroutine Anderson_ResetMix(AC)
    type  (Anderson_Context), intent(inout) :: AC

    AC%N = -1
    AC%Slot = -1

    return
  end subroutine Anderson_ResetMix

  !******************************************************************************
  !
  !  FreeAnderson - Frees all the data associated with the AC data structure
  !    
  !      AC -pointer to the Anderson context to free
  !
  !******************************************************************************

  subroutine FreeAnderson(AC)
    type (Anderson_Context), pointer :: AC  

    deallocate(AC%Xprev, AC%Fprev , AC%DX, AC%DF, AC%Matrix, AC%Gamma)
    deallocate(AC%DupMatrix, AC%U, AC%VT, AC%Work, AC%RWork, AC%IPIV, AC%S)
    deallocate(AC)

    return
  end subroutine FreeAnderson

  !******************************************************************************
  !
  !  initAnderson - initializes and Anderson_Context data structure for use
  !
  !   AC       - Anderson context created and returned
  !   Nmax     - max number of vectors to keep
  !   VecSize  - Size of each vector
  !   NewMix   - Mixing factor
  !
  !******************************************************************************

  subroutine initAnderson(AC, Nmax, VecSize, NewMix, CondNo)
    type (Anderson_Context), pointer :: AC
    integer(kind=4), intent(in) :: Nmax
    integer(kind=4), intent(in) :: VecSize
    real(kind=8), intent(in) :: NewMix
    real(kind=8), intent(in) :: CondNo

    integer :: i
    real(kind=8) :: a1,a2,a3

    allocate(AC) !** allocate the pointer

    AC%Nmax = Nmax !*** Store the contants
    AC%VecSize = VecSize
    AC%NewMix = NewMix
    AC%Err_Unit = 6

    AC%N = -1 !** init the rest of the structure
    AC%Slot = -1

    allocate(AC%Xprev(VecSize), AC%Fprev(VecSize) , AC%DX(VecSize,Nmax), &
         AC%DF(VecSize,Nmax), AC%Matrix(Nmax,Nmax) , AC%Gamma(Nmax), &
         Stat=i)

    AC%Lwork=5*Nmax*Nmax+10*Nmax
    AC%LRwork= 5*Nmax*Nmax+7*Nmax
    AC%ConditionNo= CondNo

    ! Calculate machine accuracy
    AC%Machaccur = 0
    a1 = 4.d0/3.d0
    do while (AC%Machaccur == 0.d0)
       a2 = a1 - 1.d0
       a3 = a2 + a2 + a2
       AC%Machaccur = abs(a3 - 1.d0)
    enddo

    allocate(AC%DupMatrix(Nmax,Nmax), AC%U(Nmax, Nmax), AC%VT(Nmax,Nmax), &
         AC%Work(AC%Lwork), AC%RWork(AC%LRWork), AC%IPIV(8*Nmax), &
         AC%S(Nmax), STAT=i)

    AC%Matrix = 0

    return
  end subroutine initAnderson


end module anderson_mixing
