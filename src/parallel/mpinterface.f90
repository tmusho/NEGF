!  mpinterface.f90 
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
!  PROGRAM: mpinterface(Module)
!
!  PURPOSE: Contains message passing file interface
!
!  SUBROUTINES: 
!               
!                
!
!****************************************************************************
!
module Mpinterface

use Startup, ONLY : pi,q

implicit none

integer(kind=4), parameter, private :: WAITTIME = 1
integer(kind=4), parameter, private :: MAXLOOP = 1000000000
integer(kind=4), parameter, private :: UnitCouple = 12
integer(kind=4), private :: tr = 1
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Mpinterface_Write_Header - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mpinterface_Write_Header(Lb,Lw,Lc,Nl,mat1,mat2,Lg,Lg_o,Tp,dT,debug_l)

   integer(kind=4), intent(inout) :: debug_l
   integer(kind=4), intent(in) :: mat1,mat2,Nl,Lg_o
   real(kind=8), intent(in) :: Lb,Lw,Lc,Tp,Lg,dT

   write(6,*)
   write(6,*)'+++ Coupled Run; Sending Header Data to NEGF Phonon +++'

   open(UnitCouple,file='elements.bin',form='unformatted')
   write(UnitCouple)Lb,Lw,Lc,Nl,mat1,mat2,Lg,Lg_o,Tp,dT,debug_l
   close(UnitCouple)

   open(UnitCouple,file='nodeok.mpf')
   write(UnitCouple,*)'echo nodes written'
   close(UnitCouple)

   if(debug_l>5)write(6,*)' NEGF Electron: HEADER MESSAGE written'

end subroutine Mpinterface_Write_Header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Mpinterface_Read_Header - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mpinterface_Read_Header(Lb,Lw,Lc,Nl,mat1,mat2,Lg,Lg_o,Tp,dT,debug_l)

   integer(kind=4), intent(inout) :: debug_l
   integer(kind=4), intent(out) :: mat1,mat2,Nl,Lg_o
   real(kind=8), intent(out) :: Lb,Lw,Lc,Tp,Lg,dT

   write(6,*)
   write(6,*)'+++ Coupled Run; Waiting for NEGF Electron Data Header File+++'

   call WaitOnFile('nodeok.mpf','ph_stop.mpf')

   open(UnitCouple,file='elements.bin',form='unformatted')
   read(UnitCouple)Lb,Lw,Lc,Nl,mat1,mat2,Lg,Lg_o,Tp,dT,debug_l
   close(UnitCouple,status='delete')

   close(UnitCouple)
   write(6,*)' Header File Read'

   call SafeOpen(UnitCouple,'nodeok.mpf','lag_stop.mpf')
   close(UnitCouple,status='delete')

   if(debug_l>5)write(6,*)' NEGF Phonon: HEADER MESSAGE read'

end subroutine Mpinterface_Read_Header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Mpinterface_Write_Eig - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mpinterface_Write_Eig(Eig,Nb,debug_l)

   integer(kind=4) :: r
   integer(kind=4), intent(inout) :: debug_l
   integer(kind=4), intent(in) :: Nb
   real(kind=8), intent(in), dimension(:) :: Eig

   if(debug_l>5)write(6,*)'+++ Coupled Run; Writing Eig for NEGF Electron Data +++'

   open(UnitCouple,file='lag_out.bin',form='unformatted')

   write(UnitCouple) Nb
   do r = 1, Nb
      write(UnitCouple) Eig(r)
   end do

   close(UnitCouple)

   call SafeOpen(UnitCouple,'ee_ok.mpf','ph_stop.mpf') !stop due to phonon error
   close(UnitCouple)

   call SafeOpen(UnitCouple,'ph_ok.mpf','ee_stop.mpf')
   write(UnitCouple,*)'ph_out.mpf done'
   close(UnitCouple)

end subroutine Mpinterface_Write_Eig

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Mpinterface_Read_Eig - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mpinterface_Read_Eig(Eig,debug_l)

   integer(kind=4) :: r, isopen, err
   integer(kind=4) :: Nb
   integer(kind=4), intent(in) :: debug_l
   real(kind=8), allocatable, intent(out), dimension(:) :: Eig

   if(debug_l>5)write(6,*)'+++ Coupled Run; Waiting for Eig Data from NEGF Phonon +++'
   call WaitOnFile('ee_ok.mpf','ph_stop.mpf')
   call RemoveWaitFile(UnitCouple,'ee_ok.mpf','ph_stop.mpf')

   open(UnitCouple,FILE='lag_out.bin',form='unformatted',IOSTAT=isopen)

   read(UnitCouple) Nb
   allocate(Eig(Nb), stat=err) ! Input - Diagonal Values, Output - Allocation Eigenvalue Vector
   Eig(:) = 0.0
   do r = 1, Nb
      read(UnitCouple) Eig(r)
   end do

   close(UnitCouple)

   call SafeOpen(UnitCouple,'ph_ok.mpf','ph_stop.mpf')
   close(UnitCouple,status='delete')
   if(debug_l>5)write(6,*)' NEGF Electron: Phonon Eigenvalues read'

end subroutine Mpinterface_Read_Eig

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Mpinterface_Write_Gn - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mpinterface_Write_Gn(GG,NE,Eo,Ef,Np,Nb,debug_l,chng_tprt)

   integer(kind=4) :: r, j, m
   integer(kind=4), intent(in) :: debug_l
   integer(kind=4), intent(in) :: NE, Np, Nb
   real(kind=8), intent(in) :: chng_tprt
   real(kind=8), intent(in) :: Eo, Ef
   real(kind=8), intent(in), dimension(:,:,:) :: GG

   real(kind=8) :: Eo_t, Ef_t

   if(debug_l>5)write(6,*)'+++ Coupled Run; Writing Gn for NEGF Electron Data +++'
   call WaitRemoveFile('ph_ok.mpf','ee_stop.mpf')

   open(UnitCouple,file='lag_out.bin',form='unformatted')

   Eo_t = Eo/q; Ef_t = Ef/q !convert to eV

   write(UnitCouple) NE,Np,Nb,Eo_t,Ef_t,chng_tprt
   do r = 1, NE
     do j = 1, Np
       do m = 1, Nb
         write(UnitCouple) GG(j,r,m)
       end do
     end do
   end do

   close(UnitCouple)

   call SafeOpen(UnitCouple,'ee_ok.mpf','ph_stop.mpf')
   close(UnitCouple)

   call SafeOpen(UnitCouple,'ph_ok.mpf','ee_stop.mpf')
   write(UnitCouple,*)'ph_out.mpf done'
   close(UnitCouple)

end subroutine Mpinterface_Write_Gn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Mpinterface_Read_Gn - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mpinterface_Read_Gn(GG,NE,Eo,Ef,Np,Nb,debug_l,chng_tprt)

   integer(kind=4) :: r, j, m, isopen, err
   integer(kind=4), intent(out) :: NE, Np, Nb
   integer(kind=4), intent(in) :: debug_l
   real(kind=8), intent(out) :: Eo, Ef,chng_tprt
   real(kind=8), allocatable, intent(out), dimension(:,:,:) :: GG

   if(debug_l>5)write(6,*)'+++ Coupled Run; Waiting for Gn Data from NEGF Phonon +++'
   call WaitOnFile('ee_ok.mpf','ph_stop.mpf')
   call RemoveWaitFile(UnitCouple,'ee_ok.mpf','ph_stop.mpf')

   open(UnitCouple,FILE='lag_out.bin',form='unformatted',IOSTAT=isopen)

   read(UnitCouple) NE,Np,Nb,Eo,Ef,chng_tprt
   allocate(GG(NE,Np,Nb), stat=err) ! Input - Diagonal Values, Output - Allocation Eigenvalue Vector
   GG(:,:,:) = 0.0
   do r = 1, NE
     do j = 1, Np
       do m = 1, Nb
         read(UnitCouple) GG(r,j,m)
       end do
     end do
   end do

   close(UnitCouple)

   call SafeOpen(UnitCouple,'ph_ok.mpf','lagstop.mpf')
   close(UnitCouple,status='delete')
   if(debug_l>5)write(6,*)' NEGF Electron: Phonon Gn read'

end subroutine Mpinterface_Read_Gn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Mpinterface_Read_Gn - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mpinterface_Read_Gn2(GG,NE,Eo,Ef,Np,Nb,debug_l,chng_tprt)

   integer(kind=4) :: r, j, m, isopen, err
   integer(kind=4), intent(out) :: NE, Np, Nb
   integer(kind=4), intent(in) :: debug_l
   real(kind=8), intent(out) :: Eo, Ef,chng_tprt
   real(kind=8), intent(out), dimension(:,:,:) :: GG

   if(debug_l>5)write(6,*)'+++ Coupled Run; Waiting for Gn Data from NEGF Phonon +++'
   call WaitOnFile('ee_ok.mpf','ph_stop.mpf')
   call RemoveWaitFile(UnitCouple,'ee_ok.mpf','ph_stop.mpf')

   open(UnitCouple,FILE='lag_out.bin',form='unformatted',IOSTAT=isopen)

   read(UnitCouple) NE,Np,Nb,Eo,Ef,chng_tprt
   GG(:,:,:) = 0.0
   do r = 1, NE
     do j = 1, Np
       do m = 1, Nb
         read(UnitCouple) GG(r,j,m)
       end do
     end do
   end do

   close(UnitCouple)

   call SafeOpen(UnitCouple,'ph_ok.mpf','lagstop.mpf')
   close(UnitCouple,status='delete')
   if(debug_l>5)write(6,*)' NEGF Electron: Phonon Gn read'

end subroutine Mpinterface_Read_Gn2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Mpinterface_Write_Gp - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mpinterface_Write_Gp(GG,NE,Np,Nb,debug_l,chng_tprt)

   integer(kind=4) :: r, j, m
   integer(kind=4), intent(in) :: debug_l
   integer(kind=4), intent(in) :: NE, Np, Nb
   real(kind=8), intent(in) :: chng_tprt
   real(kind=8), intent(in), dimension(:,:,:) :: GG

   if(debug_l>5)write(6,*)'+++ Coupled Run; Writing Gp for NEGF Electron Data +++'
   call WaitRemoveFile('ph_ok.mpf','ee_stop.mpf')

   open(UnitCouple,file='lag_out.bin',form='unformatted')

   write(UnitCouple) NE,Np,Nb,chng_tprt
   do r = 1, NE
     do j = 1, Np
       do m = 1, Nb
         write(UnitCouple) GG(j,r,m)
       end do
     end do
   end do

   close(UnitCouple)

   call SafeOpen(UnitCouple,'ee_ok.mpf','ph_stop.mpf')
   close(UnitCouple)

   call SafeOpen(UnitCouple,'ph_ok.mpf','ee_stop.mpf')
   write(UnitCouple,*)'lag_out.mpf done'
   close(UnitCouple)

end subroutine Mpinterface_Write_Gp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Mpinterface_Read_Gp - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mpinterface_Read_Gp(GG,NE,Np,Nb,debug_l,chng_tprt)

   integer(kind=4) :: r, j, m, isopen, err
   integer(kind=4), intent(out) :: NE, Np, Nb
   integer(kind=4), intent(in) :: debug_l
   real(kind=8), intent(out) :: chng_tprt
   real(kind=8), allocatable, intent(out), dimension(:,:,:) :: GG

   if(debug_l>5)write(6,*)'+++ Coupled Run; Waiting for Gp Data from NEGF Phonon +++'
   call WaitOnFile('ee_ok.mpf','ph_stop.mpf')
   call RemoveWaitFile(UnitCouple,'ee_ok.mpf','ph_stop.mpf')

   open(UnitCouple,FILE='lag_out.bin',form='unformatted',IOSTAT=isopen)

   read(UnitCouple) NE,Np,Nb,chng_tprt
   allocate(GG(NE,Np,Nb), stat=err) ! Input - Diagonal Values, Output - Allocation Eigenvalue Vector
   GG(:,:,:) = 0.0
   do r = 1, NE
     do j = 1, Np
       do m = 1, Nb
         read(UnitCouple) GG(r,j,m)
       end do
     end do
   end do

   close(UnitCouple)

   call SafeOpen(UnitCouple,'ph_ok.mpf','ph_stop.mpf')
   close(UnitCouple,status='delete')
   if(debug_l>5)write(6,*)' NEGF Electron: Phonon Gn read'

end subroutine Mpinterface_Read_Gp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Mpinterface_Read_Gp - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mpinterface_Read_Gp2(GG,NE,Np,Nb,debug_l,chng_tprt)

   integer(kind=4) :: r, j, m, isopen, err
   integer(kind=4), intent(out) :: NE, Np, Nb
   integer(kind=4), intent(in) :: debug_l
   real(kind=8), intent(out) :: chng_tprt
   real(kind=8), intent(out), dimension(:,:,:) :: GG

   if(debug_l>5)write(6,*)'+++ Coupled Run; Waiting for Gp Data from NEGF Phonon +++'
   call WaitOnFile('ee_ok.mpf','ph_stop.mpf')
   call RemoveWaitFile(UnitCouple,'ee_ok.mpf','ph_stop.mpf')

   open(UnitCouple,FILE='lag_out.bin',form='unformatted',IOSTAT=isopen)

   read(UnitCouple) NE,Np,Nb,chng_tprt
   GG(:,:,:) = 0.0
   do r = 1, NE
     do j = 1, Np
       do m = 1, Nb
         read(UnitCouple) GG(r,j,m)
       end do
     end do
   end do

   close(UnitCouple)

   call SafeOpen(UnitCouple,'ph_ok.mpf','ph_stop.mpf')
   close(UnitCouple,status='delete')
   if(debug_l>5)write(6,*)' NEGF Electron: Phonon Gn read'

end subroutine Mpinterface_Read_Gp2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Mpinterface_Write_Fermi - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mpinterface_Write_Fermi(f_ph,Nb,debug_l)

   integer(kind=4) :: r, j, m
   integer(kind=4), intent(in) :: debug_l
   integer(kind=4), intent(in) :: Nb
   real(kind=8), intent(in), dimension(:) :: f_ph

   if(debug_l>5)write(6,*)'+++ Coupled Run; Writing Fermi Data for NEGF Electron Data +++'
   call WaitRemoveFile('ph_ok.mpf','ee_stop.mpf')

   open(UnitCouple,file='lag_out.bin',form='unformatted')

   write(UnitCouple) Nb
   do r = 1, Nb
     write(UnitCouple) f_ph(r)
   end do

   close(UnitCouple)
   call SafeOpen(UnitCouple,'ee_ok.mpf','ph_stop.mpf')

   close(UnitCouple)

   call SafeOpen(UnitCouple,'ph_ok.mpf','ee_stop.mpf')
   write(UnitCouple,*)'lag_out.mpf done'
   close(UnitCouple)

end subroutine Mpinterface_Write_Fermi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Mpinterface_Read_Fermi - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mpinterface_Read_Fermi(f_ph,Nb,debug_l)

   integer(kind=4) :: r, j, m, isopen, err
   integer(kind=4), intent(out) :: Nb
   integer(kind=4), intent(in) :: debug_l
   real(kind=8), allocatable, intent(out), dimension(:) :: f_ph

   if(debug_l>5)write(6,*)'+++ Coupled Run; Waiting for Fermi Data from NEGF Phonon +++'
   call WaitOnFile('ee_ok.mpf','ph_stop.mpf')
   call RemoveWaitFile(UnitCouple,'ee_ok.mpf','ph_stop.mpf')

   open(UnitCouple,FILE='lag_out.bin',form='unformatted',IOSTAT=isopen)

   read(UnitCouple) Nb
   allocate(f_ph(Nb), stat=err) ! Input - Diagonal Values, Output - Allocation Eigenvalue Vector
   f_ph(:) = 0.0
   do r = 1, Nb
     read(UnitCouple) f_ph(r)
   end do

   close(UnitCouple)

   call SafeOpen(UnitCouple,'ph_ok.mpf','ph_stop.mpf')
   close(UnitCouple,status='delete')
   if(debug_l>5)write(6,*)' NEGF Electron: Phonon Fermi read'

end subroutine Mpinterface_Read_Fermi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Mpinterface_Read_Fermi - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mpinterface_Read_Fermi2(f_ph,Nb,debug_l)

   integer(kind=4) :: r, j, m, isopen, err
   integer(kind=4), intent(out) :: Nb
   integer(kind=4), intent(in) :: debug_l
   real(kind=8), intent(out), dimension(:) :: f_ph

   if(debug_l>5)write(6,*)'+++ Coupled Run; Waiting for Fermi Data from NEGF Phonon +++'
   call WaitOnFile('ee_ok.mpf','ph_stop.mpf')
   call RemoveWaitFile(UnitCouple,'ee_ok.mpf','ph_stop.mpf')

   open(UnitCouple,FILE='lag_out.bin',form='unformatted',IOSTAT=isopen)

   read(UnitCouple) Nb
   f_ph(:) = 0.0
   do r = 1, Nb
     read(UnitCouple) f_ph(r)
   end do

   close(UnitCouple)

   call SafeOpen(UnitCouple,'ph_ok.mpf','ph_stop.mpf')
   close(UnitCouple,status='delete')
   if(debug_l>5)write(6,*)' NEGF Electron: Phonon Fermi read'

end subroutine Mpinterface_Read_Fermi2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Mpinterface_Write_Scph - Write scattered phonons
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mpinterface_Write_Scph(Gph,Np_p,Nb_p,sc_ee,debug_l)

   integer(kind=4) :: r, j, m, sc
   integer(kind=4), intent(in) :: debug_l, sc_ee
   integer(kind=4), intent(in) :: Np_p, Nb_p
   real(kind=8), intent(in), dimension(:,:,:) :: Gph

   if(debug_l>5)write(6,*)'+++ Coupled Run; Writing Scattering Data +++'
   call WaitRemoveFile('ph_ok.mpf','ee_stop.mpf')

   open(UnitCouple,file='lag_out.bin',form='unformatted')

   write(UnitCouple) sc_ee !phonon break loop
   write(UnitCouple) Np_p,Nb_p
   do j = 1, Np_p
     do m = 1, Nb_p
       write(UnitCouple) Gph(j,1,m)
     end do
   end do

   close(UnitCouple)

   call SafeOpen(UnitCouple,'ee_ok.mpf','ph_stop.mpf')
   close(UnitCouple)

   call SafeOpen(UnitCouple,'cl_ok.mpf','ee_stop.mpf')
   write(UnitCouple,*)'cl_out.mpf done'
   close(UnitCouple)

   call WaitRemoveFile('cl_ok.mpf','ee_stop.mpf')

end subroutine Mpinterface_Write_Scph

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Mpinterface_Read_Scph - Read Scattered Phonons
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mpinterface_Read_Scph(Gph,Np_p,Nb_p,sc_ee,debug_l)

   integer(kind=4) :: r, j, m, isopen, err
   integer(kind=4), intent(out) :: Np_p, Nb_p, sc_ee
   integer(kind=4), intent(in) :: debug_l
   real(kind=8), intent(inout), dimension(:,:) :: Gph

   if(debug_l>5)write(6,*)'+++ Coupled Run; Waiting Scattered Data from NEGF Phonon +++'
   call WaitOnFile('cl_ok.mpf','ph_stop.mpf')
   call RemoveWaitFile(UnitCouple,'cl_ok.mpf','ph_stop.mpf')
   call WaitOnFile('ee_ok.mpf','ee_stop.mpf')
   call RemoveWaitFile(UnitCouple,'ee_ok.mpf','ph_stop.mpf')

   open(UnitCouple,FILE='lag_out.bin',form='unformatted',IOSTAT=isopen)

   read(UnitCouple) sc_ee
   if(sc_ee .gt. 0) then
   read(UnitCouple) Np_p,Nb_p
   do j = 1, Np_p
     do m = 1, Nb_p
       read(UnitCouple) Gph(j,m)
     end do
   end do
   end if

   close(UnitCouple)

   call SafeOpen(UnitCouple,'ph_ok.mpf','ph_stop.mpf')
   close(UnitCouple,status='delete')

   if(debug_l>5)write(6,*)'+++ Coupled Run; Read Scattered Data from NEGF Phonon +++'

end subroutine Mpinterface_Read_Scph


subroutine Mpinterface_Clean
   call RemoveWaitFile(UnitCouple,'ee_ok.mpf','ph_stop.mpf')
   call RemoveWaitFile(UnitCouple,'ph_ok.mpf','ph_stop.mpf')
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Mpinterface_Wait - 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WaitOnFile(filename,stopfile)
!  Wait on File to appear

!   ----- arguments -----
   character(*), intent(in) :: filename    ! name of file to be opened
   character(*), intent(in) :: stopfile    ! FE code stopped file

!   ----- local variable -----
   logical(4) :: here,done
   integer(4) :: ict,iend
   character(len=80) :: line

   inquire(file=filename,exist=here)

   ict=0
   done=.false.
   do while((.NOT.here).AND.(ict.le.MAXLOOP).AND.(.NOT.done))
      ict=ict+1
      inquire(file=filename,exist=here)
        call sleep(waittime)   !unix comment
      inquire(file=stopfile,exist=done)
   enddo

   if(.not.here)then
      call CouplerIO_Abort('CouplerIO WaitOnFile :: timed out waiting for file '//filename)
   else if(done)then
      call CouplerIO_Abort('CouplerIO WaitOnFile ::  Stop File Found; file name is '//stopfile)
   endif

   return

end subroutine WaitOnFile

subroutine RemoveWaitFile(filenum,filename,stopfile)
!recover for file open error and try again

!   ----- arguments -----
   integer(4), intent(in) :: filenum        ! unit number of file to be opened
   character(*), intent(in) :: filename    ! name of file to be opened
   character(*), intent(in) :: stopfile    ! FE code stopped file

   integer(kind=4) :: isopen, err

   open(filenum,file=filename,form='unformatted',IOSTAT=isopen)
   close(filenum,status='delete')

end subroutine RemoveWaitFile

subroutine SafeOpen(filenum,filename,stopfile)
!recover for file open error and try again

!   ----- arguments -----
   integer(4), intent(in) :: filenum        ! unit number of file to be opened
   character(*), intent(in) :: filename    ! name of file to be opened
   character(*), intent(in) :: stopfile    ! FE code stopped file


!  ----- local variables -----
   logical(4) :: opened,done
   integer(4) :: ict

   ict=0
   done=.false.
   opened=.false.
   do while ((ict.le.MAXLOOP).AND.(.NOT.done).AND.(.NOT.opened))
      ict=ict+1
      open(filenum,file=filename,err=100)
      opened=.true.
      cycle

100   continue
      if(ict.eq.1)write(6,*)'open waiting for file ',filename
        call sleep(waittime)
      inquire(file=stopfile,exist=done)
   enddo

   if(.not.opened)then
      call CouplerIO_Abort('Coupler SafeOpen :: timed out waiting for file '//filename)
   else if(done)then
      call CouplerIO_Abort('Coupler SafeOpen :: Fluid Stop File Found; file name is  '//stopfile)
   endif

   return
end subroutine SafeOpen

subroutine WaitRemoveFile(filename,stopfile)
!   Wait for file to be removed

!   ----- arguments -----
   character(*), intent(in) :: filename    ! name of file to be opened
   character(*), intent(in) :: stopfile    ! FE code stopped file

!   ----- local variable -----
   logical(4) :: here,done
   character(len=80) :: line
   integer(4) :: ict

   inquire(file=filename,exist=here)

   ict=0
   done=.false.
   do while((here).AND.(ict.le.MAXLOOP).AND.(.NOT.done))

      ict=ict+1
      inquire(file=filename,exist=here)
        call sleep(waittime)   !unix comment
      inquire(file=stopfile,exist=done)

   enddo

   if(here)then
      call CouplerIO_Abort('CouplerIO WaitRemoveFile :: timed out waiting for file '//filename)
   else if(done)then
      call CouplerIO_Abort('CouplerIO WaitRemoveFile ::  Stop File Found; file name is '//stopfile)
   endif

   return

end subroutine WaitRemoveFile

subroutine CouplerIO_Abort(string)
   character(*), intent(in) :: string

   write(6,'(a)')string
   stop 'CouplerIO_Abort'
      
end subroutine CouplerIO_Abort

end module


