!  ReadInput.f90 
!
!  Authors: G. Walker - Vanderbilt University
!           T. Musho - Vanderbilt University
!
!  Created: 1/20/11
!  Last Modified: 1/20/11
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
!  PROGRAM: InputFormat(Module)
!
!  PURPOSE: Checks Input Format
!
!  SUBROUTINES: 
!              
!            
!
!****************************************************************************
!
module ReadInput
use Utility, only :  UnitOut, &
                     UnitIn, &
                     StdOut, &
                     UnitLog, &
                     LINE_LENGTH, &
                     FluidID, &
                     Utility_ErrorMsg
use Input, only : Input_FindSection, &
                  Input_GetLine, &
                  Input_FindKeyWord, &
                  Input_Iread,   &
                  Input_Fread
use Debug
implicit none

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ReadInput_Case - Read Case
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadInput_Case(Title,coupled)    

!  ----- Argument Lines -----iprotect
   character(40), intent(out) :: Title             ! initial step
   integer, intent(out) :: coupled

!  ----- Local variables -----
   integer(4), parameter :: NUM_INT_LINES=3
   integer(4) :: NumLines,LineNoStart
   character(LINE_LENGTH) :: NewLine

   if(FluidID==0)then
      !write(StdOut,*)
      write(StdOut,*)'+++ Reading Case +++'
   endif   

   call Input_FindSection(.true.,'<case>','<end case>',NUM_INT_LINES, &
                          NUM_INT_LINES,NumLines,LineNoStart,.true.)

   NewLine    =Input_GetLine(LineNoStart+2)
   Title      =NewLine(1:40)
   NewLine    =Input_GetLine(LineNoStart+3)

   if(index(NewLine,'pf')/=0)then
      coupled=1
   else if(index(NewLine,'iv')/=0)then
      coupled=2
   else
      call Utility_ErrorMsg('ReadInput_Case:: run type must be either PF or IV')
   endif

   write(UnitLog,*)'+++ Case Read +++'   
 
   return
end subroutine ReadInput_Case

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ReadInput_Integration - Read Integration
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadInput_Integration(amr_o,amr_d,cop,coe,sc_dv,sc_sb,sc_po,sc_ip,NE,Eo,Ef,an,sc_sc,sc_tp,sc_num)    

!  ----- Argument Lines -----
   real(8), intent(out) :: sc_dv,sc_sb,sc_po,sc_ip,sc_sc,sc_tp
   real(8), intent(out) :: an, Eo, Ef
   integer(4), intent(out) :: amr_o, amr_d, sc_num            
   integer(4), intent(out) :: cop,coe 
   integer(4), intent(out) :: NE  

!  ----- Local variables -----
   integer(4), parameter :: MIN_INT_LINES=15,MAX_INT_LINES=15,LARGEINT=1000000
   integer(4) :: NumLines,LineNoStart,LinePtr

   if(FluidID==0)then
      !write(StdOut,*)
      write(StdOut,*)'+++ Reading Integration +++'
   endif  

   call Input_FindSection(.true.,'<integration>','<end integration>',MIN_INT_LINES, &
                          MAX_INT_LINES,NumLines,LineNoStart,.true.)

   amr_o    =Input_Iread(LineNoStart+1,'AMR        ',0,1) 
   amr_d    =Input_Iread(LineNoStart+2,'AMR Div    ',1,1000)
   cop      =Input_Iread(LineNoStart+3,'Potnt Type ',0,1)
   coe      =Input_Iread(LineNoStart+4,'Int Cutoff ',0,1) 
   sc_tp    =Input_Fread(LineNoStart+5,'Tp Conv  ',0.D0,1.D0)     
   sc_dv    =Input_Fread(LineNoStart+6,'Volt Conv  ',0.D0,1.D0)   
   sc_sb    =Input_Fread(LineNoStart+7,'Subbd Conv ',0.D0,1.D0)
   sc_po    =Input_Fread(LineNoStart+8,'Potnt Conv ',0.D0,1.D0)
   sc_ip    =Input_Fread(LineNoStart+9,'Init Pot   ',0.D0,1.D0)
   sc_sc    =Input_Fread(LineNoStart+10,'Scat Conv   ',0.D0,1.D0)
   NE       =Input_Iread(LineNoStart+11,'Num E      ',0,LARGEINT)
   Eo       =Input_Fread(LineNoStart+12,'Ecut Low   ',-1.D0,1.D0)
   Ef       =Input_Fread(LineNoStart+13,'Ecut High  ',-1.D0,1.D0)
   an       =Input_Fread(LineNoStart+14,'Cell Size  ',0.D0,1.D0)
   sc_num   =Input_Iread(LineNoStart+15,'Sub Num    ',0,LARGEINT)
 
   write(UnitLog,*)'+++ Integration Read +++'  

   return
end subroutine ReadInput_Integration

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ReadInput_Potential - Read Potential
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadInput_potential(coc,smix,wmix,nank)    

!  ----- Argument Lines -----
   real(8), intent(out) :: smix,wmix
   integer(4), intent(out) :: coc,nank            

!  ----- Local variables -----
   integer(4), parameter :: MIN_INT_LINES=4,MAX_INT_LINES=4,LARGEINT=1000000
   integer(4) :: NumLines,LineNoStart,LinePtr

   if(FluidID==0)then
      !write(StdOut,*)
      write(StdOut,*)'+++ Reading Potential +++'
   endif  

   call Input_FindSection(.true.,'<potential>','<end potential>',MIN_INT_LINES, &
                          MAX_INT_LINES,NumLines,LineNoStart,.true.)

   coc      =Input_Iread(LineNoStart+1,'Potnt Type ',0,2) 
   smix     =Input_Fread(LineNoStart+2,'Simple Mix ',0.D0,10.D0)
   wmix     =Input_Fread(LineNoStart+3,'Andr Mix   ',0.D0,100.D0)
   nank     =Input_Iread(LineNoStart+4,'Andr Prev  ',0,100)   
 
   write(UnitLog,*)'+++ Reading Potential +++'  

   return
end subroutine ReadInput_potential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ReadInput_strain - Read Strain
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadInput_strain(sct)

!  ----- Argument Lines -----
   integer(4), intent(out) :: sct

!  ----- Local variables -----
   integer(4), parameter :: MIN_INT_LINES=1,MAX_INT_LINES=1,LARGEINT=1000000
   integer(4) :: NumLines,LineNoStart,LinePtr

   if(FluidID==0)then
      !write(StdOut,*)
      write(StdOut,*)'+++ Reading Strain +++'
   endif

   call Input_FindSection(.true.,'<strain>','<end strain>',MIN_INT_LINES, &
                          MAX_INT_LINES,NumLines,LineNoStart,.true.)

   sct      =Input_Iread(LineNoStart+1,'Strain    ',0,1)

   write(UnitLog,*)'+++ Reading Strain +++'

   return
end subroutine ReadInput_strain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ReadInput_temp - Read Temperature
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadInput_temp(Tp,dT)

!  ----- Argument Lines -----
   real(8), intent(out) :: Tp,dT

!  ----- Local variables -----
   integer(4), parameter :: MIN_INT_LINES=2,MAX_INT_LINES=2,LARGEINT=1000000
   integer(4) :: NumLines,LineNoStart,LinePtr

   if(FluidID==0)then
      !write(StdOut,*)
      write(StdOut,*)'+++ Reading Temperature +++'
   endif

   call Input_FindSection(.true.,'<temperature>','<end temperature>',MIN_INT_LINES, &
                          MAX_INT_LINES,NumLines,LineNoStart,.true.)

   Tp      =Input_Fread(LineNoStart+1,'Avg Temp    ',0.D0,2500.D0)
   dT      =Input_Fread(LineNoStart+2,'Diff Temp    ',-2500.D0,2500.D0)


   write(UnitLog,*)'+++ Reading Temperature +++'

   return
end subroutine ReadInput_temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ReadInput_restart - Read Restart
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadInput_restart(tor,rsr,rss)

!  ----- Argument Lines -----
   real(8), intent(out) :: rsr
   integer(4), intent(out) :: tor,rss

!  ----- Local variables -----
   integer(4), parameter :: MIN_INT_LINES=3,MAX_INT_LINES=3,LARGEINT=1000000
   integer(4) :: NumLines,LineNoStart,LinePtr

   if(FluidID==0)then
      !write(StdOut,*)
      write(StdOut,*)'+++ Reading Restart +++'
   endif

   call Input_FindSection(.true.,'<restart>','<end restart>',MIN_INT_LINES, &
                          MAX_INT_LINES,NumLines,LineNoStart,.true.)

   tor =Input_Iread(LineNoStart+1,'Typ Rst     ',0,3)
   rsr =Input_Fread(LineNoStart+2,'Tme Rst     ',0.D0,100.D8)
   rss =Input_Iread(LineNoStart+3,'Stp Rst     ',0,LARGEINT)

   write(UnitLog,*)'+++ Reading Restart +++'

   return
end subroutine ReadInput_restart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ReadInput_Geometry - Read Geometry
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadInput_geometry(Nl,Lg,Lg_o,Lc)

!  ----- Argument Lines -----
   real(8), intent(out) :: Lg,Lc
   integer(4), intent(out) :: Nl, Lg_o

!  ----- Local variables -----
   integer(4), parameter :: MIN_INT_LINES=4,MAX_INT_LINES=4,LARGEINT=1000000
   integer(4) :: NumLines,LineNoStart,LinePtr

   if(FluidID==0)then
      !write(StdOut,*)
      write(StdOut,*)'+++ Reading Geometry +++'
   endif

   call Input_FindSection(.true.,'<geometry>','<end geometry>',MIN_INT_LINES, &
                          MAX_INT_LINES,NumLines,LineNoStart,.true.)

   Nl       =Input_Iread(LineNoStart+1,'Num Layr    ',0,LARGEINT)
   Lg       =Input_Fread(LineNoStart+2,'Grading     ',0.D0,100.D0)
   Lg_o     =Input_Iread(LineNoStart+3,'Mat Grad    ',1,2)
   Lc       =Input_Fread(LineNoStart+4,'Contct Leng ',0.D0,100.D0)

   write(UnitLog,*)'+++ Reading Geometry +++'

   return
end subroutine ReadInput_geometry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ReadInput_Solver - Read Solver
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadInput_solver(csp)

!  ----- Argument Lines -----
   integer(4), intent(out) :: csp

!  ----- Local variables -----
   integer(4), parameter :: MIN_INT_LINES=1,MAX_INT_LINES=1,LARGEINT=1000000
   integer(4) :: NumLines,LineNoStart,LinePtr

   if(FluidID==0)then
      !write(StdOut,*)
      write(StdOut,*)'+++ Reading Solver +++'
   endif

   call Input_FindSection(.true.,'<solver>','<end solver>',MIN_INT_LINES, &
                          MAX_INT_LINES,NumLines,LineNoStart,.true.)

   csp       =Input_Iread(LineNoStart+1,'Solver Type    ',1,4)

   write(UnitLog,*)'+++ Reading Solver +++'

   return
end subroutine ReadInput_solver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ReadInput_material1 - Read Material
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadInput_material1(mat,Lb,Nd_Si,Nc_Si,epsilrS,ms_e,a0Si,D001Si,G001Si,EdUSi,Eph_si,op_cut_si)

!  ----- Argument Lines -----
   real(8), intent(out) :: Lb,Nd_Si,Nc_Si,epsilrS,ms_e,a0Si,D001Si,G001Si,EdUSi
   real(8), intent(out) :: Eph_si, op_cut_si
   integer(4), intent(out) :: mat

!  ----- Local variables -----
   real(8), parameter :: LARGEREAL = 1E30
   integer(4), parameter :: MIN_INT_LINES=12,MAX_INT_LINES=12,LARGEINT=1000000
   integer(4) :: NumLines,LineNoStart,LinePtr

   if(FluidID==0)then
      !write(StdOut,*)
      write(StdOut,*)'+++ Reading Materials 1 +++'
   endif

   call Input_FindSection(.true.,'<material1>','<end material1>',MIN_INT_LINES, &
                          MAX_INT_LINES,NumLines,LineNoStart,.true.)

   mat      =Input_Iread(LineNoStart+1,'Mat1 ',1,2)
   Lb       =Input_Fread(LineNoStart+2,'Layr Thck  ',0.D0,100.D0)
   Nd_Si    =Input_Fread(LineNoStart+3,'Dopnt Conv ',0.D0,LARGEREAL)
   Nc_Si    =Input_Fread(LineNoStart+4,'Effect Dopnt   ',0.D0,LARGEREAL)
   epsilrS  =Input_Fread(LineNoStart+5,'Rel Perm  ',0.D0,100.D0)
   ms_e     =Input_Fread(LineNoStart+6,'Effct Mass  ',0.D0,100.D0)
   a0Si     =Input_Fread(LineNoStart+7,'Equi Latt  ',0.D0,100.D0)
   D001Si   =Input_Fread(LineNoStart+8,'Def Strn Pot  ',0.D0,100.D0)
   G001Si   =Input_Fread(LineNoStart+9,'Shr Mod  ',0.D0,100.D0)
   EdUSi    =Input_Fread(LineNoStart+10,'Def Pot  ',0.D0,100.D0)
   Eph_si    =Input_Fread(LineNoStart+11,'Eph  ',0.D0,100.D0)
   op_cut_si    =Input_Fread(LineNoStart+12,'Eph Op Cut  ',0.D0,100.D0)


   write(UnitLog,*)'+++ Reading Material 1 +++'

   return
end subroutine ReadInput_material1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ReadInput_material2 - Read Material
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadInput_material2(mat,Lb,Nd_Ge,Nc_Ge,epsilrG,ms_e,a0Ge,D001Ge,G001Ge,EdUGe,Eph_ge,op_cut_ge)

!  ----- Argument Lines -----
   real(8), intent(out) :: Lb,Nd_Ge,Nc_Ge,epsilrG,ms_e,a0Ge,D001Ge,G001Ge,EdUGe
   real(8), intent(out) :: Eph_ge, op_cut_ge
   integer(4), intent(out) :: mat

!  ----- Local variables -----
   real(8), parameter :: LARGEREAL = 1E30
   integer(4), parameter :: MIN_INT_LINES=12,MAX_INT_LINES=12,LARGEINT=1000000
   integer(4) :: NumLines,LineNoStart,LinePtr

   if(FluidID==0)then
      !write(StdOut,*)
      write(StdOut,*)'+++ Reading Materials 2 +++'
   endif

   call Input_FindSection(.true.,'<material2>','<end material2>',MIN_INT_LINES, &
                          MAX_INT_LINES,NumLines,LineNoStart,.true.)

   mat      =Input_Iread(LineNoStart+1,'Mat2 ',1,2)
   Lb       =Input_Fread(LineNoStart+2,'Layr Thck  ',0.D0,100.D0)
   Nd_Ge    =Input_Fread(LineNoStart+3,'Dopnt Conv ',0.D0,LARGEREAL)
   Nc_Ge    =Input_Fread(LineNoStart+4,'Effect Dopnt   ',0.D0,LARGEREAL)
   epsilrG  =Input_Fread(LineNoStart+5,'Rel Perm  ',0.D0,100.D0)
   ms_e     =Input_Fread(LineNoStart+6,'Effct Mass  ',0.D0,100.D0)
   a0Ge     =Input_Fread(LineNoStart+7,'Equi Latt  ',0.D0,100.D0)
   D001Ge   =Input_Fread(LineNoStart+8,'Def Strn Pot  ',0.D0,100.D0)
   G001Ge   =Input_Fread(LineNoStart+9,'Shr Mod  ',0.D0,100.D0)
   EdUGe    =Input_Fread(LineNoStart+10,'Def Pot  ',0.D0,100.D0)
   Eph_ge    =Input_Fread(LineNoStart+11,'Eph  ',0.D0,100.D0)
   op_cut_ge    =Input_Fread(LineNoStart+12,'Eph Op Cut  ',0.D0,100.D0)


   write(UnitLog,*)'+++ Reading Material 2 +++'

   return
end subroutine ReadInput_material2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ReadInput_Debug - Read Debug
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadInput_Debug(NumFluidProcs,nstartOld,FluidID)

   integer(4), intent(in) :: NumFluidProcs
   integer(4), intent(in) :: nstartOld
   integer(4), intent(in) :: FluidID
   integer(4), parameter :: MIN_INT_LINES=1,MAX_INT_LINES=1
   integer(4) :: NumLines,LineNoStart,i
   character(LINE_LENGTH) :: NewLine
   character(80) :: Outfile

   if(FluidID==0)then
      !write(StdOut,*)
      write(StdOut,*)'+++ Reading Debug Control +++'
   endif 
 
   call Input_FindSection(.true.,'<debug control>','<end debug control>',MIN_INT_LINES,MAX_INT_LINES,NumLines,LineNoStart,.true.)

   debug_l      =Input_Iread(LineNoStart+1,'Debug Level ',1,10)

   !if(debug_l .gt. 5)then
   !  open Output File
   !   Outfile=Utility_MakeFileName('debug',NumFluidProcs,nstartOld,FluidID,'out')   
   !   open(UnitDebug,file=OutFile)
   !   write(UnitDebug,'(30x,a)')version
   !endif

   write(UnitLog,*)'+++ Reading Debug Control +++'

end subroutine ReadInput_Debug

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ReadInput_Voltage - Read Voltage
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadInput_voltage(NV,Vo,Vf,Vneg,Vpos,Ineg,Ipos,Vmax,Imax,Imin,Nb)

!  ----- Argument Lines -----
   real(8), intent(out) :: Vo,Vf,Vneg,Vpos,Ineg,Ipos,Vmax,Imax,Imin
   integer(4), intent(out) :: NV, Nb

!  ----- Local variables -----
   real(8), parameter :: LARGEREAL = 1E30
   integer(4), parameter :: MIN_INT_LINES=11,MAX_INT_LINES=11,LARGEINT=1000000
   integer(4) :: NumLines,LineNoStart,LinePtr

   if(FluidID==0)then
      !write(StdOut,*)
      write(StdOut,*)'+++ Reading Voltage +++'
   endif

   call Input_FindSection(.true.,'<voltage>','<end voltage>',MIN_INT_LINES, &
                          MAX_INT_LINES,NumLines,LineNoStart,.true.)

   NV       =Input_Iread(LineNoStart+1,'NV  ',0,1000)
   Vo       =Input_Fread(LineNoStart+2,'Vo  ',-LARGEREAL,LARGEREAL)
   Vf       =Input_Fread(LineNoStart+3,'Vf  ',-LARGEREAL,LARGEREAL)
   Vneg     =Input_Fread(LineNoStart+4,'Vneg   ',-LARGEREAL,LARGEREAL)
   Vpos     =Input_Fread(LineNoStart+5,'Vpos  ',-LARGEREAL,LARGEREAL)
   Ineg     =Input_Fread(LineNoStart+6,'Ineg  ',-LARGEREAL,LARGEREAL)
   Ipos     =Input_Fread(LineNoStart+7,'Ipos  ',-LARGEREAL,LARGEREAL)
   Vmax     =Input_Fread(LineNoStart+8,'Vmax  ',-LARGEREAL,LARGEREAL)
   Imax     =Input_Fread(LineNoStart+9,'Imax  ',-LARGEREAL,LARGEREAL)
   Imin     =Input_Fread(LineNoStart+10,'Imin  ',-LARGEREAL,LARGEREAL)
   Nb       =Input_Iread(LineNoStart+11,'Nb  ',0,10000)

   write(UnitLog,*)'+++ Reading Voltage +++'

   return
end subroutine ReadInput_voltage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ReadInput_scatter - Read Scatter
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadInput_scatter(tst,psc,coc_s,smix_s,wmix_s,Nank_s)

!  ----- Argument Lines -----
   real(8), intent(out) :: smix_s, wmix_s
   integer(4), intent(out) :: tst,psc,coc_s,Nank_s

!  ----- Local variables -----
   real(8), parameter :: LARGEREAL = 1E30
   integer(4), parameter :: MIN_INT_LINES=6,MAX_INT_LINES=6,LARGEINT=1000000
   integer(4) :: NumLines,LineNoStart,LinePtr

   if(FluidID==0)then
      !write(StdOut,*)
      write(StdOut,*)'+++ Reading Scattering +++'
   endif

   call Input_FindSection(.true.,'<scattering>','<end scattering>',MIN_INT_LINES, &
                          MAX_INT_LINES,NumLines,LineNoStart,.true.)

   tst       =Input_Iread(LineNoStart+1,'tst  ',1,8)
   psc       =Input_Iread(LineNoStart+2,'psc  ',0,1)
   coc_s     =Input_Iread(LineNoStart+3,'coc_s  ',0,LARGEINT)
   smix_s     =Input_Fread(LineNoStart+4,'smix_s   ',-LARGEREAL,LARGEREAL)
   wmix_s     =Input_Fread(LineNoStart+5,'wmix_s  ',-LARGEREAL,LARGEREAL)
   Nank_s     =Input_Iread(LineNoStart+6,'Nank_s   ',0,LARGEINT)
   write(UnitLog,*)'+++ Reading Scattering +++'

   return
end subroutine ReadInput_scatter
end module
