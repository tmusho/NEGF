!
!    Copyright 2012 Terence Musho
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
!
!  sltc.f90 
!
!  Authors: G. Walker - Vanderbilt University
!           T. Musho  - Vanderbilt University
!           
! 
!  Program was originally written in MathLab(Octave/GNU) Converted to Fortran 90/95
!
!****************************************************************************
!
!  program: slpf
!
!  PURPOSE:  Calculate phonon transport characteristics for a strained Si/Ge/Si superlattice. 
!            Code was written to couple to NEGF Electron
!
!****************************************************************************
!
program sltc

    !Use global variables from module
    use Startup
    use Debug
    use Errors
    use Hamil
    use Eigen
    use Init     
    use Negf

implicit none

    !! Define Variables !!
    integer(kind=4) :: iargc, n, err, ni
    real(kind=8) :: t_end, t_start

    character*100 buffer ! Command Line Buffer
    character*100 :: op

    call MPI_INIT(err)
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,err)
    if(NPROC .eq. 1) then
       write(*,*) 'Must have >1 ncpu'
       stop 'abort'
    end if

    !! Title screen
    if(RANK .eq. 0) call Debug_Title

    if(rank .eq. 0) call Mpinterface_Read_Header(Lb,Lw,Lc,Nl,mat1,mat2,Lg,Lg_o,Tp,dT,debug_l)
    !send parameters to clients
    call MPI_Bcast(Lb,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
    call MPI_Bcast(Lw,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
    call MPI_Bcast(Lc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
    call MPI_Bcast(Nl,1,MPI_INTEGER,0,MPI_COMM_WORLD,err)
    call MPI_Bcast(mat1,1,MPI_INTEGER,0,MPI_COMM_WORLD,err)
    call MPI_Bcast(mat2,1,MPI_INTEGER,0,MPI_COMM_WORLD,err)
    call MPI_Bcast(Lg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
    call MPI_Bcast(Lg_o,1,MPI_INTEGER,0,MPI_COMM_WORLD,err)
    call MPI_Bcast(Tp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
    call MPI_Bcast(debug_l,1,MPI_INTEGER,0,MPI_COMM_WORLD,err)
    call MPI_Bcast(dT,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,err)
    call MPI_BARRIER(MPI_COMM_WORLD,err) 

    !! Start Time !!
    t_start = MPI_Wtime()

    !! Main execution sequence
    ! Setting up variables
    call Startup_Main !! Call startup module
    if(RANK .eq. 0) then
      110 format (A)
      write(6,110) '** Setting up variables'
      write(op,110) '** Setting up variables'
      call debug_output(op)
      call flush(6)
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,err) 

    ! Setting up variables
    call Grad_Main !! Call startup module
    if(RANK .eq. 0) then
      write(6,110) '** Setting up gradient'
      write(op,110) '** Setting up gradient'
      call debug_output(op)
      call flush(6)
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,err)

    ! Constructing Hamiltonian matrix
    if(RANK .eq. 0) then
      write(*,110) '** Constructing Hamiltonian matrix'
      write(op,110) '** Constructing Hamiltonian matrix'
      call debug_output(op)
      call flush(6)
    end if
    call Hamil_Main !! Call hamiltonian module
    call MPI_BARRIER(MPI_COMM_WORLD,err) 

    ! Calculating eigenvalues
    if(RANK .eq. 0) then
      write(*,110) '** Calculating eigenvalues'
      write(op,110) '** Calculating eigenvalues'
      call debug_output(op)
      call flush(6)
    end if
    call Eigen_Main !! Call eignvalues
    call MPI_BARRIER(MPI_COMM_WORLD,err)

    if(RANK .eq. 0) then
      write(*,110) '** Redefine Energy Range'
      write(op,110) '** Redefine Energy Range'
      call debug_output(op)
      call flush(6)
    end if
    if(usr .eq. 0 .and. RANK .eq. 0) call Startup_Energy_Redefine(DH)
    call MPI_BARRIER(MPI_COMM_WORLD,err)

    ! Initallizing Variables
    if(RANK .eq. 0) then
      write(*,110) '** Initializing matricies'
      write(op,110) '** Initializing matricies'
      call debug_output(op)
      call flush(6)
    end if
    call Init_Main !! Call initalization main loop
    call MPI_BARRIER(MPI_COMM_WORLD,err)

    ! Entering main calculation
    if(RANK .eq. 0) then
      write(*,110) '** Entering main calculation'
      write(op,110) '** Entering main calculation'
      call debug_output(op)
      call flush(6)
    end if
    call Negf_Main !! Call NEGF loop
    call MPI_BARRIER(MPI_COMM_WORLD,err)

    if(RANK .eq. 0) then
      write(*,110) '** Deinitializing matricies'
      write(op,110) '** Deinitializing matricies'
      call debug_output(op)
      call flush(6)
    end if
    call Init_Deallocate !! Deallocate dynamic memory

    !! End Time !!
    t_end = MPI_Wtime(); t_end = t_end-t_start

    if(rank .eq. 0) then
    112 format (A,E15.6,E15.6,E15.6,I3)
    111 format (A,E15.6,E15.6,E15.6,E15.6)
    120 format (A,F12.6,A)
    121 format (A,ES15.6,A)

    write(*,*) !Write to screen
    write(*,112) 'Lb, Lw, Lc, Nl = ',Lb,Lw,Lc,Nl
    write(*,111) 'Th, Tc, L, Jtot = ',Tp+dT/2,Tp-dT/2,Lt,Itot
    write(*,121) 'Computed Thermal Energy = ',abs(Qtot),' J/m^2'
    write(*,121) 'Computed Thermal Conductivity = ',abs(Itot)*Lt/dT,' W/K-m'
    write(*,121) 'Computed Thermal Conductance = ',abs(Itot)/dT,' W/K-m^2'
    write(*,121) 'Analytic Thermal Conductance = ',pi**2*kbm**2*Tp/(3*hbar*2*pi)*NE,' W/K'
    !write(*,121) 'Average Diffusivity = ',Df,' ?'
    !write(*,121) 'Average Specific Heat = ',Cv,' ?'
    write(*,120) 'Elapse Time (time) = ',t_end,' sec'
    write(*,*)

    write(op,*) '--------------------------------------------------'  !Write to file
    call debug_output(op)
    write(op,111) 'Th, Tc, L, Jtot = ',Tp+dT/2,Tp-dT/2,Lt,Itot
    call debug_output(op)   
    write(op,111) 'Thermal Conductivity = ',abs(Itot)*Lt/dT
    call debug_output(op)
    write(op,111) 'Thermal Conductance = ',Itot/dT
    call debug_output(op)
    write(op,111) 'Thermal Power = ',abs(Itot)
    call debug_output(op)
    write(op,111) 'Thermal Energy = ',abs(Qtot)
    call debug_output(op)
    write(op,111) 'Average Temperature = ',Tp
    call debug_output(op)
    !write(op,111) 'Average Diffusivity = ',Df
    !call debug_output(op)
    !write(op,111) 'Average Specific Heat = ',Cv
    !call debug_output(op)
    write(op,120) 'Elapse Time (time) = ',t_end,' sec' !Write to file
    call debug_output(op)
    end if

    call MPI_FINALIZE(err)

end program sltc
