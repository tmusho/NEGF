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
!  slpf.f90 
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
!  PURPOSE:  Calculate current density vs. voltage characteristics for a strained Si/Ge/Si superlattice. 
!            Ballistic transport is assumed. Current is summed over conduction band and 21 subbands.
!            This version is parallelized by splitting the integration loop across procs
!
!****************************************************************************
!
program slpf

    !Use global variables from module
    use Startup
    use Debug
    use Errors
    use Device
    use Hamil
    use Eigen
    use Init     
    use Negf
    use Mpinterface
    use Negf_Scatter

implicit none

    !! Define Variables !!
    integer(kind=4) :: iargc, n, err, ni
    real(kind=8) :: t_end, t_start

    character*120 buffer ! Command Line Buffer
    character*120 :: op

    call MPI_INIT(err)
    call MPI_COMM_RANK(MPI_COMM_WORLD,RANK,err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,err)
    if(NPROC .eq. 1) then
       write(*,*) 'Must have >1 ncpu'
       stop 'abort'
    end if

    !! Read Command Line Args !!
    n = iargc()
    if(n .ne. 0) then
    call getarg(1,buffer); ni = 1 ! Read first command line determine debug value
    if(buffer .eq. '-c') then !process command line arguments
       itype = 1; ni = ni + 1
       call getarg(ni,buffer) ! Get debug values from command line
       if(buffer .eq. '-d' .and. n .eq. 10) then ! Check for debug flag
         ni = ni + 1
         call getarg(ni,buffer) ! Get debug values from command line
         read(buffer,*) debug_l; ni = ni + 1
         call getarg(ni,buffer)! Read Lb
       else if(n .ne. 8) then
         write(*,*) 'slpf [-d debug_level] Bar_Len(m) Well_Len(m) Con_Len(m) Num_Lay(#) Dop_Con(cm-3) &
                   & mat1(#) mat2(#)'
         stop 'Command Line Arguments Error'
       end if
      
       read(buffer,*) Lb; ni = ni + 1
       call getarg(ni,buffer)
       read(buffer,*) Lw; ni = ni + 1
       call getarg(ni,buffer)
       read(buffer,*) Lc; ni = ni + 1
       call getarg(ni,buffer)
       read(buffer,*) Nl; ni = ni + 1
       call getarg(ni,buffer)
       read(buffer,*) Nd_Ge; ni = ni + 1
       call getarg(ni,buffer)
       read(buffer,*) mat1; ni = ni + 1
       call getarg(ni,buffer)
       read(buffer,*) mat2; ni = ni + 1
    else if(buffer .eq. '-i' .and. n .eq. 2) then !process input script
      itype = 2; ni = ni + 1
      call getarg(ni,buffer)
      read(buffer,*) input_fname; ni = ni + 1
    else if(buffer .eq. '-r' .and. n .eq. 2) then ! Check for restart flag
      itype = 2; ni = ni + 1
      call getarg(ni,buffer) ! Get debug values from command line
      read(buffer,*) restart_fname; ni = ni + 1
      rst = 1
    else if(n .ne. 0) then
      write(*,*) 'Usage: slpf [OPTION...] [FILE]...'
      write(*,*) 'slpf calculates the electrical transport of superlattice thermoelectric devices &
                  & using a non-equilibrium Greens function formalism'
      write(*,*) ' '
      write(*,*) 'Examples:'
      write(*,*) '  slpf -c -d 1 1e-9 1e-9 0 1 1e18 1 2'
      write(*,*) '  slpf -i bi2te3.inp'
      write(*,*) '  slpf'
      write(*,*) ' '
      write(*,*) ' Main input modes:'
      write(*,*) ' -c       command line input'
      write(*,*) ' -d       debug level [1-10] 10=lots'
      write(*,*) '   ARG1       thickness material 1'
      write(*,*) '   ARG2       thickness material 2'
      write(*,*) '   ARG3       thickness contacts'
      write(*,*) '   ARG4       number of bilayers'
      write(*,*) '   ARG5       dopant concentration'
      write(*,*) '   ARG6       material 1 selection [1or2]'
      write(*,*) '   ARG6       material 2 selection [1or2]'
      write(*,*) ' '
      write(*,*) ' -i       user specified input file'
      write(*,*) ' -r       user specified restart file' 
      write(*,*) ' '
      write(*,*) 'Report bugs to <greg.walker@vanderbilt.edu>.'
      write(*,*) ' '
      stop 'Command Line Arguments Error'
    end if
    else
      itype = 2
      if(RANK .eq. 0) write(6,*) 'Using default input script slpf.inp. To specify a different script use -i flag.'
    end if

    !! Start Time !!
    t_start = MPI_Wtime()

    !! Main execution sequence
    ! Setting up variables
    call Startup_Main !! Call startup module
    if(RANK .eq. 0) then
      110 format (A)
    write(6,110) ' +++ Setting up variables +++'
    write(op,110) ' +++ Setting up variables +++'
      call debug_output(op)
      call flush(6)
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,err) 

    ! Constructing device templates
    if(RANK .eq. 0) then
    write(*,110) ' +++ Constructing device templates +++'
    write(op,110) ' +++ Constructing device templates +++'
      call debug_output(op)
      call flush(6)
    end if
    call Device_Main !! Call device module
    call MPI_BARRIER(MPI_COMM_WORLD,err) 

    ! Constructing Hamiltonian matrix
    if(RANK .eq. 0) then
    write(*,110) ' +++ Constructing Hamiltonian matrix +++'
    write(op,110) ' +++ Constructing Hamiltonian matrix +++'
      call debug_output(op)
      call flush(6)
    end if
    call Hamil_Main !! Call hamiltonian module
    call MPI_BARRIER(MPI_COMM_WORLD,err) 

    ! Calculating eigenvalues
    if(RANK .eq. 0) then
    write(*,110) ' +++ Calculating eigenvalues +++'
    write(op,110) ' +++ Calculating eigenvalues +++'
      call debug_output(op)
      call flush(6)
    end if
    call Eigen_Main !! Call eignvalues
    call MPI_BARRIER(MPI_COMM_WORLD,err) 

    ! Grap Phonon Data, Phonon code runs along side this exe
    if(RANK .eq. 0) then
      if(tst .lt. 6) then
        call Mpinterface_Read_Eig(DH_p,debug_l)
        call Mpinterface_Read_Fermi(f_p,NE_p,debug_l)
        call Mpinterface_Read_Gn(Gn_p,NE_p,Eo_p,Ef_p,Np_p,Nb_p,debug_l,chng_tprt)
        call Mpinterface_Read_Gp(Gp_p,NE_p,Np_p,Nb_p,debug_l,chng_tprt)
        call Mpinterface_Clean
        call MPI_BARRIER(MPI_COMM_WORLD,err)
        call MPI_Bcast(Nb_p,1,MPI_INTEGER,0,MPI_COMM_WORLD,err) !send number of phonons to clients
      else
        Nb_p = 1 !single phonon freq
        !call Mpinterface_killph
      end if
    else
      Nb_p = 1 !single phonon freq
      if(tst .lt. 6) then
        call MPI_BARRIER(MPI_COMM_WORLD,err)
        call MPI_Bcast(Nb_p,1,MPI_INTEGER,0,MPI_COMM_WORLD,err) !send number of phonons to clients
      end if
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,err) !sync

    ! Calculating eigenvalues
    if(RANK .eq. 0) then
    write(*,110) ' +++ Initializing matricies +++'
    write(op,110) ' +++ Initializing matricies +++'
      call debug_output(op)
      call flush(6)
    end if
    call Init_Main !! Call initalization main loop
    call MPI_BARRIER(MPI_COMM_WORLD,err)

    ! Entering main calculation
    if(RANK .eq. 0) then
    write(*,110) ' +++ Entering main calculation +++'
    write(op,110) ' +++ Entering main calculation +++'
      call debug_output(op)
      call flush(6)
    end if
    call Negf_Main !! Call NEGF loop
    call MPI_BARRIER(MPI_COMM_WORLD,err)

    if(RANK .eq. 0) then
      write(*,110) ' +++ Deallocate matricies +++'
      write(op,110) ' +++ Deallocate matricies +++'
      call debug_output(op)
      call flush(6)
    end if
    call Init_Deallocate !! Deallocate dynamic memory

    !! End Time !!
    t_end = MPI_Wtime(); t_end = t_end-t_start

    if(rank .eq. 0) then
    112 format (A,E15.6,E15.6,E15.6,I3,E15.6)
    111 format (A,E15.6,E15.6,E15.6,E15.6)
    113 format (A,E15.6,E15.6,E15.6,E15.6,E15.6)
    120 format (A,F12.6,A)
    121 format (A,E15.6,A)

    write(*,*) !Write to screen
    write(*,112) 'Lb, Lw, Lc, Nl, Lt = ',Lb,Lw,Lc,Nl,Lt
    write(*,111) 'Th, Tc, V, Itot = ',Tp+dT/2,Tp-dT/2,V,Itot
    write(*,113) 'S, sigma, PF, Kel, time = ',S,sigma,PF,ke,t_end
    write(*,121) 'Lorenz Number = ',ke/(sigma*Tp),' W-Ohm/K^2'
    write(*,121) 'Overall Mobility = ',mob,' m^2/V-s'
    write(*,*)

    write(op,*) '--------------------------------------------------'  !Write to file
    call debug_output(op)    
    write(op,120) 'Seebeck Coefficient (S) = ',S,' uV/K' !Write to file
    call debug_output(op)
    write(op,121) 'Electrical Conductivity (sigma) = ',sigma,' 1/ohm-m' !Write to file
    call debug_output(op)
    write(op,121) 'Power Factor (PF) = ',PF,' uV^2/K^2-ohm-m' !Write to file
    call debug_output(op)
    write(op,121) 'Thermal Electrical Conductivity (Kel) = ',ke,' W/m-K' !Write to file
    call debug_output(op)
    write(op,121) 'Lorenz Number (L=Kel/sigmaT) = ',ke/(sigma*Tp),' W-Ohm/K^2' !Write to file
    call debug_output(op)
    write(op,121) 'Overall Mobility = ',mob,' m^2/V-s' !Write to file
    call debug_output(op)
    write(op,120) 'Elapse Time (time) = ',t_end,' sec' !Write to file
    call debug_output(op)
    end if

    call MPI_FINALIZE(err)

    end program slpf
