!
! Dean Townsley 2008
!

#include "constants.h"

subroutine Burn_init()

  use Burn_data

  use Driver_interface, ONLY : Driver_abort, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stampMessage
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use IO_interface, only : IO_getPrevScalar

  implicit none

  integer, parameter      :: infile_unit=2
  integer, parameter      :: spec_num=238  ! number of lines in SpeciesList.txt file
  character(len=4)        :: isotopeName
  real                    :: abar,zbar,bindEnergy,spinfactor
  integer                 :: isotope, nfound, nread
  integer                 :: nprocs

  logical                 :: detIgnFile
  real                    :: calib
  integer                 :: istat, i, error

  logical, save           :: restart

  call Driver_getMype(MESH_COMM, bn_meshMe)
  call Driver_getComm(MESH_COMM, bn_meshComm)
  call Driver_getNumProcs(MESH_COMM, bn_meshNumProcs)

  ! standard burn parameters (not specific to parametric burning)
  call RuntimeParameters_get("useBurn", bn_useBurn)
  call RuntimeParameters_get("useShockBurn", bn_useShockBurn)
  call RuntimeParameters_get("enucDtFactor", bn_enucDtFactor)


  call RuntimeParameters_get("bn_thermalReact", bn_thermalReact)
  call RuntimeParameters_get("bn_thermalReactInFlameThreshold", bn_thermalReactInflameThreshold)

  call RuntimeParameters_get("restart", restart)
  if (restart) then
     call IO_getPrevScalar("bn_neutLoss", bn_neutLoss, error)
     if (error == NOTFOUND) bn_neutLoss = 0.0
  else
     bn_neutLoss = 0.0
  endif
  bn_neutLossThisProcStep = 0.0

  !-------------
  !  determine whether we are detonating manually or automatically
  call RuntimeParameters_get("bn_detIgnFile", detIgnFile)
  call RuntimeParameters_get("bn_detIgnFileName", bn_detIgnFileName)
  call RuntimeParameters_get("bn_autoDDT", bn_autoDDT)

  if (detIgnFile .and. bn_autoDDT) then
     call Driver_abort("Cannot detonate manually and automatically!")
  endif

  !  read detonation ignition points from file
  if (detIgnFile) then

     call Logfile_stampMessage('[Burn/parametric init] Reading detonation ignition points from file')
     open(unit=infile_unit,file=bn_detIgnFileName,status='OLD',iostat=istat)
     if (istat /= 0) call Driver_abort("Unable to open detonation ignition points file")

     ! eat header
     read(infile_unit,*)
     read(infile_unit,*) pbIgnNum
     allocate(pbIgnTime(pbIgnNum),STAT=istat)
        if (istat /= 0) call Driver_abort("Cannot allocate pbIgnTime in paraburn init")
     allocate(pbIgnX(pbIgnNum),STAT=istat)
        if (istat /= 0) call Driver_abort("Cannot allocate pbIgnX in paraburn init")
     allocate(pbIgnY(pbIgnNum),STAT=istat)
        if (istat /= 0) call Driver_abort("Cannot allocate pbIgnY in paraburn init")
     allocate(pbIgnZ(pbIgnNum),STAT=istat)
        if (istat /= 0) call Driver_abort("Cannot allocate pbIgnZ in paraburn init")
     allocate(pbIgnR(pbIgnNum),STAT=istat)
        if (istat /= 0) call Driver_abort("Cannot allocate pbIgnR in paraburn init")
     do i = 1, pbIgnNum
        read(infile_unit,*) pbIgnTime(i), pbIgnX(i), pbIgnY(i), pbIgnZ(i), pbIgnR(i)
     enddo
     if (bn_meshMe == MASTER_PE) then
        write(*,*)  'read detonation ignition points t x y z'
        do i = 1, pbIgnNum
           write(*,*) pbIgnTime(i), pbIgnX(i), pbIgnY(i), pbIgnZ(i), pbIgnR(i)
        enddo
     endif
     close(unit=infile_unit)

  else if (bn_autoDDT) then

     if (.not. bn_thermalReact) then
        call Logfile_stampMessage('[Burn_init] bn_thermalReact is currently false but bn_autoDDT is true.')
        call Logfile_stampMessage('[Burn_init] bn_thermalReact will be reset to true&
            & when  detonation conditions are about to be met')
     endif


     call RuntimeParameters_get("pbIgnRho", pbIgnRho)
     call RuntimeParameters_get("pbIgnRhoCalib", calib)
     pbIgnRho = exp( log(pbIgnRho) - calib )
     call RuntimeParameters_get("pbIgnRhoFact", pbIgnRhoFact)
     if (pbIgnRhoFact >= 1.0)  &
        call Driver_abort("pbIgnRhoFact must be less than 1")
     pbIgnRhoFact = 10.e0**( pbIgnRhoFact * log10(pbIgnRho) ) 
     call RuntimeParameters_get("pbIgnPhfa", pbIgnPhfa)
     call RuntimeParameters_get("pbIgnDist", pbIgnDist)
     call RuntimeParameters_get("pbIgnRad", pbIgnRad)
     call RuntimeParameters_get("pbIgnSep", pbIgnSep)
     call RuntimeParameters_get("pbIgnNumMax", pbIgnNumMax)

     ! allocate global arrays for detonation points
     pbIgnNum = 0
     allocate(pbIgnTime(pbIgnNumMax),stat=istat)
        if (istat/=0) call Driver_abort("Cannot allocate pbIgnTime in paraburn init")
     allocate(pbIgnX(pbIgnNumMax),stat=istat)
        if (istat/=0) call Driver_abort("Cannot allocate pbIgnX in paraburn init")
     allocate(pbIgnY(pbIgnNumMax),stat=istat)
        if (istat/=0) call Driver_abort("Cannot allocate pbIgnY in paraburn init")
     allocate(pbIgnZ(pbIgnNumMax),stat=istat)
        if (istat/=0) call Driver_abort("Cannot allocate pbIgnZ in paraburn init")

     if (restart) then
        
        call Logfile_stampMessage('[Burn/parametric init] Reading previous detonation ignition points (if any) from file')
        open(unit=infile_unit,file=bn_detIgnFileName,status='unknown',iostat=istat)
        if (istat /= 0) call Driver_abort("Unable to open detonation ignition points file")

        istat = 0
        do i = 1, pbIgnNumMax

           read(infile_unit,*,IOSTAT=istat) pbIgnTime(i), pbIgnX(i),  &
                                               pbIgnY(i), pbIgnZ(i)
           if (istat/=0) exit

        enddo

        if (istat > 0) then
           call Driver_abort("Unable to read detonation ignition points file")
        else if (istat < 0) then !EOF reached
           pbIgnNum = i - 1
        else
           call Driver_abort("pbIgnNumMax is too small to read in previous detonation ignition points")
        endif

        close(unit=infile_unit)
     endif

  else

     pbIgnNum = 0

  endif

  ! use built-in physical constants
  call PhysicalConstants_get("electron mass", m_e)
  call PhysicalConstants_get("Avogadro", N_A)
  call PhysicalConstants_get("speed of light", c_l)
  call PhysicalConstants_get("proton mass", m_p)
  m_n = 1.67492716e-24

  ! read nuclear info from file
  ! form cribbed from SimulationComposation
  open(unit=infile_unit,file="SpeciesList.txt")
  nfound=0
  nread=0
  do while((nfound<6).and.(nread<=spec_num))
     nread=nread+1
     read(infile_unit,*)isotopeName,zbar,abar,spinfactor,bindEnergy
     if (trim(isotopeName) .eq. 'he4') then
        yi_he4 = 1.0/abar
        ye_he4 = zbar/abar
        q_he4 = bindEnergy/abar
        nfound = nfound + 1
     else if (trim(isotopeName) .eq. 'c12') then
        yi_c12 = 1.0/abar
        ye_c12 = zbar/abar
        q_c12  = bindEnergy/abar
        nfound = nfound + 1
     else if (trim(isotopeName) .eq. 'o16') then
        yi_o16 = 1.0/abar
        ye_o16 = zbar/abar
        q_o16  = bindEnergy/abar
        nfound = nfound + 1
     else if (trim(isotopeName) .eq. 'ne22') then
        yi_ne22 = 1.0/abar
        ye_ne22 = zbar/abar
        q_ne22  = bindEnergy/abar
        nfound = nfound + 1
     else if (trim(isotopeName) .eq. 'mg24') then
        yi_mg24 = 1.0/abar
        ye_mg24 = zbar/abar
        q_mg24  = bindEnergy/abar
        nfound = nfound + 1
     else if (trim(isotopeName) .eq. 'si28') then
        yi_si28 = 1.0/abar
        ye_si28 = zbar/abar
        q_si28  = bindEnergy/abar
        nfound = nfound + 1
     end if
  end do
  close(unit=infile_unit)

end subroutine Burn_init

