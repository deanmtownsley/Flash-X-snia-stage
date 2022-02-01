!!****if* source/Simulation/SimulationMain/YahilLattimerCollapse/Simulation_initBlock
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_initBlock(real,pointer :: solnData(:,:,:,:),
!!                            Grid_tile_t(IN)  :: tileDesc  )
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the YahilLattimerCollapse
!!  explosion problem.
!!
!!  References:  Yahil, A. (1983). Self-similar stellar collapse. 
!!               The Astrophysical Journal, 265, 1047
!!
!! ARGUMENTS
!!
!!  solnData  -        pointer to solution data
!!  tileDesc -         describes the tile or block to initialize
!!
!!
!! PARAMETERS
!!
!!
!!
!!***

!!REORDER(4): solnData

subroutine Simulation_initBlock(solnData,tileDesc)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getCellCoords, &
                             Grid_getCellVolumes, &
                             Grid_subcellGeometry, &
                             Grid_getDeltas
  use Grid_tile, ONLY : Grid_tile_t 
  use Driver_interface, ONLY : Driver_abort
  use ut_interpolationInterface
  use model_interp_module
 
  implicit none

#include "constants.h"
#include "Simulation.h"
  
  real,              pointer    :: solnData(:,:,:,:)
  type(Grid_tile_t), intent(in) :: tileDesc
  
  integer  ::  i, j, k, n, jLo, jHi
  integer  ::  ii, jj, kk
  real     ::  drProf
  real,allocatable,dimension(:) :: rProf, vProf, rhoProf, pProf
  real,allocatable,dimension(:) :: eProf, mProf
  real     ::  distInv, xDist, yDist, zDist
  real     ::  sumRho, sumP, sumVX, sumVY, sumVZ
  real     ::  vel, diagonal
  real     ::  xx, dxx, yy, dyy, zz, dzz, frac
  real     ::  vx, vy, vz, p, rho, e, ek, eint
  real     ::  dist
  real     ::  vSub, rhoSub, pSub, errIgnored

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  real,allocatable :: cellVolumes(:, :, :)
  integer,dimension(LOW:HIGH,MDIM) :: tileLimits
  integer,dimension(LOW:HIGH,MDIM) :: grownTileLimits
  integer,dimension(MDIM) :: axis

  real,allocatable :: dvSub(:,:)
  real     :: dvc, quotinv

  real :: deltas(1:MDIM)

  if (sim_tinitial > 0.0) &
     call Driver_abort('EORROR: not able to take tinitial > 0')

  if (sim_useProfileFromFile) then
   
     allocate(rProf(sim_initReso))
     allocate(vProf(sim_initReso))
     allocate(rhoProf(sim_initReso))
     allocate(pProf(sim_initReso))
     allocate(mProf(sim_initReso))
     allocate(eProf(sim_initReso))

  !
  !  Construct the radial samples needed for the initialization.
  !
     diagonal = (sim_xMax-sim_xMin)**2
     diagonal = diagonal + K2D*(sim_yMax-sim_yMin)**2
     diagonal = diagonal + K3D*(sim_zMax-sim_zMin)**2
     diagonal = sqrt(diagonal)
  
     drProf = diagonal / (sim_initReso-1)

     do i = 1, sim_initReso
        rProf(i)   = (i-1) * drProf
     enddo

     CALL interp1d_linear(sim_rProf, sim_rhoProf, rProf, rhoProf)
     CALL interp1d_linear(sim_rProf, sim_pProf,   rProf, pProf)
     CALL interp1d_linear(sim_rProf, sim_vProf,   rProf, vProf  )
     CALL interp1d_linear(sim_rProf, sim_mProf,   rProf, mProf  )
     CALL interp1d_linear(sim_rProf, sim_eProf,   rProf, eProf  )

  end if !useProfileFromFile

  ! get the coordinate information for the current block
  tileLimits = tileDesc%limits
  grownTileLimits = tileDesc%grownLimits
    
  call Grid_getDeltas(tileDesc%level, deltas)

  ! Find a real difference between z's if problem is >= 3D
  if (NDIM > 2) then
     dzz = deltas(KAXIS)
  ! Otherwise this problem is <= 2D, so dzz is meaningless
  else
     dzz = 0.0
  endif

  ! Find a real difference between y's if problem is >= 2D
  if (NDIM > 1) then
     dyy = deltas(JAXIS)
  ! Otherwise this problem is <= 1D, so dyy is meaningless
  else
    dyy = 0.0
  endif

  dxx = deltas(IAXIS)

  allocate(xCoord(grownTileLimits(LOW, IAXIS):grownTileLimits(HIGH, IAXIS))); xCoord = 0.0
  allocate(yCoord(grownTileLimits(LOW, JAXIS):grownTileLimits(HIGH, JAXIS))); yCoord = 0.0
  allocate(zCoord(grownTileLimits(LOW, KAXIS):grownTileLimits(HIGH, KAXIS))); zCoord = 0.0
  allocate(cellVolumes(grownTileLimits(LOW, IAXIS):grownTileLimits(HIGH, IAXIS), &
                       grownTileLimits(LOW, JAXIS):grownTileLimits(HIGH, JAXIS), &
                       grownTileLimits(LOW, KAXIS):grownTileLimits(HIGH, KAXIS)))

  call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                          grownTileLimits(LOW,  :), &
                          grownTileLimits(HIGH, :), &
                          xCoord)
#if NDIM >= 2
  call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, &
                          grownTileLimits(LOW,  :), &
                          grownTileLimits(HIGH, :), &
                          yCoord)
#endif
#if NDIM == 3
  call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, &
                          grownTileLimits(LOW,  :), &
                          grownTileLimits(HIGH, :), &
                          zCoord)
#endif
  call Grid_getCellVolumes(tileDesc%level, &
                           lbound(cellVolumes), ubound(cellVolumes), cellVolumes)

  !
  !     For each cell
  !  

  allocate(dvSub(0:sim_nSubZones-1,0:(sim_nSubZones-1)*K2D))

  do k = grownTileLimits(LOW,KAXIS), grownTileLimits(HIGH,KAXIS)
     zz = zCoord(k)

     do j = grownTileLimits(LOW, JAXIS), grownTileLimits(HIGH, JAXIS)
        yy = yCoord(j)
        
        do i = grownTileLimits(LOW,IAXIS), grownTileLimits(HIGH, IAXIS)
           xx = xCoord(i)

           dvc = cellVolumes(i, j, k)
           call Grid_subcellGeometry(   sim_nSubZones, &
                                     1+(sim_nSubZones-1)*K2D, &
                                     1+(sim_nSubZones-1)*K3D, &
                                     dvc, dvSub, &
                                     xCoord(i)-0.5*dxx, xCoord(i)+0.5*dxx)

           sumRho = 0.
           sumP   = 0.
           sumVX  = 0.
           sumVY  = 0.
           sumVZ  = 0.
           
           !
           !       Break the cell into sim_nSubZones^NDIM sub-zones, and look up the
           !       appropriate quantities along the 1d profile for that subzone.  
           !
           !       Have the final values for the zone be equal to the average of
           !       the subzone values.
           ! 

           do kk = 0, (sim_nSubZones-1)*K3D
              zz    = zCoord(k) + ((real(kk)+0.5)*sim_inSubzones-.5)*dzz 
              zDist = (zz - sim_zCenter) * K3D
              
              do jj = 0, (sim_nSubZones-1)*K2D
                 yy    = yCoord(j) + ((real(jj)+0.5)*sim_inSubzones-.5)*dyy
                 yDist = (yy - sim_yCenter) * K2D
                 
                 do ii = 0, (sim_nSubZones-1)
                    xx    = xCoord(i) + ((real(ii)+0.5)*sim_inSubzones-.5)*dxx
                    xDist = xx - sim_xCenter
                    
                    dist    = sqrt( xDist**2 + yDist**2 + zDist**2 )
                    distInv = 1. / max( dist, 1.E-10 )
                    call sim_find (rProf, sim_initReso, dist, jLo)
                    !
                    !  a point at `dist' is frac-way between jLo and jHi.   We do a
                    !  linear interpolation of the quantities at jLo and jHi and sum those.
                    ! 
                    if (jLo .eq. 0) then
                       jLo = 1
                       jHi = 1
                       frac = 0.
                    else if (jLo .eq. sim_initReso) then
                       jLo = sim_initReso
                       jHi = sim_initReso
                       frac = 0.
                    else
                       jHi = jLo + 1
                       frac = (dist - rProf(jLo)) / & 
                         (rProf(jHi)-rProf(jLo))
                    endif

                    pSub   =  pProf(jLo) + frac*(pProf(jHi)  - pProf(jLo))

                    rhoSub =  rhoProf(jLo) + frac*(rhoProf(jHi)- rhoProf(jLo))

                    vSub   = vProf(jLo) + frac*(vProf(jHi)  - vProf(jLo))

                    ! 
                    !   Now total these quantities.   Note that  v is a radial velocity; 
                    !   we multiply by the tangents of the appropriate angles to get
                    !   the projections in the x, y and z directions.
                    !
                    sumP = sumP + pSub * dvSub(ii,jj)
                    
                    sumRho = sumRho + rhoSub * dvSub(ii,jj)
                    
                    vel = vSub * dvSub(ii,jj)
                    
                    sumVX  = sumVX  + vel*xDist*distInv
                    sumVY  = sumVY  + vel*yDist*distInv
                    sumVZ  = sumVZ  + vel*zDist*distInv
                    
                 enddo ! ii
              enddo ! jj
           enddo ! kk
           
!!$           quotinv = sim_inszd
           quotinv = 1.0 / dvc
           rho = sumRho * quotinv
           p   = sumP   * quotinv
           vx  = sumVX  * quotinv
           vy  = sumVY  * quotinv
           vz  = sumVZ  * quotinv
           ek  = 0.5*(vx*vx + vy*vy + vz*vz)
           !
           !  assume gamma-law equation of state
           !
           e   = p/(sim_gamma-1.)
           eint= e/rho
           e   = e/rho + ek
           
           axis(IAXIS)=i
           axis(JAXIS)=j
           axis(KAXIS)=k


           if (NSPECIES > 0) then
              solnData(SPECIES_BEGIN,i,j,k)=1.0-(NSPECIES-1)*sim_smallX
              solnData(SPECIES_BEGIN+1:SPECIES_END,i,j,k)=sim_smallX
           end if
           solnData(DENS_VAR,i,j,k)=rho
           solnData(PRES_VAR,i,j,k)=p
           solnData(ENER_VAR,i,j,k)=e
#ifdef EINT_VAR
           solnData(EINT_VAR,i,j,k)=eint
#endif
           solnData(GAME_VAR,i,j,k)=sim_gamma
           solnData(GAMC_VAR,i,j,k)=sim_gamma
           solnData(VELX_VAR,i,j,k)=vx
           solnData(VELY_VAR,i,j,k)=vy
           solnData(VELZ_VAR,i,j,k)=vz
           solnData(TEMP_VAR,i,j,k)=sim_smallT
#ifdef BDRY_VAR
           solnData(BDRY_VAR,i,j,k)=    -1.0
#endif
        enddo
     enddo
  enddo

  deallocate(dvSub)

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  deallocate(cellVolumes)

  deallocate(rProf)
  deallocate(vProf)
  deallocate(rhoProf)
  deallocate(pProf)
  deallocate(mProf)
  deallocate(eProf)

  return
end subroutine Simulation_initBlock

!******************************************************************************

!  Routine:     sim_find()

!  Description: Given a monotonically increasing table x(N) and a test value
!               x0, return the index i of the largest table value less than
!               or equal to x0 (or 0 if x0 < x(1)).  Use binary search.

subroutine sim_find (x, N, x0, i)

  implicit none

! Arguments, LBR guessed intent on these
  integer, intent(IN) :: N
  integer, intent(OUT):: i
  real, intent(IN)    :: x(N), x0

! local variables
  integer  il, ir, im

  if (x0 .lt. x(1)) then

     i = 0

  elseif (x0 .gt. x(N)) then

     i = N

  else

     il = 1
     ir = N
10   if (ir .eq. il+1) goto 20
     im = (il + ir) / 2
     if (x(im) .gt. x0) then
        ir = im
     else
        il = im
     endif
     goto 10
20   i = il

  endif

  return
end subroutine sim_find
