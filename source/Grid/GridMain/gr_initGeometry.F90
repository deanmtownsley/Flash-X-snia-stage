!!****if* source/Grid/GridMain/gr_initGeometry
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
!!  gr_initGeometry
!!
!!
!! SYNOPSIS
!!
!!  call gr_initGeometry()
!!
!!
!! DESCRIPTION
!!
!!  Perform Grid data initializations that are related to the
!!  geometry of the simulation. 
!!
!!  Uses the following data, which must have been
!!  established before the call:
!!    gr_domainBC(LOW:HIGH,IAXIS:KAXIS)
!!
!!  Initializes the following data, which will be available
!!  after the call:
!!    gr_geometry (one of CARTESIAN, POLAR, SPHERICAL, or CYLINDRICAL)
!!    gr_geometryOverride
!!    gr_dirGeom
!!    gr_dirIsAngular
!!    gr_imin, gr_jmin, gr_kmin      (UG or Paramesh Grid only)
!!    gr_imax, gr_jmax, gr_kmax      (UG or Paramesh Grid only)
!!    dr_globalDomain(LOW:HIGH,IAXIS:KAXIS)
!!
!! NOTES
!!
!!  This subroutine is normally called from Grid_init.
!!  This code was moved into a separate file because it is
!!  essentially independent of the choice of Grid implementation.
!!
!!  Here are some initialization things that are NOT done in this
!!  subroutine, and that should therefore generally occur in later steps:
!!    - maximum refinement levels
!!    - per-level quantities, for example, deltas (cell spacings) for each level
!!
!!  CARTESIAN, POLAR, SPHERICAL, and CYLINDRICAL are defined
!!  in constants.h.
!!
!! SIDE EFFECTS
!!
!!  Loads the following runtime parameter into Grid_data module variables:
!!    geometry   into  gr_str_geometry (as string) and gr_geometry (as integer)
!!    geometryOverride  into  gr_geometryOverride
!!    xmin,xmax,ymin,ymax,zmin,zmax
!!                    - into elements of dr_globalDomain
!!                    - into gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax
!!                      (additionally), if the Grid is UG or PARAMESH.
!!
!!  On return, the components of gr_dirGeom will be set appropriately for
!!  the simulation's geometry as indicated by gr_geometry.
!!
!!  For angle coordinates, values will be scaled by multiplying with pi/180.
!!  This may affect ymin, ymax and/or zmin, zmax, depending on gr_geometry.
!!  For example, if the JAXIS grid direction stands for an angle (as is the case
!!  in POLAR and SPHERICAL coordinates),  ymin and ymax given in degrees
!!  in a flash.par file are converted here to radians.
!!
!!
!! ARGUMENTS
!!
!!  none
!!
!! SEE ALSO
!!
!!  Grid_init
!!  constants.h
!!
!!***

#include "Simulation.h"

subroutine gr_initGeometry()

  use Driver_interface, ONLY : Driver_abort
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
       RuntimeParameters_mapStrToInt
  use RuntimeParameters_interface, ONLY : RuntimeParameters_setReal
  use Logfile_interface, ONLY : Logfile_stampMessage,Logfile_stamp

  use Grid_data, ONLY : gr_geometry, gr_geometryOverride, &
                        gr_str_geometry, &
                        gr_dirGeom, gr_dirIsAngular, gr_domainBC, &
                        gr_globalDomain
#if defined(FLASH_GRID_UG) || defined(FLASH_GRID_PARAMESH)
  use Grid_data, ONLY : gr_imin, gr_imax, gr_jmin, gr_jmax, gr_kmin, gr_kmax
#endif

  implicit none

#include "constants.h"

  integer :: idir

  real    :: imin, imax, jmin, jmax, kmin, kmax

  call RuntimeParameters_get("geometry",gr_str_geometry)
  call RuntimeParameters_mapStrToInt(gr_str_geometry, gr_geometry)

  call RuntimeParameters_get("geometryOverride", gr_geometryOverride)
  call RuntimeParameters_get('xmin', imin)
  call RuntimeParameters_get('xmax', imax)
  call RuntimeParameters_get('ymin', jmin)
  call RuntimeParameters_get('ymax', jmax)
  call RuntimeParameters_get('zmin', kmin)
  call RuntimeParameters_get('zmax', kmax)

  do idir=1,NDIM
     call Logfile_stamp(idir,'[gr_initGeometry] checking BCs for idir')
     if (gr_domainBC(LOW,idir) .GE. NOT_BOUNDARY) then
        call Logfile_stamp('LOW boundary condition type not recognized','[gr_initGeometry]')
     end if
     if (gr_domainBC(HIGH,idir) .GE. NOT_BOUNDARY) then
        call Logfile_stamp('HIGH boundary condition type not recognized','[gr_initGeometry]')
     end if
     if (gr_domainBC(LOW,idir) .GE. NOT_BOUNDARY .OR. gr_domainBC(HIGH,idir) .GE. NOT_BOUNDARY) then
        call Driver_abort('A boundary condition type not recognized. &
             & Maybe [xyz][lr]_boundary_type runtime parameter specification was invalid.')
     end if
  end do

#ifdef GRID_GEOM_CARTESIAN
  if (gr_geometry .ne. CARTESIAN) then
     print *,'WARNING The geometry runtime parameter is different from&
          & CARTESIAN geometry specified at setup time.'
     call Logfile_stampMessage('WARNING The geometry runtime parameter is different from&
          & CARTESIAN geometry specified at setup time. FLASH will honor the runtime parameter.')
  end if
#endif
#ifdef GRID_GEOM_POLAR
  if (gr_geometry .ne. POLAR) then
     print *,'WARNING The geometry runtime parameter is different from&
          & POLAR geometry specified at setup time.'
     call Logfile_stampMessage('WARNING The geometry runtime parameter is different from&
          & POLAR geometry specified at setup time. FLASH will honor the runtime parameter.')
  end if
#endif
#ifdef GRID_GEOM_CYLINDRICAL
  if (gr_geometry .ne. CYLINDRICAL) then
     print *,'WARNING The geometry runtime parameter is different from&
          & CYLINDRICAL geometry specified at setup time.'
     call Logfile_stampMessage('WARNING The geometry runtime parameter is different from&
          & CYLINDRICAL geometry specified at setup time. FLASH will honor the runtime parameter.')
  end if
#endif
#ifdef GRID_GEOM_SPHERICAL
  if (gr_geometry .ne. SPHERICAL) then
     print *,'WARNING The geometry runtime parameter is different from&
          & SPHERICAL geometry specified at setup time.'
     call Logfile_stampMessage('WARNING The geometry runtime parameter is different from&
          & SPHERICAL geometry specified at setup time. FLASH will honor the runtime parameter.')
  end if
#endif


  gr_dirIsAngular = .FALSE.

  if (gr_geometry == CARTESIAN)then
     gr_dirGeom(IAXIS) = XYZ
     gr_dirGeom(JAXIS) = XYZ
     gr_dirGeom(KAXIS) = XYZ
  elseif(gr_geometry == POLAR)then
     gr_dirGeom(IAXIS) = RAD_CYL
     gr_dirGeom(JAXIS) = PHI_CYL
     gr_dirGeom(KAXIS) = XYZ
     gr_dirIsAngular(JAXIS) = .TRUE.
  elseif(gr_geometry == CYLINDRICAL) then
     gr_dirGeom(IAXIS) = RAD_CYL
     gr_dirGeom(JAXIS) = XYZ
     gr_dirGeom(KAXIS) = PHI_CYL
     gr_dirIsAngular(KAXIS) = .TRUE.
  elseif(gr_geometry == SPHERICAL) then
     gr_dirGeom(IAXIS) = RAD_SPH
     gr_dirGeom(JAXIS) = THETA
     gr_dirGeom(KAXIS) = PHI_SPH
     gr_dirIsAngular(JAXIS) = .TRUE.
     gr_dirIsAngular(KAXIS) = .TRUE.
  else
     call Driver_abort("[Grid_init] unsupported geometry ")
  end if


#ifdef DEBUG
!! DEV: The following was taken from init_mesh.F90 FLASH2, but
!! the restriction seems unnecessary. -KW
! finally make sure that the geometry is valid.  In particular, we only
! support (or plan to support) those geometries listed in the header of
! grid.F90.  

  if ( (gr_geometry == CYLINDRICAL .AND. NDIM == 1) .OR. &
       (gr_geometry == POLAR .AND. NDIM == 3) ) then
     
     print *, "ERROR: geometry invalid"
     call Driver_abort("ERROR: geometry invalid")
  endif
#endif

! If we are dealing with angular coordinates, the user specified the extrema
! in degrees.  Here we multiply by pi/180, and restore the extrema.  This
! way, when gr_createDomain etc. are called, the blocks will be created with the
! proper dimensions.  This ensures that all coordinate values and coordinate
! differences for this direction will be in radians.
  
  if (gr_geometry /= CARTESIAN) then

! The radial coordinate must always be >= 0.0.  If it is 0.0, then we should
! have a reflecting boundary there.

     if (imin < 0.0) then
        if (.NOT. gr_geometryOverride) &
             call Driver_abort("ERROR: radial coordinate cannot be < 0.0")
     endif

     if (imin == 0.0 .AND. &
     (gr_domainBC(LOW,IAXIS) /= REFLECTING .AND. gr_domainBC(LOW,IAXIS) /= AXISYMMETRIC)) then
        !! FUTURE: We could have some special treatment for a boundary at r=0.0,
        !! like using the singular_line code provided in Paramesh3 ff. - KW
        if (.NOT. gr_geometryOverride) &
             call Driver_abort("ERROR: reflecting or axisymmetric boundary required at x = 0.0 for radial coords")
     endif

  endif

  if (gr_dirIsAngular(JAXIS)) then

! Make sure the range is valid.
     if (gr_geometry == SPHERICAL) then

#if NDIM > 1
! y is the spherical theta coordinate.  It ranges from 0 to pi.  Since we
! are specifying the range in degrees, make sure that it is not > 180.
        if (jmin > 180.0 .OR. jmax > 180.0) then
           print *, 'ERROR: the theta coordinate in spherical geometry cannot be > pi.'
           print *, '       Check ymin and ymax.  The angles are assumed to be specified'
           print *, '       in degrees on input.'
           call Driver_abort("ERROR: theta coordinate range invalid")
        endif
#endif

     elseif (gr_geometry == POLAR) then

#if NDIM > 1
! In polar coordinates, y is the phi coordinate, which can range from 0
! to 2 pi.  
        if (jmin > 360.0 .OR. jmax > 360.0) then
           print *, 'ERROR: the phi coordinate in polar geometry cannot be > 2 pi.'
           print *, '       Check ymin and ymax.  The angles are assumed to be specified'
           print *, '       in degrees on input.'
           call Driver_abort("ERROR: phi coordinate range invalid")
        endif
#endif

     else
        
        call Driver_abort("ERROR: y cannot be an angular coordinate in this geometry")
        
     endif
        
     jmin = jmin*PI/180
     jmax = jmax*PI/180

     ! In case some other unit gets the coordinate range runtime parameters AFTER Grid_init
     ! is called, it will get the scaled values expressing coordinates in radians.
     call RuntimeParameters_setReal("ymin", jmin)
     call RuntimeParameters_setReal("ymax", jmax)

  endif

  if (gr_dirIsAngular(KAXIS)) then

! Make sure the range is valid.
     if (gr_geometry == SPHERICAL .OR. gr_geometry == CYLINDRICAL) then

#if NDIM > 2
! z is the phi coordinate in both spherical and cylindrical coords.
        if (kmin > 360.0 .OR. kmax > 360.0) then
           print *, 'ERROR: the phi coordinate in the current geometry cannot be > 2 pi.'
           print *, '       Check zmin and zmax.  The angles are assumed to be specified'
           print *, '       in degrees on input.'
           call Driver_abort("ERROR: phi coordinate range invalid")
        endif
#endif

     else
        
        call Driver_abort("ERROR: z cannot be an angular coordinate in this geometry")
        
     endif

     kmin = kmin*PI/180
     kmax = kmax*PI/180

     ! In case some other unit gets the coordinate range runtime parameters AFTER Grid_init
     ! is called, it will get the scaled values expressing coordinates in radians.
     call RuntimeParameters_setReal("zmin", kmin)
     call RuntimeParameters_setReal("zmax", kmax)
  endif

  ! Store computational domain limits in a convenient array.
  ! Implementations of Grid_getDomainBoundBox may use this array
  ! to return the correct domain bounds.
  gr_globalDomain(LOW,IAXIS) = imin
  gr_globalDomain(LOW,JAXIS) = jmin
  gr_globalDomain(LOW,KAXIS) = kmin
  gr_globalDomain(HIGH,IAXIS) = imax
  gr_globalDomain(HIGH,JAXIS) = jmax
  gr_globalDomain(HIGH,KAXIS) = kmax

  ! Some implementations of Grid_getDomainBoundBox have their own
  ! way to return the correct domain bounds. For UG and PARAMESH,
  ! module variables gr_imin,...,gr_kmax, are used.
#if defined(FLASH_GRID_UG) || defined(FLASH_GRID_PARAMESH)
  gr_imin = imin
  gr_jmin = jmin
  gr_kmin = kmin
  gr_imax = imax
  gr_jmax = jmax
  gr_kmax = kmax
#endif

end subroutine gr_initGeometry
