!!****if* source/Grid/GridMain/AMR/Paramesh4/interpolation/prolong/amr_1blk_cc_prol_gen_work_fun
!! NOTICE
!!  This file derived from PARAMESH - an adaptive mesh library.
!!  Copyright (C) 2003, 2004 United States Government as represented by the
!!  National Aeronautics and Space Administration, Goddard Space Flight
!!  Center.  All Rights Reserved.
!!  Copyright (C) 2016 The University of Chicago
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Use of the PARAMESH software is governed by the terms of the
!!  usage agreement which can be found in the file
!!  'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!!
!! NAME
!!
!!  amr_1blk_cc_prol_gen_work_fun
!!
!! SYNOPSIS
!!
!!  call amr_1blk_cc_prol_gen_work_fun(real(inout) :: recvt,
!!                                     integer(in) :: ia,
!!                                     integer(in) :: ib,
!!                                     integer(in) :: ja,
!!                                     integer(in) :: jb,
!!                                     integer(in) :: ka,
!!                                     integer(in) :: kb,
!!                                     integer(in) :: idest,
!!                                     integer(in) :: ioff,
!!                                     integer(in) :: joff,
!!                                     integer(in) :: koff,
!!                                     integer(in) :: mype,
!!                                     integer(in) :: lb,
!!                                     integer(in) :: pe_p,
!!                                     integer(in) :: lb_p,
!!                                     integer(in) :: interp)
!!
!! DESCRIPTION
!!
!!  This routine is a wrapper routine which calls the functions
!!  which prolong data for WORK. The argument idest selects the layer
!!  within WORK on which prolongation is required.
!!  The argument interp should select the actual interpolation function to be
!!  called, but is ignored here.
!!
!! ARGUMENTS
!!
!!   recvt : 
!!
!!   ia : 
!!
!!   ib : 
!!
!!   ja : 
!!
!!   jb : 
!!
!!   ka : 
!!
!!   kb : 
!!
!!   idest : 
!!
!!   ioff : 
!!
!!   joff : 
!!
!!   koff : 
!!
!!   mype : my Processor Number
!!
!!   lb : 
!!
!!   pe_p : 
!!
!!   lb_p : 
!!
!!   interp : 
!!
!! HISTORY
!!
!!  Original:     Peter MacNeice          January 2002
!!  Modified:     Klaus Weide             September 2006
!!  Modified:     Klaus Weide             September 2016 - handle interp=30..32
!!
!!***

#include "paramesh_preprocessor.fh"


      subroutine amr_1blk_cc_prol_gen_work_fun(recvt, & 
     &       ia,ib,ja,jb,ka,kb, & 
     &       idest,ioff,joff,koff,mype,lb,pe_p,lb_p,interp)




      Use paramesh_interfaces, Only :                  & 
                        amr_1blk_cc_prol_work_mg

      implicit none

      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype
      integer, intent(in) :: lb,lb_p,pe_p
      integer, intent(in) :: interp
      real,    intent(inout) :: recvt(:,:,:)

!------------------------------------

      if (interp .GE. 30 .AND. interp .LE. 32) then
! Call special version developed for MC multigrid
         call amr_1blk_cc_prol_work_mg(recvt,       &
             ia,ib,ja,jb,ka,kb,                     &
             idest,ioff,joff,koff,lb,interp-30)

      else
! Call the minimally changed subroutine for Paramesh2
         call amr_prolong_gen_work1_fun &
     &     (recvt,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff, &
     &     mype,lb)

      end if

      end subroutine amr_1blk_cc_prol_gen_work_fun
