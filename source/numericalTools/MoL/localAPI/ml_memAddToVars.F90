!!****f* source/numericalTools/MoL/localAPI/ml_memAddToVars
!!
!!  NAME
!!
!!      ml_memAddToVars
!!
!!  SYNOPSIS
!!
!!      call ml_memAddToVars(integer, intent(in) :: dst
!!                           real,    intent(in) :: dstFac
!!                           integer, intent(in) :: nsrcs
!!                           integer, intent(in) :: srcs(nsrcs)
!!                           real,    intent(in) :: facs(nsrcs))
!!
!!      call ml_memAddToVars(integer, intent(in) :: dst
!!                           real,    intent(in) :: dstFac
!!                           real,    intent(in) :: val)
!!
!!      call ml_memAddToVars(integer, intent(in) :: dst
!!                           real,    intent(in) :: dstFac
!!                           integer, intent(in) :: src1, src2, ...
!!                           real,    intent(in) :: fac1, fac2, ...)
!!
!!  DESCRIPTION
!!
!!      Perform a linear combination of source terms into the specified destination
!!      for all variables evolved by MoL.  The destination can either be the evolved
!!      variables in UNK or their corresponding locations in one of the MoL-specific
!!      scratch memory locations.  The source terms will only-ever be taken from
!!      MoL-specific scratch memory locations.  The linear combinations will take one
!!      of the following forms:
!!
!!          dst = dstFac*dst + val
!!
!!              or
!!
!!          dst = dstFac*dst + fac1*src1 + ... + facN*srcN
!!
!!      Valid locations include (defined in MoL.h):
!!          - MOL_EVOLVED : Evolved variables in UNK
!!          - MOL_INITIAL : Copy of the evolved variables at the start of a timestep
!!          - MOL_RHS     : The currently-being-calculated RHS terms
!!          - other       : Each integrator may specify some additional number of
!!                          of scratch-memory for intermediate stages/RHS terms
!!
!!  ARGUMENTS
!!
!!      dst    : Index of the destination location to store the linear combination
!!      dstFac : Scaling factor for the destination - set this to zero to overwrite
!!               the existing value
!!      val    : Scalar value to set all variables to in dst
!!      nsrcs  : Number of source terms for the general N-src implementation
!!      srcs   : Array of source index locations in MoL scratch memory
!!      facs   : Array of scaling factors for each source term
!!      src1
!!       |     : Index of each source term location to use
!!      srcN
!!      fac1
!!       |     : Scaling factor for each source term location to use
!!      facN
!!
!!  NOTES
!!
!!      For optimal memory use, there are specific implementations for up to twelve
!!      source terms - this will need to increase if any methods beyond 4th-order
!!      IMEX(-MRI) are added in the future.  For now, twelve source terms are
!!      sufficient for the 4rd-order IMEX-MRI-GARK multi-rate integrator
!!***
subroutine ml_memAddToVarsN(dst, dstFac, nsrcs, srcs, facs)
   implicit none

   integer, intent(in) :: dst
   real, intent(in) :: dstFac
   integer, intent(in) :: nsrcs
   integer, intent(in) :: srcs(nsrcs)
   real, intent(in) :: facs(nsrcs)

   return
end subroutine ml_memAddToVarsN

subroutine ml_memAddToVars0(dst, dstFac, val)
   implicit none

   integer, intent(in) :: dst
   real, intent(in) :: dstFac
   real, intent(in) :: val

   return
end subroutine ml_memAddToVars0

subroutine ml_memAddToVars1(dst, dstFac, &
                            src1, &
                            fac1)
   implicit none

   integer, intent(in) :: dst
   real, intent(in) :: dstFac
   integer, intent(in) :: src1
   real, intent(in) :: fac1

   return
end subroutine ml_memAddToVars1

subroutine ml_memAddToVars2(dst, dstFac, &
                            src1, src2, &
                            fac1, fac2)
   implicit none

   integer, intent(in) :: dst
   real, intent(in) :: dstFac
   integer, intent(in) :: src1, src2
   real, intent(in) :: fac1, fac2

   return
end subroutine ml_memAddToVars2

subroutine ml_memAddToVars3(dst, dstFac, &
                            src1, src2, src3, &
                            fac1, fac2, fac3)
   implicit none

   integer, intent(in) :: dst
   real, intent(in) :: dstFac
   integer, intent(in) :: src1, src2, src3
   real, intent(in) :: fac1, fac2, fac3

   return
end subroutine ml_memAddToVars3

subroutine ml_memAddToVars4(dst, dstFac, &
                            src1, src2, src3, src4, &
                            fac1, fac2, fac3, fac4)
   implicit none

   integer, intent(in) :: dst
   real, intent(in) :: dstFac
   integer, intent(in) :: src1, src2, src3, src4
   real, intent(in) :: fac1, fac2, fac3, fac4

   return
end subroutine ml_memAddToVars4

subroutine ml_memAddToVars5(dst, dstFac, &
                            src1, src2, src3, src4, src5, &
                            fac1, fac2, fac3, fac4, fac5)
   implicit none

   integer, intent(in) :: dst
   real, intent(in) :: dstFac
   integer, intent(in) :: src1, src2, src3, src4, src5
   real, intent(in) :: fac1, fac2, fac3, fac4, fac5

   return
end subroutine ml_memAddToVars5

subroutine ml_memAddToVars6(dst, dstFac, &
                            src1, src2, src3, src4, src5, src6, &
                            fac1, fac2, fac3, fac4, fac5, fac6)
   implicit none

   integer, intent(in) :: dst
   real, intent(in) :: dstFac
   integer, intent(in) :: src1, src2, src3, src4, src5, src6
   real, intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6

   return
end subroutine ml_memAddToVars6

subroutine ml_memAddToVars7(dst, dstFac, &
                            src1, src2, src3, src4, src5, src6, &
                            src7, &
                            fac1, fac2, fac3, fac4, fac5, fac6, &
                            fac7)
   implicit none

   integer, intent(in) :: dst
   real, intent(in) :: dstFac
   integer, intent(in) :: src1, src2, src3, src4, src5, src6, &
                          src7
   real, intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, &
                       fac7

   return
end subroutine ml_memAddToVars7

subroutine ml_memAddToVars8(dst, dstFac, &
                            src1, src2, src3, src4, src5, src6, &
                            src7, src8, &
                            fac1, fac2, fac3, fac4, fac5, fac6, &
                            fac7, fac8)
   implicit none

   integer, intent(in) :: dst
   real, intent(in) :: dstFac
   integer, intent(in) :: src1, src2, src3, src4, src5, src6, &
                          src7, src8
   real, intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, &
                       fac7, fac8

   return
end subroutine ml_memAddToVars8

subroutine ml_memAddToVars9(dst, dstFac, &
                            src1, src2, src3, src4, src5, src6, &
                            src7, src8, src9, &
                            fac1, fac2, fac3, fac4, fac5, fac6, &
                            fac7, fac8, fac9)
   implicit none

   integer, intent(in) :: dst
   real, intent(in) :: dstFac
   integer, intent(in) :: src1, src2, src3, src4, src5, src6, &
                          src7, src8, src9
   real, intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, &
                       fac7, fac8, fac9

   return
end subroutine ml_memAddToVars9

subroutine ml_memAddToVars10(dst, dstFac, &
                             src1, src2, src3, src4, src5, src6, &
                             src7, src8, src9, src10, &
                             fac1, fac2, fac3, fac4, fac5, fac6, &
                             fac7, fac8, fac9, fac10)
   implicit none

   integer, intent(in) :: dst
   real, intent(in) :: dstFac
   integer, intent(in) :: src1, src2, src3, src4, src5, src6, &
                          src7, src8, src9, src10
   real, intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, &
                       fac7, fac8, fac9, fac10

   return
end subroutine ml_memAddToVars10

subroutine ml_memAddToVars11(dst, dstFac, &
                             src1, src2, src3, src4, src5, src6, &
                             src7, src8, src9, src10, src11, &
                             fac1, fac2, fac3, fac4, fac5, fac6, &
                             fac7, fac8, fac9, fac10, fac11)
   implicit none

   integer, intent(in) :: dst
   real, intent(in) :: dstFac
   integer, intent(in) :: src1, src2, src3, src4, src5, src6, &
                          src7, src8, src9, src10, src11
   real, intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, &
                       fac7, fac8, fac9, fac10, fac11

   return
end subroutine ml_memAddToVars11

subroutine ml_memAddToVars12(dst, dstFac, &
                             src1, src2, src3, src4, src5, src6, &
                             src7, src8, src9, src10, src11, src12, &
                             fac1, fac2, fac3, fac4, fac5, fac6, &
                             fac7, fac8, fac9, fac10, fac11, fac12)
   implicit none

   integer, intent(in) :: dst
   real, intent(in) :: dstFac
   integer, intent(in) :: src1, src2, src3, src4, src5, src6, &
                          src7, src8, src9, src10, src11, src12
   real, intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, &
                       fac7, fac8, fac9, fac10, fac11, fac12

   return
end subroutine ml_memAddToVars12
