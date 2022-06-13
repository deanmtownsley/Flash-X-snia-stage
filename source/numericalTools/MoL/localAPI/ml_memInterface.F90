!!****h* source/numericalTools/MoL/localAPI/ml_memInterface
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
!!  NAME
!!
!!      ml_memInterface
!!
!!  SYNOPSIS
!!
!!      use ml_memInterface
!!
!!  DESCRIPTION
!!
!!      This is the header file for the method of lines time integration unit
!!      that defines memory storage and operations procedures
!!
!!***
module ml_memInterface

    implicit none

    !! ================================ !!
    !!  Memory allocation/deallocation  !!
    !! ================================ !!

    interface
        subroutine ml_memAlloc
        end subroutine ml_memAlloc
    end interface

    interface
        subroutine ml_memFree
        end subroutine ml_memFree
    end interface


    !! ===================================== !!
    !!  Directly set or copy to/from memory  !!
    !! ===================================== !!

    interface
        subroutine ml_memZero(dst)
            integer, intent(in) :: dst
        end subroutine ml_memZero
    end interface

    interface
        subroutine ml_memCopy(dst, src)
            integer, intent(in) :: dst, src
        end subroutine ml_memCopy
    end interface


    !! ============================== !!
    !!  Linear combination operators  !!
    !! ============================== !!

    interface ml_memAddToVars
        subroutine ml_memAddToVarsN(dst, dstFac, nsrcs, srcs, facs)
            integer, intent(in) :: dst
            real,    intent(in) :: dstFac
            integer, intent(in) :: nsrcs
            integer, intent(in) :: srcs(nsrcs)
            real,    intent(in) :: facs(nsrcs)
        end subroutine ml_memAddToVarsN

        subroutine ml_memAddToVars0(dst, dstFac, val)
            integer, intent(in) :: dst
            real,    intent(in) :: dstFac
            real,    intent(in) :: val
        end subroutine ml_memAddToVars0

        subroutine ml_memAddToVars1(dst, dstFac, &
                                    src1,        &
                                    fac1)
            integer, intent(in) :: dst
            real,    intent(in) :: dstFac
            integer, intent(in) :: src1
            real,    intent(in) :: fac1
        end subroutine ml_memAddToVars1

        subroutine ml_memAddToVars2(dst, dstFac, &
                                    src1, src2,  &
                                    fac1, fac2)
            integer, intent(in) :: dst
            real,    intent(in) :: dstFac
            integer, intent(in) :: src1, src2
            real,    intent(in) :: fac1, fac2
        end subroutine ml_memAddToVars2

        subroutine ml_memAddToVars3(dst, dstFac,      &
                                    src1, src2, src3, &
                                    fac1, fac2, fac3)
            integer, intent(in) :: dst
            real,    intent(in) :: dstFac
            integer, intent(in) :: src1, src2, src3
            real,    intent(in) :: fac1, fac2, fac3
        end subroutine ml_memAddToVars3

        subroutine ml_memAddToVars4(dst, dstFac,            &
                                    src1, src2, src3, src4, &
                                    fac1, fac2, fac3, fac4)
            integer, intent(in) :: dst
            real,    intent(in) :: dstFac
            integer, intent(in) :: src1, src2, src3, src4
            real,    intent(in) :: fac1, fac2, fac3, fac4
        end subroutine ml_memAddToVars4

        subroutine ml_memAddToVars5(dst, dstFac,                  &
                                    src1, src2, src3, src4, src5, &
                                    fac1, fac2, fac3, fac4, fac5)
            integer, intent(in) :: dst
            real,    intent(in) :: dstFac
            integer, intent(in) :: src1, src2, src3, src4, src5
            real,    intent(in) :: fac1, fac2, fac3, fac4, fac5
        end subroutine ml_memAddToVars5

        subroutine ml_memAddToVars6(dst, dstFac,                        &
                                    src1, src2, src3, src4, src5, src6, &
                                    fac1, fac2, fac3, fac4, fac5, fac6)
            integer, intent(in) :: dst
            real,    intent(in) :: dstFac
            integer, intent(in) :: src1, src2, src3, src4, src5, src6
            real,    intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6
        end subroutine ml_memAddToVars6

        subroutine ml_memAddToVars7(dst, dstFac,                        &
                                    src1, src2, src3, src4, src5, src6, &
                                    src7,                               &
                                    fac1, fac2, fac3, fac4, fac5, fac6, &
                                    fac7)
            integer, intent(in) :: dst
            real,    intent(in) :: dstFac
            integer, intent(in) :: src1, src2, src3, src4, src5, src6, &
                                   src7
            real,    intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, &
                                   fac7
        end subroutine ml_memAddToVars7

        subroutine ml_memAddToVars8(dst, dstFac,                        &
                                    src1, src2, src3, src4, src5, src6, &
                                    src7, src8,                         &
                                    fac1, fac2, fac3, fac4, fac5, fac6, &
                                    fac7, fac8)
            integer, intent(in) :: dst
            real,    intent(in) :: dstFac
            integer, intent(in) :: src1, src2, src3, src4, src5, src6, &
                                   src7, src8
            real,    intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, &
                                   fac7, fac8
        end subroutine ml_memAddToVars8

        subroutine ml_memAddToVars9(dst, dstFac,                        &
                                    src1, src2, src3, src4, src5, src6, &
                                    src7, src8, src9,                   &
                                    fac1, fac2, fac3, fac4, fac5, fac6, &
                                    fac7, fac8, fac9)
            integer, intent(in) :: dst
            real,    intent(in) :: dstFac
            integer, intent(in) :: src1, src2, src3, src4, src5, src6, &
                                   src7, src8, src9
            real,    intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, &
                                   fac7, fac8, fac9
        end subroutine ml_memAddToVars9

        subroutine ml_memAddToVars10(dst, dstFac,                        &
                                     src1, src2, src3, src4, src5, src6, &
                                     src7, src8, src9, src10,            &
                                     fac1, fac2, fac3, fac4, fac5, fac6, &
                                     fac7, fac8, fac9, fac10)
            integer, intent(in) :: dst
            real,    intent(in) :: dstFac
            integer, intent(in) :: src1, src2, src3, src4, src5, src6, &
                                   src7, src8, src9, src10
            real,    intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, &
                                   fac7, fac8, fac9, fac10
        end subroutine ml_memAddToVars10

        subroutine ml_memAddToVars11(dst, dstFac,                        &
                                     src1, src2, src3, src4, src5, src6, &
                                     src7, src8, src9, src10, src11,     &
                                     fac1, fac2, fac3, fac4, fac5, fac6, &
                                     fac7, fac8, fac9, fac10, fac11)
            integer, intent(in) :: dst
            real,    intent(in) :: dstFac
            integer, intent(in) :: src1, src2, src3, src4, src5, src6, &
                                   src7, src8, src9, src10, src11
            real,    intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6, &
                                   fac7, fac8, fac9, fac10, fac11
        end subroutine ml_memAddToVars11

        subroutine ml_memAddToVars12(dst, dstFac,                           &
                                     src1, src2, src3, src4, src5, src6,    &
                                     src7, src8, src9, src10, src11, src12, &
                                     fac1, fac2, fac3, fac4, fac5, fac6,    &
                                     fac7, fac8, fac9, fac10, fac11, fac12)
            integer, intent(in) :: dst
            real,    intent(in) :: dstFac
            integer, intent(in) :: src1, src2, src3, src4, src5, src6,    &
                                   src7, src8, src9, src10, src11, src12
            real,    intent(in) :: fac1, fac2, fac3, fac4, fac5, fac6,    &
                                   fac7, fac8, fac9, fac10, fac11, fac12
        end subroutine ml_memAddToVars12
    end interface ml_memAddToVars
end module ml_memInterface
