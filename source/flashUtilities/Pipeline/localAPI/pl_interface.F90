!!****ih* source/flashUtilities/Pipeline/localAPI/pl_interface
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
!!  pl_interface
!!
!! SYNOPSIS
!!
!!   use pl_interface
!!
!! DESCRIPTION
!!
!!  This is the header file for the Pipeline unit that defines its
!!  private interfaces.
!!
!!***

Module pl_interface

  interface
     subroutine pl_localPipelineSetup () 
     end subroutine pl_localPipelineSetup
  end interface

  interface
     subroutine pl_localPostRecvMsg (channel)
       integer, intent (in) :: channel
     end subroutine pl_localPostRecvMsg
  end interface

  interface
     subroutine pl_localPostSendMsg (channel)
       integer, intent (in) :: channel
     end subroutine pl_localPostSendMsg
  end interface

  interface
     subroutine pl_localSaveRecvItems (channel, isSaved)
       integer, intent (in)  :: channel
       logical, intent (out) :: isSaved
     end subroutine pl_localSaveRecvItems
  end interface

  interface
     subroutine pl_printGlobalStatusVector (fileUnit)
       integer, intent (in) :: fileUnit
     end subroutine pl_printGlobalStatusVector
  end interface

end Module pl_interface
