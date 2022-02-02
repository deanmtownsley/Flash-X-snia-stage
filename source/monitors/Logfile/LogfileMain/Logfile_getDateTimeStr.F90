!!****if* source/monitors/Logfile/LogfileMain/Logfile_getDateTimeStr
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
!!  Logfile_getDateTimeStr
!!
!! SYNOPSIS
!!  call Logfile_getDateTimeStr(character(len=28)(OUT) :: dateTimeStr)
!!
!!
!! DESCRIPTION
!!  Get a string that contains the current date and time in the same
!!  format used in Logfie_stamp messages.
!!
!!  This interface is meant for use by code outside the Logfile unit
!!  that wants to generate messages in the format of the Logfile unit
!!  but circumventing  the Logfile_stamp interface.
!!
!!
!! ARGUMENTS
!!  dateTimeStr - the date and time string is returned here.
!!
!! NOTES
!!  Currently, the caller could just call current_date_time directly
!!  to the same effect as calling this subroutine, with the following
!!  difference:
!!  o The string returned by this implementation is enclosed in
!!    [ square brackets ].
!!
!!***

subroutine Logfile_getDateTimeStr(dateTimeStr)

  implicit none

  character(len=28), intent(OUT)        :: dateTimeStr

  character(len=40), save               :: localdateTimeStr

  
  call current_date_time(localDateTimeStr)

9 format('[ ',a,' ]')
  write(dateTimeStr,9) trim(localDateTimeStr)
  

end subroutine Logfile_getDateTimeStr
