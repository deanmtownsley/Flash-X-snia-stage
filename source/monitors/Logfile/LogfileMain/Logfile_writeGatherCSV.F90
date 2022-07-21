!!****if* source/monitors/Logfile/LogfileMain/Logfile_writeGatherCSV
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
!!   Logfile_writeGatherCSV
!!
!! SYNOPSIS
!!
!!   Logfile_writeGatherCSV(character(len=*):: strArr(length,dim),
!!                        integer(in)     :: length,
!!                        integer(in)     :: dim,
!!                        integer(in)     :: strLen,
!!                        integer(in)     :: numHeaders,
!!                   logical(in),optional :: reduced,
!!                   logical(in),optional :: separateFiles)
!!
!! DESCRIPTION
!!
!!   Logfile_writeGatherCSV writes the data handed to it by the Timers unit.
!!   Logfile_writeGatherCSV does not do any of the performance calculations.  Its only
!!   role is to format the data handed to it by Timers and write the data neatly to the
!!   Logfile.csv.
!!
!!   The summary implemented here  holds the 
!!   amount of indentation for each timer, and strArr(numHeadears+i:, 7) holds the indentation
!!   strArr(numHeadears+i:, 8) holds index and strArr(numHeadears+i:, 8+numProcs) holds timer o/p
!!   gathered from each proc.
!!
!!   This routine writes <logfile>.csv: 
!!    - a summary with timing information gathered from all processors is written by processor 0
!!    - no header information is written.
!!
!! ARGUMENTS
!!
!!   strArr  - array holding all the run summary information (like evolved zones, seconds
!!             in monitoring period, etc).  The first rows hold header information,
!!             the next hold timer summary information. 
!!   length  - first dimension of strArr; the number of headers + the number 
!!             of lines in the summary, + 1 for the summary column names
!!   dim     - the number of columns in the timer summary + 1 for the indentation 
!!             of the timers
!!   strLen  - length of each string entry (likely MAX_STRING_LENGTH)
!!   numHeaders - the number of name/value pairs in the header of the summary, after this timers start
!!   reduced - if present and .TRUE., generate summary of reduced timer data;
!!             otherwise generate a normal local-processor summary.
!!   separateFiles - if true, every processor writes its summary to its own file named
!!                   timer_summary_<process id>
!!
!! EXAMPLE
!!
!! initialization, 6, 0, 0.398, 0.397, 0.397
!! eos, 7, 1, 0.002, 0.002, 0.002
!! guardcell internal, 8, 1, 0.072, 0.065, 0.052
!! amr_guardcell, 9, 2, 0.069, 0.061, 0.048
!! eos gc, 10, 2, 0.002, 0.003, 0.002
!! writeCheckpoint, 11, 1, 0.022, 0.022, 0.022
!! restrictAll, 12, 2, 0.004, 0.004, 0.003
!! ...
!!
!! NOTES
!!
!!  The user will likely never call this routine directly.  Developers will need to
!!  understand it if output format to logfile is changed.
!!
!!***

subroutine Logfile_writeGatherCSV(strArr, length, dim, strLen, numHeaders, reduced, separateFiles)

  use Logfile_data, ONLY : log_globalMe,  log_lun, log_fileOpen, log_globalNumProcs, log_fileName   
  use Logfile_interface, ONLY : Logfile_break, Logfile_close, &
    Logfile_open
!   use Driver_interface, ONLY : Driver_getMype, Driver_getComm, Driver_getNumProcs


  implicit none

#include "constants.h"

  integer, intent(in)                                  :: length, dim, strLen, numHeaders
  character(len=MAX_STRING_LENGTH), intent(in), dimension(length,dim)  :: strArr
  logical, optional, intent(IN)                        :: reduced
  logical, optional, intent(IN)                        :: separateFiles
  character(len=MAX_STRING_LENGTH)                     :: indentStr, tag, filename

  integer  :: i, j, tmpLen
  logical  :: doreduced, doseparate
  integer :: logUnit
  logical :: logUnitLocal=.false.

!   call Driver_getNumProcs(GLOBAL_COMM,log_globalNumProcs)
   filename = trim(log_fileName) // trim(".csv")
   open(20, file = trim(filename), action = "write")
   if (log_globalMe .eq. MASTER_PE) then
      do i=0, length-numHeaders
         write(20 ,fmt="(G0, a)", advance="no") adjustl(trim(strArr(numHeaders+i, 2)))
         write(20 ,fmt="(G0, a)", advance="no") ","
         !  all values below are populated in tmr_buildSummary.F90
         write(20 ,fmt="(1x,a)", advance="no") trim(adjustl(strArr(numHeaders+i, 7)))
         write(20 ,fmt="(G0, a)", advance="no") ","
         write(20 ,fmt="(1x,a)", advance="no") trim(adjustl(strArr(numHeaders+i, 8)))
         write(20 ,fmt="(G0, a)", advance="no") ","

         do j=1, log_globalNumProcs
            if (j .eq. log_globalNumProcs) then
               write(20,fmt="(1x,a)", advance="no") trim(adjustl(strArr(numHeaders+i, 8+j)))
            else
               write(20,fmt="(1x, a)", advance="no") trim(adjustl(strArr(numHeaders+i, 8+j)))
         write(20 ,fmt="(G0, a)", advance="no") ","
            endif
         enddo
         write(20,*)
      enddo
   endif

end subroutine Logfile_writeGatherCSV