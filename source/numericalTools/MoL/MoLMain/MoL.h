/**
   @copyright Copyright 2023 UChicago Argonne, LLC and contributors

   @par License
   @parblock
      Licensed under the Apache License, Version 2.0 (the "License");
      you may not use this file except in compliance with the License.

      Unless required by applicable law or agreed to in writing, software
      distributed under the License is distributed on an "AS IS" BASIS,
      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
      See the License for the specific language governing permissions and
      limitations under the License.
   @endparblock

   @file
   @ingroup MoL

   @brief Commonly used defines for MoL data structures, procedures, and
          verbosity setttings
*/

#ifndef __FLASH_MOL_H
#define __FLASH_MOL_H

#define MOL_INVALID -1 /**< @brief Used to refer to any invalid and/or unused structure */

/**
   @name MoL function types
   @{
*/
#define MOL_RHS_EXPLICIT 1     /**< @brief RHS for (slow) explicit terms */
#define MOL_RHS_IMPLICIT 2     /**< @brief RHS for (slow) implicit terms */
#define MOL_RHS_FAST 3         /**< @brief RHS for (fast) explicit terms */
#define MOL_IMPLICIT_UPDATE 4  /**< @brief Implicit updates */
#define MOL_POST_UPDATE 5      /**< @brief Post-update (slow) per-stage */
#define MOL_POST_UPDATE_FAST 6 /**< @brief Post-update (fast) per-stage */
/** @} */

/**
   @name MoL data-structs
   @{
*/
#define MOL_EVOLVED 0 /**< @brief Current evolved-variable state (forwards to UNK) */
#define MOL_INITIAL 1 /**< @brief Evolved variable state at the start of a timestep */
#define MOL_RHS 2     /**< @brief Current RHS to fill for evolved variables */
/** @} */

/**
   @name MoL messaging verbosity
   @{
*/
#define MOL_VERBOSE_ERROR 1  /**< @brief Only error messages */
#define MOL_VERBOSE_WARN 2   /**< @brief Only error and warning messages */
#define MOL_VERBOSE_STATUS 3 /**< @brief All messages */
/** @} */

#endif
