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
   @ingroup physics_Spacetime_SpacetimeMain

   @brief Re-named variable defines for the ADM variables provided in
   SpacetimeMain
*/

#ifndef __FLASH_ADM_H
#define __FLASH_ADM_H

#include "Simulation.h"

#define ALP_VAR SP00_VAR

#define BETAX_VAR SP01_VAR
#define BETAY_VAR SP02_VAR
#define BETAZ_VAR SP03_VAR

#define GXX_VAR SP04_VAR
#define GXY_VAR SP05_VAR
#define GXZ_VAR SP06_VAR
#define GYY_VAR SP07_VAR
#define GYZ_VAR SP08_VAR
#define GZZ_VAR SP09_VAR

#define KXX_VAR SP10_VAR
#define KXY_VAR SP11_VAR
#define KXZ_VAR SP12_VAR
#define KYY_VAR SP13_VAR
#define KYZ_VAR SP14_VAR
#define KZZ_VAR SP15_VAR

#define TE_VAR SP16_VAR

#define TSX_VAR SP17_VAR
#define TSY_VAR SP18_VAR
#define TSZ_VAR SP19_VAR

#define TSXX_VAR SP20_VAR
#define TSXY_VAR SP21_VAR
#define TSXZ_VAR SP22_VAR
#define TSYY_VAR SP23_VAR
#define TSYZ_VAR SP24_VAR
#define TSZZ_VAR SP25_VAR

#endif
