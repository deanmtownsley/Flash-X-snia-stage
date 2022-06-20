#ifndef __FLASH_MOL_H
#define __FLASH_MOL_H

#define MOL_INVALID -1

/*
    MoL function types

    MOL_RHS_EXPLICIT     : RHS for (slow) explicit terms
    MOL_RHS_IMPLICIT     : RHS for (slow) implicit terms
    MOL_RHS_FAST         : RHS for (fast) explicit terms
    MOL_IMPLICIT_UPDATE  : Implicit updates
    MOL_POST_UPDATE      : Post-update (slow) per-stage
    MOL_POST_UPDATE_FAST : Post-update (fast) per-stage
*/
#define MOL_RHS_EXPLICIT 1
#define MOL_RHS_IMPLICIT 2
#define MOL_RHS_FAST 3
#define MOL_IMPLICIT_UPDATE 4
#define MOL_POST_UPDATE 5
#define MOL_POST_UPDATE_FAST 6

/*
    MoL data-structs

    MOL_EVOLVED : Current evolved-variable state (forwards to UNK)
    MOL_INITIAL : Evolved variable state at the start of a timestep
    MOL_RHS     : Current RHS to fill for evolved variables
*/
#define MOL_EVOLVED 0
#define MOL_INITIAL 1
#define MOL_RHS 2

/*
    MoL messaging verbosity

    MOL_VERBOSE_ERROR  : Only error messages
    MOL_VERBOSE_WARN   : Only error and warning messages
    MOL_VERBOSE_STATUS : All messages
*/
#define MOL_VERBOSE_ERROR 1
#define MOL_VERBOSE_WARN 2
#define MOL_VERBOSE_STATUS 3

#endif
