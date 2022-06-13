#ifndef __FLASH_MOL_H
#define __FLASH_MOL_H

#define MOL_INVALID -1

/*
    MoL memory data-structs

    MOL_EVOLVED : Current evolved-variable state (forwards to UNK)
    MOL_INITIAL : Evolved variable state at the start of a timestep
    MOL_RHS     : Current RHS to fill for evolved variables
*/
#define MOL_EVOLVED 0
#define MOL_INITIAL 1
#define MOL_RHS 2

#endif
