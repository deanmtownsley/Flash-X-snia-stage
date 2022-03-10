#
# Use
#    source THIS_FILE
# or
#    . THIS_FILE
# from bash; other Bourne shells or bash invoked as /bin/sh are probably
# lacking some required features.

pathmungeany () {
# Usage: pathmungeany COMPONENT [WHERE] [VARNAME]
#  COMPONENT :- the component(s) to (maybe) add to the path
#  WHERE :- "after" means: append unless present anywhere
#           "first" means: prepend unless already present at beginning
#           anything else means: prepend unless already present anywhere
#  VARNAME :- the name of the environment variable to test and (maybe) modify
#             if empty, defaults to PATH
# Examples:
#    pathmungeany ${MPI}/bin
#    pathmungeany ${MPI}/lib first LD_LIBRARY_PATH
#    pathmungeany ${MPI}/share/man before MANPATH
# Notes:
#  This functions expects that COMPONENT and VARNAME (if given and nonempty)
#  have sane values.  Shit may well happen if special characters like
#  '$', ';', '(', ')', '|', etc. etc. are present.
    if [ -z "$3" ]; then
	varname="PATH"
    else
	varname="$3"
    fi
    if [ -f /bin/egrep -a -x /bin/egrep ]; then
	EGREP=/bin/egrep
    else
	EGREP=/usr/bin/egrep
    fi

    unset previous_pathstring
    eval "declare -p '${varname}'>/dev/null 2>&1" \
	 "&&previous_pathstring="'"${'"${varname}"'}"'
    eval 'test -z "${'"${varname}"'+set}"' && unset previous_pathstring

    if [ -z "${previous_pathstring+set}" ] ; then
        eval "export ${varname}"'=$1'
    elif [ "$2" = "first" ] && 
	! echo "${previous_pathstring}" |($EGREP -q  "^$1($|:)") ; then
        eval "export ${varname}"'=$1:${'"${varname}"'}'
    elif ! echo "${previous_pathstring}" |($EGREP -q "(^|:)$1($|:)") ; then
        if [ "$2" = "after" ] ; then
            eval "export ${varname}"'=${'"${varname}"'}:$1'
        else
            eval "export ${varname}"'=$1:${'"${varname}"'}'
        fi
    fi
}
