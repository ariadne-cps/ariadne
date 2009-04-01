#!/bin/sh

EXECUTABLE=$1

BASENAME=`basename $EXECUTABLE`

# The top-level directory name
TOPDIR=..

# Set the profile directory name
PROFILEDIR=$TOPDIR/profile

# Make the profile directory if it doesn't already exist
if test ! -d $PROFILEDIR; then mkdir $PROFILEDIR; fi

# Get the current working revision
REVISION=`svn info -r HEAD | grep 'Revision' | sed 's/Revision: //'`

# Test if the current working revision is untouched
if svn status ../ | grep '^M' > /dev/null
then
    # Revision has been changed since last checked in
    # Write output to file with name $EXECUTABLE-w$REVISION.NUMBER.log
    NUMBER=1
    LOGFILENAME=$PROFILEDIR/$BASENAME-w$REVISION.0$NUMBER.log

    # Test for previous log file
    while test -f $LOGFILENAME
    do
        NUMBER=$[$NUMBER+1]
        if test "$NUMBER" -lt 10; then
            # Pad number with leading zero to preserve order of filenames
            LOGFILENAME=$PROFILEDIR/$BASENAME-w$REVISION.0$NUMBER.log
        else
            LOGFILENAME=$PROFILEDIR/$BASENAME-w$REVISION.$NUMBER.log
        fi
    done
else
    # Write output to file with name $EXECUTABLE-r$REVISION.log
    LOGFILENAME=$PROFILEDIR/$BASENAME-r$REVISION.log
fi

echo $LOGFILENAME


(echo $HOSTNAME; date; echo; bash $EXECUTABLE) 1>  $LOGFILENAME

