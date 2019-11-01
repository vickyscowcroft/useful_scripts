#!/bin/csh -f

## prints record 9 onwards - useful for filenames with spaces in

set outname = $1

ls -ltr | gawk '{for (i=9; i<NF; i++) printf $i " "; print $NF}' > $outname
