#!/bin/bash
# simple script.  Give it a TAB delimited file $2 with the first column holding
# indexes.  This will pick out the line that matches $1, and return a command
# line assigning the column header names as shell variables whose values are the
# follows in the file $2.

if [ $# -ne 2 ]; then
  echo "Wrong number of arguments in $0 
Expecting two arguments, with syntax like this:

    line_assign.sh IndexValue  TSV-file
    
Where 
    IndexValue is the number in the first column of the row you want to
       extract values from.
    TSV-file is a tab-separated-values file with the first column named
       index, and the other columns named according to shell variables whose
       values you want to set to the values in the row IndexValue
       
Typical usage:

    eval \$(line_assign.sh 7 numbered-units.tsv)

That will set, in the current shell, a series of shell variables whose
names are the column names to values in the
row IndexValue
" > /dev/stderr

exit 1

fi

awk -F"\t" -v LINE=$1 '
  $1 == "index" {for(i=1; i<=NF; i++) vars[i]=$i; next}
  $1 == LINE {for(i=1; i<=NF; i++) printf("%s=\"%s\"; ", vars[i], $i)}
' $2


