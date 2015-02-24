# Count up occurrences of each functional category in a file
# (or anything else in the second column of a tab-delimed file
cut -f2 $1 | sort | uniq -c
