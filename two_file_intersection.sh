# Print all lines in file2 whose 1st column is also in the 1st column of file 1
# Should not be sensitive to order, as the first files contents are stored in
# an array

awk 'NR==FNR{a[$1];next}($1 in a)' $1 $2
