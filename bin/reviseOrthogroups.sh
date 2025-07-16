#!/usr/bin/env bash

prot_reps=$1
cds_reps=$2
prefix=$3

#Define functions
transpose_columns() {
    awk '
        { 
            for (i=1; i<=NF; i++)  {
                a[NR,i] = $i
            }
        }
        NF>p { p = NF }
        END {    
            for(j=1; j<=p; j++) {
                str=a[1,j]
                for(i=2; i<=NR; i++){
                    str=str"\t"a[i,j];
                }
                print str
            }
        }' "$1"
}

filter_fasta() {
    local output_file="$3"
    print_sequence=false
    while IFS= read -r line; do
        if [[ $line == ">"* ]]; then
            accession=$(echo "$line")
            accession=${accession#>}
            if grep -Fxq "$accession" "$2"; then
                print_sequence=true
                echo "$line" >> "$output_file"
            else
                print_sequence=false
            fi
        else
            if $print_sequence; then
                echo "$line" >> "$output_file"
            fi
        fi
    done < "$1"
}

filter_zero_rows() {
    local input_file="$1"
    local output_file="$2"
    local zero_value="$3" # Add a parameter for the zero value (0 or .)
    awk -v zero="$zero_value" '{ # Pass zero_value as an awk variable
        all_zero = 1
        for (i = 2; i <= NF; i++) {
            if ($i != zero) { # Use the awk variable 'zero'
                all_zero = 0
                break
            }
        }
        if (!all_zero) {
            print
        }
    }' "$input_file" > "$output_file"
}

summary=PAV_counts.txt
pav=PAV_table.txt

#create list of orthogroups to search for
list=$(grep "^>"${prefix}"_" ${prot_reps} | sed 's/^>//')
echo $list

#generate summary files
header="isolate"
for name in $list; do
    header="$header\t$name"
done
echo -e "$header" > $summary
echo -e "$header" > $pav

cat *.PAV_counts.txt >> $summary
cat *.PAV_table.txt >> $pav

transpose_columns $summary > "transposed_$summary"
transpose_columns $pav > "transposed_$pav"

#Find non-empty ortholog groups
filter_zero_rows "transposed_$summary" "revised_$summary" "0"
filter_zero_rows "transposed_$pav" "revised_$pav" "."

#representative panel output (generated first run)
tail -n +2 "revised_$summary" | awk '{print $1}' > "pattern.txt"
filter_fasta ${prot_reps} "pattern.txt" "revised_representatives.prot.fasta"
filter_fasta ${cds_reps} "pattern.txt" "revised_representatives.cds.fasta"