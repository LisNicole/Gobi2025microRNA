#!/bin/sh

set -e

while getopts d:a:r:l:f:m:g: flag; do
    case "${flag}" in
        d) directory=${OPTARG};;
        a) adapter_file=${OPTARG};;
        f) family=${OPTARG};;
        r) reference_directory=${OPTARG};;
        m) miRBase_species_abbreviation=${OPTARG};;
        g) genome_species_abbreviation=${OPTARG};;
    esac
done

c_values=(1 2 3 4 5 6 7 8 9 10)
l_values=(1 2 3)

for minimum_libraries in "${l_values[@]}"; do
    for minimum_reads in "${c_values[@]}"; do
        echo "Running script with c=$minimum_reads and l=$minimum_libraries"
        bash /u/halle/ge59puj/home_at/Lotta_Gobi2025/bin/Peak_Calling.sh \
            -d "$directory" \
            -a "$adapter_file" \
            -c "$minimum_reads" \
            -l "$minimum_libraries" \
            -f "$family" \
            -r "$reference_directory" \
            -m "$miRBase_species_abbreviation" \
            -g "$genome_species_abbreviation"

    done
done

for l in "${l_values[@]}"; do

    touch "${directory}l${l}_peaks_per_minread.csv"
    chmod 666 "${directory}l${l}_peaks_per_minread.csv"

    for c in "${c_values[@]}"; do
      peaks=$(grep "Total peaks:"  "${directory}peakcalling.family.summary_c${c}_l${l}.txt" | awk '{print $3}')
      echo "${c},${peaks}" >> "${directory}l${l}_peaks_per_minread.csv"
    done

    awk -F, 'NF==2 {gsub(/^ +| +$/, "", $1); gsub(/^ +| +$/, "", $2); print $1, $2}' "${directory}l${l}_peaks_per_minread.csv" > "${directory}l${l}_cleaned_peaks_per_minread.dat"

    gnuplot -persist <<-EOF
  set terminal png size 800,600
  set output '${directory}l${l}_peaks_histogram.png'
  set title "Peaks per Minimum Reads for l = ${l}"
  set xlabel "Minimum Reads (c)"
  set ylabel "Total Peaks"
  set grid
  set style fill solid
  set boxwidth 0.5
  set xtics 1
  plot '${directory}l${l}_cleaned_peaks_per_minread.dat' using 1:2 with boxes title 'Total Peaks'
EOF

    echo "Histogram saved as ${directory}l${l}_peaks_histogram.png"

done