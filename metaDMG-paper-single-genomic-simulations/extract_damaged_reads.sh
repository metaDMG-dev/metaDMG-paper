

species=$1
# FILES=$1/*.fq
FILES="./fastq/$species/*.fq"
# FILES="./fastq/sim-*-*-??-*-*.fq"

out=damaged_reads_$species.txt
echo "HEADER" > $out

for file in $FILES
do
  echo "Processing $file file..."
    name="$(basename -- $file)"
    echo $name >> $out
    # samtools view $file | cut -f1 | grep -o 'mod.*' | sort | uniq -c >> $out
    awk 'NR % 4 == 1' $file | grep -o 'mod.*' | sort | uniq -c >> $out
done


