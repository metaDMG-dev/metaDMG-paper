
species=$1

damages=("0.0" "0.96")
metaDMG_dirs=("lca" "mismatches" "fit_results" "results")

echo "Moving files for damage $species"

for damage in "${damages[@]}"
do
    echo "Moving old_bam/$species/sim-$species-$damage"
    mv old_bam/$species/sim-$species-$damage-* bam/$species/

    echo "Moving old_fastq/$species/sim-$species-$damage"
    mv old_fastq/$species/sim-$species-$damage-* fastq/$species/

    for metaDMG_dir in "${metaDMG_dirs[@]}"
    do
        echo "Moving old_data/$species/$metaDMG_dir/sim-$species-$damage"
        mv old_data/$species/$metaDMG_dir/sim-$species-$damage-* data/$species/$metaDMG_dir/
    done

done
