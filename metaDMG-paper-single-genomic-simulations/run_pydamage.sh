# bowtie2-build Mycobacterium_leprae.fa Mycobacterium_leprae.fa

species=$1
# species="homo"
# species="betula"
# species="GC-low"
# species="GC-mid"
# species="GC-high"

# quick=false
quick=true

# cores=4
cores=10


if [ "$quick" = true ] ; then
    declare -a damages=("0.0" "0.31")
    # declare -a Nreads=("100" "1000" "1000" "10000" "100000")
    declare -a Nreads=("100" "1000")
    declare -a length_means=("60")
    ids=`seq 0 1`

else
    declare -a damages=("0.0" "0.035" "0.065" "0.162" "0.31" "0.472" "0.633" "0.96")
    declare -a Nreads=("10" "25" "50" "100" "250" "500" "1000" "2500" "5000" "10000" "25000" "50000" "100000")
    declare -a length_means=("60")
    #declare -a length_means=("35" "60" "90")
    ids=`seq 0 99`
fi

bam_dir=bam/$species
pydamage_dir=pydamage/$species

mkdir -p $bam_dir
mkdir -p $pydamage_dir

COUNTER=0


function run_pydamage {
    if ! [[ -f $pydamage_out  &&  -s $pydamage_out  ]] # does not exist or is empty
    then
        samtools sort $bam > $bam_sorted
        samtools index $bam_sorted
        poetry run pydamage --outdir $pydamage_tmp_dir analyze $bam_sorted

        mv $pydamage_tmp_dir/pydamage_results.csv $pydamage_out
        rm -r $pydamage_tmp_dir
        rm $bam_sorted
        rm $bam_sorted.bai
    fi

}


# initialize a semaphore with a given number of tokens
open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done

}

# run the given command asynchronously and pop/push tokens
run_with_lock(){
    local x
    # this read waits until there is something to read
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
     ( "$@"; )
    # push the return code of the command to the semaphore
    printf '%.3d' $? >&3
    )&
}



if ! [[ "$cores" -eq 1 ]];
then
    open_sem $cores
fi

for length_mean in "${length_means[@]}"
do
    for id in $ids;
    do
        for damage in "${damages[@]}"
        do
            for Nread in "${Nreads[@]}"
            do

                basename=sim-$species-$damage-$Nread-$length_mean-$id
                bam=$bam_dir/$basename.bam
                bam_sorted=$bam_dir/$basename.sorted.bam

                pydamage_tmp_dir=$pydamage_dir/tmp/$basename
                pydamage_out=$pydamage_dir/$basename.pydamage.csv


                if ! [[ "$cores" -eq 1 ]];
                then
                    # echo "$basename $bam_dir $pydamage_dir"
                    run_with_lock run_pydamage

                else
                    # echo "$basename $bam_dir $pydamage_dir"
                    run_pydamage
                    # echo "$damage, $Nread, $id"
                fi

                let COUNTER++

            done
        done
    done
done

wait

echo $COUNTER

echo "FINISHED"