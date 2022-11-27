# bowtie2-build Mycobacterium_leprae.fa Mycobacterium_leprae.fa

species=$1
# species="homo"
# species="betula"
# species="GC-low"
# species="GC-mid"
# species="GC-high"
# species="contig1k"
# species="contig10k"
# species="contig100k"

threads=10

cores=1
# cores=10


if [[ "$species" == "homo" ]]; then
    genome=./genome/NC_012920.1_genomic.fna
elif [[ "$species" == "betula" ]]; then
    genome=./genome/KX703002.1.fa
elif [[ "$species" == "GC-low" ]]; then
    genome=./genome/GCF_002763915.1_ASM276391v1_genomic.fna
elif [[ "$species" == "GC-mid" ]]; then
    genome=./genome/GCA_900475315.1_43721_G02_genomic.fna
elif [[ "$species" == "GC-high" ]]; then
    genome=./genome/GCA_001929375.1_ASM192937v1_genomic.fna
elif [[ "$species" == "contig1k" ]]; then
    genome=./genome/contig1000.fna
elif [[ "$species" == "contig10k" ]]; then
    genome=./genome/contig10000.fna
elif [[ "$species" == "contig100k" ]]; then
    genome=./genome/contig100000.fna
else
    echo "bad species."
    exit 1
fi



declare -a damages=("0.0" "0.035" "0.065" "0.162" "0.31" "0.472" "0.633" "0.96")
declare -a Nreads=("100000")
# declare -a length_means=("35" "60" "90")
declare -a length_means=("60")
ids=`seq 0 0`

length_std=10

data_dir=true-damage

mkdir -p $data_dir

COUNTER=0


function compute_lognormal_mean {
    mu=$length_mean
    sigma=$length_std
    echo "l($mu^2 / sqrt($sigma^2 + $mu^2))" | bc -l | awk '{printf "%f", $0}'
}

function compute_lognormal_std {
    mu=$length_mean
    sigma=$length_std
    echo "sqrt(l(1 + $sigma^2 / $mu^2))" | bc -l | awk '{printf "%f", $0}'
}


function simulate_fasta {


    if ! [[ -f $fasta.fa  &&  -s $fasta.fa  ]] # does not exist or is empty
    then
        if (( $(echo "$damage > 0.0001" |bc -l) ));
        then
            briggs="-m b,0.024,0.36,$damage,0.0097"
        else
            briggs=""
        fi

        lognorm_mean=$(compute_lognormal_mean)
        lognorm_std=$(compute_lognormal_std)

        args="-i $genome -t $threads -r $Nread -ld LogNorm,$lognorm_mean,$lognorm_std -seq SE -f fa $briggs -o $fasta"
        ./ngsngs $args
        # echo ./ngsngs $args
    fi
}

function make_bam_file {
    if ! [[ -f $bam  &&  -s $bam  ]] # does not exist or is empty
    then
        bowtie2 -x $genome -f $fasta.fa --threads $threads --no-unal | samtools view -bS -@ $threads - | samtools sort -n -@ $threads -> $bam
    fi

}




function simulate_fasta_and_bam {
    simulate_fasta
    make_bam_file
    # ./metaDMG-cpp getdamage --minlength 10 --printlength 50 --threads 8 $bam | cut -f7-8 > $txt
    ./metaDMG-cpp getdamage --minlength 10 --printlength 15 --threads 8 $bam | cut -f7-8 | sed -n -e 1,2p -e 16,18p -e 32,35p
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

                basename=truedamage-$species-$damage-$Nread-$length_mean-$id
                fasta=./$data_dir/$basename
                bam=./$data_dir/$basename.bam
                txt=./$data_dir/$basename.txt

                if ! [[ "$cores" -eq 1 ]];
                then
                    run_with_lock simulate_fasta_and_bam
                else
                    echo "$damage, $Nread, $id"
                    simulate_fasta_and_bam
                    # simulate_fasta
                    # make_bam_file
                fi

                let COUNTER++

            done
        done
    done
done

wait

rm meta.stat
rm meta.bdamage.gz
rm meta.res.gz

echo $COUNTER
