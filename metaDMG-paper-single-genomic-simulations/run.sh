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

# quick=false
quick=true
threads=1

# cores=4
cores=10


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



# "0.035": 1%
# "0.065": 2%
# "0.162": 5%
# "0.31": 10%
# "0.472": 15%
# "0.633": 20%
# "0.96": 30%

if [ "$quick" = true ] ; then
    declare -a damages=("0.162")
    # declare -a Nreads=("100" "1000")
    declare -a Nreads=("1000000")
    declare -a length_means=("60")
    declare -a length_means=("35" "60" "90")
    ids=`seq 0 1`
    cores=1

else
    declare -a damages=("0.0" "0.035" "0.065" "0.162" "0.31" "0.472" "0.633" "0.96")
    # declare -a damages=("0.035" "0.065" "0.162" "0.31" "0.472" "0.633")
    declare -a Nreads=("10" "25" "50" "100" "250" "500" "1000" "2500" "5000" "10000" "25000" "50000" "100000")
    # declare -a Nreads=("10" "100" "1000" "10000" "100000")
    declare -a length_means=("60")
    #declare -a length_means=("35" "60" "90")
    ids=`seq 0 99`
fi

length_std=10
quality_scores=Test_Examples/AccFreqL150R1.txt

if [ ! -f $genome.1.bt2 ]; then
    echo "got here"
    #bowtie2-build $genome $genome
fi



fastq_dir=fastq/$species
bam_dir=bam/$species
data_dir=data/$species

mkdir -p $fastq_dir
mkdir -p $bam_dir

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


function simulate_fastq {


    if ! [[ -f $fastq.fq  &&  -s $fastq.fq  ]] # does not exist or is empty
    then
        if (( $(echo "$damage > 0.0001" |bc -l) ));
        then
            # briggs="-b 0.024,0.36,$damage,0.0097"
            briggs="-m b,0.024,0.36,$damage,0.0097"
        else
            briggs=""
        fi

        lognorm_mean=$(compute_lognormal_mean)
        lognorm_std=$(compute_lognormal_std)

        # args="-i $genome -t $threads -r $Nread -ld LogNorm,$lognorm_mean,$lognorm_std -s $id -seq SE -f fq -q1 $quality_scores $briggs -o $fastq"
        # args="-i $genome -t $threads -r $Nread -ld LogNorm,$lognorm_mean,$lognorm_std -seq SE -f fa $briggs -o $fastq"
        args="-i $genome -t $threads -r $Nread -ld LogNorm,$lognorm_mean,$lognorm_std -seq SE -f fq -q1 $quality_scores $briggs -o $fastq"
        ./ngsngs $args
        # echo ./ngsngs $args
    fi
}

function make_bam_file {
    if ! [[ -f $bam  &&  -s $bam  ]] # does not exist or is empty
    then
        bowtie2 -x $genome -q $fastq.fq --threads $threads --no-unal | samtools view -bS -@ $threads - | samtools sort -n -@ $threads -> $bam
    fi

}




function simulate_fastq_and_bam {
    simulate_fastq
    make_bam_file
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
                fastq=./$fastq_dir/$basename
                bam=./$bam_dir/$basename.bam

                if ! [[ "$cores" -eq 1 ]];
                then
                    run_with_lock simulate_fastq_and_bam
                else
                    echo "$damage, $Nread, $id"
                    simulate_fastq
                    make_bam_file
                fi

                let COUNTER++

            done
        done
    done
done

wait

echo $COUNTER



# #%%
cores=5
configfile=config-$species-mini.yaml
poetry run metaDMG config ./$bam_dir/sim-$species-*-*-*-?.bam --damage-mode global --bayesian --long-name --parallel-samples $cores --overwrite --config-file $configfile --output-dir $data_dir
poetry run metaDMG compute $configfile

# #%%
cores=5
configfile=config-$species-mini.yaml
poetry run metaDMG config ./$bam_dir/sim-$species-*-*-*-1?.bam --damage-mode global --bayesian --long-name --parallel-samples $cores --overwrite --config-file $configfile --output-dir $data_dir
poetry run metaDMG compute $configfile

# #%%
cores=5
configfile=config-$species-mini.yaml
poetry run metaDMG config ./$bam_dir/sim-$species-*-*-*-2?.bam --damage-mode global --bayesian --long-name --parallel-samples $cores --overwrite --config-file $configfile --output-dir $data_dir
poetry run metaDMG compute $configfile


# #%%
cores=5
configfile=config-$species-mini.yaml
poetry run metaDMG config ./$bam_dir/sim-$species-*-*-*-3?.bam --damage-mode global --bayesian --long-name --parallel-samples $cores --overwrite --config-file $configfile --output-dir $data_dir
poetry run metaDMG compute $configfile


# #%%
cores=5
configfile=config-$species.yaml
poetry run metaDMG config ./$bam_dir/sim-$species*.bam --damage-mode global --bayesian --long-name --parallel-samples $cores --overwrite --config-file $configfile --output-dir $data_dir
poetry run metaDMG compute $configfile


# #%%
cores=5
configfile=config-$species.yaml
poetry run metaDMG config ./$bam_dir/sim-$species*.bam --damage-mode global --bayesian --long-name --parallel-samples $cores --overwrite --config-file $configfile --output-dir $data_dir
poetry run metaDMG compute $configfile


echo "FINISHED"


# # mkdir -p mapDamage
# # mapDamage -i $bamfile -r $db --no-stats -d mapDamage/$(basename $file).$(basename $db).AFG --merge-libraries
# # ./metaDMG-cpp getdamage --minlength 10 --printlength 50 --threads 8 $bam | cut -f7-8 | sed -n -e 1,2p -e 51,53p  -e 102,110p
