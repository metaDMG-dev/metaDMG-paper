# This repo is a workflow for setting up a metagenomic analysis of shotgun data and mapping this against the refseq plastid, mitochondrion and phylo norway chloroplast databases.

RefSeq release number 211 downloaded 04 april 2022

## download plastids
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/*genomic.fna.gz

## download mitos
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/*genomic.fna.gz

## collect PhyloNorway cpDNAs
cat /willerslev/edna/database/phyloN_refseq_cpDNA/norwaydb_4.fa > /willerslev/edna/database/organelles_04042022/PhyloNorway_cpDNA_Wang_etal.fa &
## download fungi
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/*genomic.fna.gz
## download viruses
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/*genomic.fna.gz
## download archaea
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/archaea/*genomic.fna.gz

# building taxonomy files

## list of custom arctic cpDNAs accession 2 taxid (ncbi format)
cp /willerslev/edna/database/phyloN_refseq_cpDNA/norway_taxid_accessionNo_4.txt .
## downloading the NCBI taxonomy files
mkdir tax_dmp
cd tax_dmp/
## list of custom arctic cpDNAs accession 2 taxid (ncbi format)
cp /willerslev/edna/database/phyloN_refseq_cpDNA/norway_taxid_accessionNo_4.txt.gz .
## download taxonomy from ncbi
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip
unzip new_taxdump.zip
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
## concatenate files
zcat nucl_gb.accession2taxid.gz norway_taxid_accessionNo_4.txt.gz | gzip > nucl_gb.accession2taxid_norway4.gz
## merge and unzip refseq fastas, remove original files
zcat *fna.gz >> refeq_sequences.fa
rm *fna.gz
## merge all fastas
cat *fa > refseq211_small.fa
## as identical accessions (mitos) can be found in the fungi and mitochondrial databases, we here deduplicate the merged fasta
perl dedup_fasta.pl refseq211_small.fa > refseq211_small_dedup.fa
## build bowtie2 databases
bowtie2-build refseq211_small_dedup.fa refseq211_small_dedup.fa --threads 15 --large-index
cd ../

## merge only cp and mtDNA fastas
cat mitochondrion_refseq.fa plastid_refseq.fa norwaydb_4.fa > mito_cpDNArefseq211.fa &
## build bowtie2 database
bowtie2-build mito_cpDNArefseq211.fa mito_cpDNArefseq211.fa --threads 10
cd ../

# Mapping data against database
for file in $(pwd)/*.fq
do
source /home/npl206/miniconda3/etc/profile.d/conda.sh
conda activate metagen_session

for DB in /willerslev/edna/database/organelles_04042022/refseq211_small_dedup.fa
do
echo Mapping $file against $DB
nice -2 bowtie2 --threads 80 -k 100 -x $DB -U $file --no-unal | samtools view -bS - > $file.$(basename $DB).bam
done
done &> refseq211_small_dedup.log.txt

for DB in /willerslev/edna/database/Homo_sapiens/homo_sapiens.fa
do
echo Mapping $file against $DB
nice -n 5 bowtie2 --threads 50 -x $DB -U $file --no-unal | samtools view -bS - > $file.$(basename $DB).bam
samtools sort $file.$(basename $DB).bam > $file.$(basename $DB).sort.bam
samtools index
done

done &> HomoS_refseq211_small.log.txt

# Mapping data against database
for file in $(pwd)/*.fq
do
source /home/npl206/miniconda3/etc/profile.d/conda.sh
conda activate metagen_session

for DB in /willerslev/users-shared/science-snm-willerslev-npl206/Mex_Cave/model_simulations/reference_dbs/Homo_sapiens_mtDNA
do
echo Mapping $bname against $DB
nice -n 5 bowtie2 --threads 10 -x $DB -U $file --no-unal | samtools view -bS - > $file.$(basename $DB).bam
done

done

conda activate bam-filter2

for file in *refseq211_small_dedup.fa.bam
do
samtools sort $file > $file.sort.bam
filterBAM --min-read-count 50 --min-expected-breadth-ratio 0.60 --min-read-ani 95 --sort-memory 2G --threads 40 $file.sort.bam
done &> filterBAM.log.txt

#### activate conda environment
conda activate metaDMG8
#### generate config.yaml file and specify files, taxonomy paths etc. also in the config file other lca and metaDMG settings can be set.
nam=/willerslev/edna/database/organelles_04042022/tax_dmp/names.dmp
nod=/willerslev/edna/ncbi_taxonomy3dec2020/nodes.dmp
acc=/willerslev/edna/ncbi_taxonomy3dec2020/combined_taxid_accssionNO_20201120.gz
# for the dandy system
nam=/projects/lundbeck/people/npl206/databases/ncbi_taxonomy3dec2020/names.dmp
nod=/projects/lundbeck/people/npl206/databases/ncbi_taxonomy3dec2020/nodes.dmp
acc=/projects/lundbeck/people/npl206/databases/ncbi_taxonomy3dec2020/combined_taxid_accssionNO_20201120.gz

# create config file
metaDMG config *.sam.gz --names $nam --nodes $nod --acc2tax $acc

change path to metaDMG-lca: to /willerslev/edna/metadamage/metadamage or /willerslev/edna/programmes/metadamage/metadamage for metaDMG7

#### OBS remember to add the lca path to metadamage thorfinns metadamage /willerslev/edna/metadamage/metadamage in the config file
metaDMG compute config.yaml

metaDMG convert --output Devlin_metaDMGout.csv --add-fit-predictions

for file in $(pwd)/*.fq
do
source /home/npl206/miniconda3/etc/profile.d/conda.sh
conda activate metagen_session

for DB in /willerslev/users-shared/science-snm-willerslev-npl206/Mex_Cave/model_simulations/reference_dbs/Homo_sapiens_mtDNA
do
echo Mapping $file against $DB
bowtie2 --threads 10 -x $DB -U $file --no-unal | samtools view -bS - > $file.$(basename $DB).bam
done

done &> human_mtDNA.log.txt
