# names=/willerslev/edna/database/organelles_04042022/tax_dmp/names.dmp
names=/willerslev/edna/ncbi_taxonomy3dec2020/names.dmp
nodes=/willerslev/edna/ncbi_taxonomy3dec2020/nodes.dmp
acc2tax=/willerslev/edna/ncbi_taxonomy3dec2020/combined_taxid_accssionNO_20201120.gz

metaDMG config /willerslev/edna/christian_data/NEW_data/*/*.ppm.fq.merged.sort.sam.gz \
--names $names \
--nodes $nodes \
--acc2tax $acc2tax \
--metaDMG-cpp ../thorfinn_metadamage/metadamage \
--custom-database


# edit config file
metaDMG compute config.yaml
