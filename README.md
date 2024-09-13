# Sunflover CMS
## 1.Получение маркеров
### 1.1 Фильтрация VCF файла
Активация среды, где есть vcftools:
```
source ~/.bashrc
conda activate GWAS-PIPELINE
```
Фильтрация c параметрами maf 0.01, max-missing 0.95 
```
vcftools --vcf /mnt2/results/20231019_WGS_helianthus/vcfs/annotation.vcf --maf 0.01 --max-missing 0.95 --recode --out filt_vcf
```
Лог:
Parameters as interpreted:
        --vcf /mnt2/results/20231019_WGS_helianthus/vcfs/annotation.vcf
        --maf 0.01
        --max-missing 0.95
        --out filt_vcf
        --recode
After filtering, kept 82 out of 82 Individuals
Outputting VCF file...
After filtering, kept 13722663 out of a possible 52147915 Sites
Run Time = 2162.00 seconds
### 1.2 Фильтрация по LD блокам
plink2 --vcf filt_vcf.recode.vcf --set-all-var-ids @:#\$r\$a --recode vcf --out updated_vcf --allow-extra-chr
plink2 --vcf filt_vcf.recode.vcf --set-all-var-ids @:#\$r\$a --indep-pairwise 100 5 0.5 --out ld --allow-extra-chr
plink2 --vcf updated_vcf.vcf  --extract ld.prune.in --allow-extra-chr --make-pgen --out pruneddata
plink2 --pfile pruneddata --recode vcf --out pruned --allow-extra-chr

Rscript ./run_GAPIT.R pruned.vcf phenotype.CMS.tsv /mnt/users/erofeevan/Drafts/Drafts

scp -P 48226 erofeevan@188.170.2.5:/mnt/users/erofeevan/Drafts/Drafts/Help_CMS/* /home


отфитровываю только 13 хромосому (обозначение можно узнать в референсном геноме: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002127325.2/ и https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002127325.1/):
vcftools --vcf pruned.vcf --chr NC_035445.2 --recode --out 13_chr

Parameters as interpreted:
        --vcf pruned.vcf
        --chr NC_035445.2
        --out 13_chr
        --recode

After filtering, kept 82 out of 82 Individuals
Outputting VCF file...
After filtering, kept 344266 out of a possible 5780645 Sites
Run Time = 24.00 seconds

# Input files
export OLD_GENOME="/mnt/reference/genomes/heliantus_annuus/GCF_002127325.2/GCF_002127325.2_HanXRQr2.0-SUNRISE_genomic.fna"
export NEW_GENOME="/mnt/users/grigorieval/helianthus/ncbi_dataset/data/GCF_002127325.1/GCF_002127325.1_HanXRQr1.0_genomic.fna"

# Output files
export INPUT_FILE="/mnt/users/erofeevan/13_chr.recode.vcf"
export OUTPUT_FILE="/mnt/users/erofeevan/13_chr_ver1.0.recode.vcf"
Поменять флаг - если геном не проиндексирован, то надо флаг 0!
chmod +x snplift
#Запустить скрипт
time ./snplift 02_infos/snplift_config.sh

tmux
time ./snplift 02_infos/snplift_config.sh
