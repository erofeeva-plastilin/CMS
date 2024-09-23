# Sunflover Mildew
## 1.Получение маркеров
### 1.1 Фильтрация VCF файла
**Активация среды, где есть vcftools**:
```
source ~/.bashrc
conda activate GWAS-PIPELINE
```
**Фильтрация c параметрами maf 0.01, max-missing 0.95** 
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
**Активация среды, где есть plink:**
```
source ~/.bashrc
conda activate GWAS-PIPELINE
```
**Фильтрация с размером окна 100кб, шагом 5 и r^2 0,5  
Добавление уникального ID в VCF файл, чтобы избежать ошибок**
```
plink2 --vcf filt_vcf.recode.vcf --set-all-var-ids @:#\$r\$a --recode vcf --out updated_vcf --allow-extra-chr
plink2 --vcf filt_vcf.recode.vcf --set-all-var-ids @:#\$r\$a --indep-pairwise 100 5 0.5 --out ld --allow-extra-chr
plink2 --vcf updated_vcf.vcf  --extract ld.prune.in --allow-extra-chr --make-pgen --out pruneddata
plink2 --pfile pruneddata --recode vcf --out pruned --allow-extra-chr
```
Лог:
515639 MiB RAM detected, ~328125 available; reserving 257819 MiB for main  
workspace.  
Using up to 40 threads (change this with --threads).  
--vcf: 13722663 variants scanned.  
--vcf: ld-temporary.pgen + ld-temporary.pvar.zst + ld-temporary.psam written.  
82 samples (0 females, 0 males, 82 ambiguous; 82 founders) loaded from  
ld-temporary.psam.  
13722663 variants loaded from ld-temporary.pvar.zst.  
Note: No phenotype data present.  
Calculating allele frequencies... done.  
--indep-pairwise (17 compute threads): 7942018/13722663 variants removed.  
**Осталось 5780645 значимых SNP**  

## 1.3 GWAS - Gapit
**Активация среды, где есть gapit:**
```
source ~/.bashrc
conda activate /mnt/users/grigorieval/miniconda3/envs/gapit
```
**Создание файла с фенотипами:**
* **Исходный vcf  файл получен на основе данных образцов: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP095974&o=acc_s%3Aa, метаданные взяты оттуда**
* **Метаданные соотнесены и отфильтрованны с используемыми в анализе образцами**
* **Отфильтрованные метаданные соотнесены с данными фенотипирования: https://www.ag.ndsu.edu/fss/ndsu-varieties/usda-sunflower-inbred-lines**
* **Составлена сводная таблица с фенотипами: https://docs.google.com/spreadsheets/d/1B0p3dcaxF79B5V6aW1dlWlSZO0Mbgq-D9E_rFblMnaY/edit?gid=437490187#gid=437490187**

**Запуск gapit**
```
Rscript ./run_GAPIT.R pruned.vcf phenotype.Mildew.tsv /mnt/users/erofeevan/Drafts/Drafts/Help_Mildew
```
Файл с фенотипами: https://github.com/erofeeva-plastilin/CMS/blob/Mildew/GWAS_phenotype/Mildew.tsv, файл с легендой: https://github.com/erofeeva-plastilin/CMS/blob/Mildew/GWAS_phenotype/legend.txt
**Результаты gapit:**  
![image](https://github.com/user-attachments/assets/d1878c79-0f50-435a-a9d9-816f23379361)


# Sunflover Rust and Mildew: Sunflower pan-genome analysis shows that hybridization altered gene content and disease resistance + USDA Sunflower Inbred Lines
VCF файл взят тут: https://www.sunflowergenome.org/pangenome-data/  
В нем 493 образца, фильтрация, чтобы оставить только культивируемые (PPN001-PPN289, кроме PPN042, 046, 083, 085, 136, 172, 285):
```
 bcftools view -S line.txt HelianthusVariants.vcf.gz -o cultivar.vcf.gz -O z
```
## Mildew
### Индексация  
```
bcftools index cultivar.vcf.gz
```
### Фильтрация только нужных snp. 
Нужные snp взяты в таблице 9: "Table S9: Summary table of GWAS results for downy mildew resistance in the SAM
population (n = 239 biologically independent samples) based on 675,291 SNP dataset. For
each SNP, attribution to QTL next to the SNP name, its chromosome and position, the effect
size of the major allele (beta) and its standard error (SE), the p-value and significance
following the simpleM correction (*) and permutation test (**), candidate gene based on the
HA412-HO.v1.1 annotation (gene name), likely candidate gene function, and candidate gene
coordinates on the HA412-HO reference genome are provided.":  https://static-content.springer.com/esm/art%3A10.1038%2Fs41477-018-0329-0/MediaObjects/41477_2018_329_MOESM1_ESM.pdf
```
bcftools view -R snp_pos.txt cultivar.vcf.gz -o filtered_cultivar.vcf
```
### Подсчет частот аллелей  
```
source ~/.bashrc
conda activate GWAS-PIPELINE
vcftools --vcf filtered_cultivar.vcf --freq --out allele_frequencies
```
### Перенос на сборку HanXRQr2.0 значимых SNP
```
# Скачивание исходного генома:
wget https://www.heliagene.org/HA412.v1.1.bronze.20141015/static/sequences_files/HA412.v1.1.bronze.20141015.fasta.gz
gunzip HA412.v1.1.bronze.20141015.fasta.gz
# Активация среды, где есть snplift:
source ~/.bashrc
conda activate snplift_env
#Меняем snplift_config.sh
# Input files
export OLD_GENOME="/mnt/users/erofeevan/Mildew/HA412.v1.1.bronze.20141015.fasta"
export NEW_GENOME="/mnt/reference/genomes/heliantus_annuus/GCF_002127325.2/GCF_002127325.2_HanXRQr2.0-SUNRISE_genomic.fna"
# Output files
export INPUT_FILE="/mnt/users/erofeevan/Mildew/filtered_cultivar.vcf"
export OUTPUT_FILE="/mnt/users/erofeevan/Mildew/filtered_cultivar_HanXRQr2.0-SUNRISE.vcf"
# Поменять флаг - если геном не проиндексирован, то надо флаг 0!
#Запустить скрипт
time ./snplift 02_infos/snplift_config.sh
```
Лог:
```
SNPLift: Number of SNPs treated at each step

43      Positions
43      Scores
41      Transferable

SNPLift: Percentage of transferred SNPs:

95.3488%

SNPLift: Run completed

End: 20240923_104510

real    1m28.226s
user    0m25.978s
sys     0m39.056s
```
Результаты: https://docs.google.com/spreadsheets/d/11D7-WJqfNM710nf4aFO29UhfeMguLUZKwjyEQFpdTlo/edit?gid=1270894803#gid=1270894803
## Rust: Sunflower pan-genome analysis shows that hybridization altered gene content and disease resistance + USDA Sunflower Inbred Lines
VCF файл получен из статьи 'Sunflower pan-genome analysis shows that hybridization altered gene content and disease resistance': https://www.sunflowergenome.org/pangenome-data/, в нем оставлены только культивируемые линии: ```/mnt/users/erofeevan/Mildew/cultivar.vcf.gz```
Фенотипы собраны тут:https://docs.google.com/spreadsheets/d/1zaGA4b8CNhZcoSkttB3r_oZKMer3AJ2GTdywpt1HsNc/edit?gid=18979323#gid=18979323   
Файл с фенотипами: https://github.com/erofeeva-plastilin/CMS/blob/Mildew/GWAS_phenotype/Rust.tsv, файл с легендой: https://github.com/erofeeva-plastilin/CMS/blob/Mildew/GWAS_phenotype/legend.txt
Фильтрация vcf файла:   
В файле 5830734 snp, он уже отфильтрован
