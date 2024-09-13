# Sunflover CMS
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

## 1.3 Gapit
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
Rscript ./run_GAPIT.R pruned.vcf phenotype.CMS.tsv /mnt/users/erofeevan/Drafts/Drafts/Help_CSM
```
**Результаты gapit:**  

![image](https://github.com/user-attachments/assets/4cfb54de-8016-4e13-92d6-4359753dc511)
**Выделены маркеры на 10 и 13 хромосомах, что соответствует статьям.
Отфильтрованные значимые маркеры (желтым выделены маркеры на 13 хромосоме):**  

![image](https://github.com/user-attachments/assets/e014bc80-c6b2-49c6-94e6-948e92d845fa)

**Все промежуточные файлы хранятся по адресу:**
```
/mnt/users/erofeevan/Drafts/Drafts/Help_CSM
```
## 2.Соотнесение координат
### 2.1 Фильтрация только 13 хромосомы
**Активация среды, где есть vcftools:**
```
source ~/.bashrc
conda activate GWAS-PIPELINE
```
**Отфитровываю только 13 хромосому (обозначение можно узнать в референсном геноме: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002127325.2/ и https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002127325.1/):**
```
vcftools --vcf pruned.vcf --chr NC_035445.2 --recode --out 13_chr
```
Лог:
Parameters as interpreted:  
        --vcf pruned.vcf  
        --chr NC_035445.2  
        --out 13_chr  
        --recode  

After filtering, kept 82 out of 82 Individuals  
Outputting VCF file...  
After filtering, kept 344266 out of a possible 5780645 Sites  
Run Time = 24.00 seconds  

### 2.2 Snplift
**Активация среды, где есть snplift:**
```
source ~/.bashrc
conda activate snplift_env
```
**Меняем snplift_config.sh**
```
# Input files
export OLD_GENOME="/mnt/reference/genomes/heliantus_annuus/GCF_002127325.2/GCF_002127325.2_HanXRQr2.0-SUNRISE_genomic.fna"
export NEW_GENOME="/mnt/users/grigorieval/helianthus/ncbi_dataset/data/GCF_002127325.1/GCF_002127325.1_HanXRQr1.0_genomic.fna"
# Output files
export INPUT_FILE="/mnt/users/erofeevan/13_chr.recode.vcf"
export OUTPUT_FILE="/mnt/users/erofeevan/13_chr_ver1.0.recode.vcf"
# Поменять флаг - если геном не проиндексирован, то надо флаг 0!
```
**Запустить скрипт**
```
time ./snplift 02_infos/snplift_config.sh
```
Лог:
SNPLift: Number of SNPs treated at each step 
339299  Positions 
339320  Scores  
327999  Transferable  
SNPLift: Percentage of transferred SNPs:   
96.6696%   
SNPLift: Run completed   
End: 20240913_103710   
real    13m31.426s   
user    11m13.651s  
sys     1m54.512s  
