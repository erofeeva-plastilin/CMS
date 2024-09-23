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
**Результаты gapit:**  
![image](https://github.com/user-attachments/assets/d1878c79-0f50-435a-a9d9-816f23379361)


# Sunflover Rust and Mildew
VCF файл взят тут: https://www.sunflowergenome.org/pangenome-data/  
В нем 493 образца, фильтрация, чтобы оставить только культивируемые (PPN001-PPN289, кроме PPN042, 046, 083, 085, 136, 172, 285):
```
 bcftools view -S line.txt HelianthusVariants.vcf.gz -o cultivar.vcf.gz -O z
```
## Mildew
Индексация  
```
bcftools index cultivar.vcf.gz
```
Фильтрация только нужных snp  
```
bcftools view -R snp_pos.txt cultivar.vcf.gz -o filtered_cultivar.vcf
```
Подсчет частот аллелей  
```
source ~/.bashrc
conda activate GWAS-PIPELINE
vcftools --vcf filtered_cultivar.vcf --freq --out allele_frequencies
```
## Rust
Фенотипы собраны тут:https://docs.google.com/spreadsheets/d/1zaGA4b8CNhZcoSkttB3r_oZKMer3AJ2GTdywpt1HsNc/edit?gid=18979323#gid=18979323
