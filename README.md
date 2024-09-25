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

### 1.3 GWAS - Gapit
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
## Mildew: Sunflower pan-genome analysis shows that hybridization altered gene content and disease resistance
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
VCF файл получен из статьи 'Sunflower pan-genome analysis shows that hybridization altered gene content and disease resistance': https://www.sunflowergenome.org/pangenome-data/, в нем оставлены только культивируемые линии 
Фенотипы собраны тут:https://docs.google.com/spreadsheets/d/1zaGA4b8CNhZcoSkttB3r_oZKMer3AJ2GTdywpt1HsNc/edit?gid=18979323#gid=18979323     
Файл с фенотипами: https://github.com/erofeeva-plastilin/CMS/blob/Mildew/GWAS_phenotype/Rust.tsv, файл с легендой: https://github.com/erofeeva-plastilin/CMS/blob/Mildew/GWAS_phenotype/legend.txt    
Фильтрация vcf файла:   
В файле 5830734 snp, дополнительно отфильтруем:
**Фильтрация c параметрами maf 0.01, max-missing 0.95**
```
source ~/.bashrc
conda activate GWAS-PIPELINE
vcftools --vcf /mnt/users/erofeevan/Drafts/Drafts/cultivar.vcf --maf 0.01 --max-missing 0.95 --recode --out cultivar_filt_maf
```
Лог
```
Parameters as interpreted:
        --vcf /mnt/users/erofeevan/Drafts/Drafts/cultivar.vcf
        --maf 0.01
        --max-missing 0.95
        --out cultivar_filt_maf
        --recode

After filtering, kept 282 out of 282 Individuals
Outputting VCF file...
After filtering, kept 1238001 out of a possible 5830734 Sites
Run Time = 1501.00 seconds
```
**Добавление уникального ID в VCF файл, чтобы избежать ошибок, и запись генотипа в упрощенном формате**
```
plink2 --vcf cultivar_filt_maf.recode.vcf --set-all-var-ids @:#\$r\$a --recode vcf --out cultivar_filt_maf_with_id --allow-extra-chr --vcf-half-call m
```
Лог:
```
Options in effect:
  --allow-extra-chr
  --export vcf
  --out cultivar_filt_maf_with_id
  --set-all-var-ids @:#$r$a
  --vcf cultivar_filt_maf.recode.vcf
  --vcf-half-call m

Start time: Tue Sep 24 13:05:13 2024
515639 MiB RAM detected, ~250685 available; reserving 250621 MiB for main
workspace.
Using up to 40 threads (change this with --threads).
--vcf: 1238001 variants scanned.
--vcf: cultivar_filt_maf_with_id-temporary.pgen +
cultivar_filt_maf_with_id-temporary.pvar.zst +
cultivar_filt_maf_with_id-temporary.psam written.
282 samples (0 females, 0 males, 282 ambiguous; 282 founders) loaded from
cultivar_filt_maf_with_id-temporary.psam.
1238001 variants loaded from cultivar_filt_maf_with_id-temporary.pvar.zst.
Note: No phenotype data present.
--export vcf to cultivar_filt_maf_with_id.vcf ... done.
End time: Tue Sep 24 13:05:24 2024
```
### GWAS
**Активация среды gapit:**
```
source ~/.bashrc
conda activate /mnt/users/grigorieval/miniconda3/envs/gapit
```
**Запуск gapit**
```
Rscript ./run_GAPIT.R cultivar_filt_maf_with_id.vcf phenotype.Rust.tsv /mnt/users/erofeevan/Drafts/Drafts/Help_Rust
```
Вышло 0 значимых результатов. Никаких SNP нет.  
![image](https://github.com/user-attachments/assets/fbcf7836-72fd-4e08-a0b5-887ac1b56201)

Проводила GWAS 5 раз:   
* На отфильтрованном cultivar_filt_maf_with_id.vcf с 4 разными фенотипами (R, MR, MS, S)
* На отфильтрованном cultivar_filt_maf_with_id.vcf с 3 разными фенотипами (R, MR, MS, S) - MS включен в состав S
* На отфильтрованном cultivar_filt_maf_with_id.vcf с 2 разными фенотипами (R, MR, MS, S) - MS включен в состав S, MR включен в состав R
* На отфильтрованном cultivar_filt_maf_with_id.vcf с 2 разными фенотипами (R, S) - MS и MR удалены
* На нефильтрованном cultivar.vcf с 2 разными фенотипами (R, S) - MS и MR удалены   
Возможные причины: R 51 линия, S 21 линия, неравновестное соотношение; vcf файл плохого качества (исходный содержал 5 млн SNP, но было непонятно, какую обработку он прошел, так как после фильтрации по MAF и maf 0.01 --max-missing 0.95 осталось 1 млн SNP; в исходном файле были проблемы с записью формата, которые были исправлены)
Попробовала перенести отфильтрованный vcf файл ```cultivar_filt_maf.recode.vcf``` со сборки HA412.v1.1 на сборку HanXRQr2.0 и поставить GWAS на нем (чтобы исправить возможные проблемы в конфигурации файла)   
Лог:  
```
1238001 Positions
1236278 Scores
1184931 Transferable

SNPLift: Percentage of transferred SNPs:

95.7132%

SNPLift: Run completed

End: 20240925_103255

real    37m44.479s
user    34m44.782s
sys     3m30.352s
```
GWAS: 
```
Rscript ./run_GAPIT.R cultivar_HanXRQr2.0-SUNRISE.vcf phenotype.Rust.tsv /mnt/users/erofeevan/Drafts/Drafts/Rust # с 2 разными фенотипами (R, S) - MS и MR удалены
```
Вышло 0 значимых результатов. Никаких SNP нет.  
![image](https://github.com/user-attachments/assets/d2c247ea-bfb4-41c4-a75f-58c0a33c6ff3)

Из статей:
* На сегодняшний день у подсолнечника обнаружено в общей сложности 17 генов устойчивости к ржавчине, R1–R5, R10–R12, R13a, R13b, R14–R18, Pu6 и Radv, у 15 из них были обнаружены нанесен на карту в различных областях генома подсолнечника: хромосома 2 (R5); хромосома 8 (R1 и R15); хромосома 11 (R12 и R14); хромосома 13 (R4, R11, R13a, R13b, R16–R18, Pu6 и Radv); и хромосома 14 ( R2): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9455867/
* Устойчивость к ржавчине связана также с генами на 13 хромосоме, которые ассоциированы с большим кластером генов ```The rust resistance genes, R13a and R16, were previously mapped to a 3.4 Mb region at the lower end of sunflower chromosome 13```: https://pubmed.ncbi.nlm.nih.gov/36076914/, ```Previous mapping placed Rf5 at an interval of 5.8 cM on sunflower chromosome 13, distal to a rust resistance gene R11 at a 1.6 cM genetic distance in an SSR map.```:https://pubmed.ncbi.nlm.nih.gov/33437028/, ```Resistance to downy mildew (DM) disease in sunflower is mediated by Pl genes which are known to be effective against the causal fungus, Plasmopara halstedii. Two DM resistance genes, Pl Arg and Pl 8 , are highly effective against P. halstedii races in the USA, and have been previously mapped to the sunflower linkage groups (LGs) 1 and 13, respectively, using simple sequence repeat (SSR) markers.```: https://pubmed.ncbi.nlm.nih.gov/28160079/
* В этих статьях: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4094432/#pone.0098628.s004, https://pubmed.ncbi.nlm.nih.gov/36076914/ и https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7804242/ выявили маркеры, но нет vcf, чтобы их восстановить, но есть позиции на разных сборках. Можно собрать генотипы интересующих линий, сделать свой vcf и попробовать проанализировать.
