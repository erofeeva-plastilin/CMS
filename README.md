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
Rscript ./run_GAPIT.R pruned.vcf phenotype.CMS.tsv /mnt/users/erofeevan/Drafts/Drafts/Help_CSM
```
**Результаты gapit:**  

![image](https://github.com/user-attachments/assets/4cfb54de-8016-4e13-92d6-4359753dc511)

**Выделены маркеры на 10 и 13 хромосомах, что соответствует статьям.**  
**Отфильтрованные значимые маркеры (желтым выделены маркеры на 13 хромосоме):**    

![image](https://github.com/user-attachments/assets/e014bc80-c6b2-49c6-94e6-948e92d845fa)

**Все промежуточные файлы хранятся по адресу:**
```
/mnt/users/erofeevan/Drafts/Drafts/Help_CSM
```
* GAPIT.Association.Filter_GWAS_results.csv
* GAPIT.Association.QQ.MLM.CMS.pdf
* GAPIT.Association.GWAS_Results.MLM.CMS.csv
* GAPIT.Association.Significant_SNPs.MLM.CMS.pdf
* GAPIT.Association.GWAS_StdErr.MLM.CMS.csv
* GAPIT.Association.Vairance_markers.MLM.CMS.csv
* GAPIT.Association.Manhattan_Chro.MLM.CMS.pdf
* GAPIT.Genotype.Distance_R_Chro.pdf
* GAPIT.Association.Manhattan_Geno.MLM.CMS.pdf
* GAPIT.Genotype.Distance.Rsquare.csv
* GAPIT.Association.Manhattans_Symphysic_Legend.pdf
* GAPIT.Genotype.Frequency_MAF.csv
* GAPIT.Association.Manhattans_Symphysic.pdf
* GAPIT.Genotype.Frequency.pdf
* GAPIT.Association.Manhattans_Symphysic_Traitsnames.csv
* GAPIT.Genotype.Kin_Zhang.csv
* GAPIT.Association.Optimum.MLM.CMS.pdf
* GAPIT.Genotype.Kin_Zhang.pdf
* GAPIT.Association.Prediction_results.MLM.CMS.csv
* GAPIT.Genotype.MAF_Heterozosity.pdf
* GAPIT.Association.PVE.MLM.CMS.csv
* GAPIT.Phenotype.View.CMS.pdf

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

**По файлу transfered_positions.tsv перенесла позиции со сборки 2.0 на 1.0, жирным шрифтом выделены позиции значимых SNP на 13 хромосоме в сборке 1.0 (отфильтрованы в порядке возрастания)**  
![image](https://github.com/user-attachments/assets/7e9a9bfe-7c72-44dc-8eaa-e886d69306bf)

**SNP из статьи https://www.mdpi.com/2073-4395/9/2/49#app1-agronomy-09-00049**
![image](https://github.com/user-attachments/assets/4b070f39-17be-4565-9cc9-6281da9bb52b)

**Сравнение**
![image](https://github.com/user-attachments/assets/f7f80f22-8b5e-465d-a5c0-2ac876105e2d)

Красными точками выделены снипы из статьи, зеленым - снипы, обнаруженные в результате анализа. Сверху подписаны позиции снипов, снизу подписаны расстояния между ближайшими снипами из статьи и нашего анализа

**Сравнение в табличном варианте, цветовые обозначения аналогичны**
![image](https://github.com/user-attachments/assets/e0e3b4d7-8f26-49a1-a0ef-d2fc3ed1eb2c)

## 3. Финализация результатов

* Статейных снипов нет в полученном vcf файле :(
* Фильтрация VCF файла для получения только интересующих позиций
* Скрипт для фильтрации pos.py:
```
positions = {
    170763933, 171027789, 172614973, 173681208, 174122679, 
    174663187, 174873488, 174887174, 174959698, 175303027, 
    175476244, 178255803
}

with open("13_chr_ver1.0.recode.vcf", "r") as vcf_file, open("output.txt", "w") as output_file:
    header_found = False  
    
    for line in vcf_file:
        if line.startswith('#'):
            if line.startswith("#CHROM"):
                output_file.write(line)
                header_found = True
            continue  
        if header_found:
            columns = line.split('\t')
            position = int(columns[1])
            if position in positions:
                output_file.write(line)
```
Запускаем скрипт:
```
python3 pos.py
```
Получаем отфильтрованный файл: https://github.com/erofeeva-plastilin/CMS/blob/main/output.txt  
Сводная таблица: https://github.com/erofeeva-plastilin/CMS/blob/main/final_table.txt  
Если обобщить, то у устойчивой линии 1/1, у неустойчивой линии - 0/0. Они очень четко разделились

## 4. Предсказание наших маркеров
Обновленный pos.py для фильтрации архивированных файлов (тут фильтрация по позициям сборки 2.0):

```
import gzip
positions = {
    153949714, 155465938, 155087415, 152602920, 152161520,
    156160357, 156370666, 156384352, 156456877, 155847945,
    159241056, 156940410, 153583253
}

with gzip.open("/mnt2/fastq2vcf_results/helianthus/2024-08-30_wgs_sunflower.vcf.gz", "rt", encoding="utf-8") as vcf_file, open("2024-08-30>
    header_found = False
    for line in vcf_file:
        if line.startswith('#'):
            if line.startswith("#CHROM"):
                output_file.write(line)
                header_found = True
            continue
        if header_found:
            columns = line.split('\t')
            position = int(columns[1])
            if position in positions:
                output_file.write(line)
```
Запускаем скрипт:
```
python3 pos.py
```
Получаем отфильтрованный файл: https://github.com/erofeeva-plastilin/CMS/blob/main/2024-08-30_wgs_sunflower_output.txt  
Финальная таблица: https://docs.google.com/spreadsheets/d/1B0p3dcaxF79B5V6aW1dlWlSZO0Mbgq-D9E_rFblMnaY/edit?gid=403362067#gid=403362067  
