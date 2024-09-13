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
Активация среды, где есть plink:
```
source ~/.bashrc
conda activate GWAS-PIPELINE
```
Фильтрация с размером окна 100кб, шагом 5 и r^2 0,5
Добавление уникального ID в VCF файл, чтобы избежать ошибок
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

## Gapit
Активация среды, где есть gapit:
```
source ~/.bashrc
conda activate /mnt/users/grigorieval/miniconda3/envs/gapit
```
Запуск gapit 
```
Rscript ./run_GAPIT.R pruned.vcf phenotype.CMS.tsv /mnt/users/erofeevan/Drafts/Drafts
```
Результаты gapit:
![image](https://github.com/user-attachments/assets/4cfb54de-8016-4e13-92d6-4359753dc511)
Выделены маркеры на 10 и 13 хромосомах, что соответствует статьям.
Отфильтрованные значимые маркеры:
	SNP	Chr	Pos	P.value	MAF	traits
5780636	NC_007977.1:76370AC	NC_007977.1	76370	7.34E-11	0.4594595	MLM.CMS
5780644	NC_023337.1:169028GT	NC_023337.1	169028	7.34E-11	0.4594595	MLM.CMS
121005	NC_035433.2:74500460AC	NC_035433.2	74500460	7.34E-11	0.4594595	MLM.CMS
123467	NC_035433.2:76166671GT	NC_035433.2	76166671	7.34E-11	0.4594595	MLM.CMS
372073	NC_035434.2:72733198AG	NC_035434.2	72733198	7.34E-11	0.4594595	MLM.CMS
539361	NC_035434.2:164133021TC	NC_035434.2	164133021	1.46E-09	0.2567568	MLM.CMS
878762	NC_035436.2:901328AC	NC_035436.2	901328	7.34E-11	0.4594595	MLM.CMS
893729	NC_035436.2:14115668TG	NC_035436.2	14115668	6.39E-09	0.4662162	MLM.CMS
1295825	NC_035436.2:189985344CT	NC_035436.2	189985344	7.34E-11	0.4594595	MLM.CMS
1868205	NC_035438.2:118764636TC	NC_035438.2	118764636	1.35E-10	0.4527027	MLM.CMS
2279874	NC_035440.2:9506603AG	NC_035440.2	9506603	1.35E-10	0.4527027	MLM.CMS
2321915	NC_035440.2:42876567GA	NC_035440.2	42876567	5.11E-10	0.4459459	MLM.CMS
2966298	NC_035442.2:42721703TC	NC_035442.2	42721703	7.18E-09	0.4864865	MLM.CMS
2966541	NC_035442.2:42870270CT	NC_035442.2	42870270	5.35E-09	0.4797297	MLM.CMS
2966589	NC_035442.2:42916873AT	NC_035442.2	42916873	5.07E-09	0.4662162	MLM.CMS
2967165	NC_035442.2:43547484GA	NC_035442.2	43547484	3.27E-09	0.472973	MLM.CMS
2967972	NC_035442.2:44022537GA	NC_035442.2	44022537	4.75E-09	0.4121622	MLM.CMS
4262659	NC_035445.2:152161520TA	NC_035445.2	152161520	3.96E-09	0.4797297	MLM.CMS
4263452	NC_035445.2:152602920TG	NC_035445.2	152602920	8.63E-09	0.4864865	MLM.CMS
4264507	NC_035445.2:153583253CT	NC_035445.2	153583253	3.29E-09	0.4594595	MLM.CMS
4264662	NC_035445.2:153949714CA	NC_035445.2	153949714	3.46E-09	0.4527027	MLM.CMS
4266460	NC_035445.2:155087415GA	NC_035445.2	155087415	7.81E-09	0.4594595	MLM.CMS
4266667	NC_035445.2:155465938CT	NC_035445.2	155465938	6.00E-09	0.4594595	MLM.CMS
4267041	NC_035445.2:155847945CT	NC_035445.2	155847945	1.36E-09	0.4797297	MLM.CMS
4267317	NC_035445.2:156160357GT	NC_035445.2	156160357	6.42E-09	0.4932432	MLM.CMS
4267531	NC_035445.2:156370666GT	NC_035445.2	156370666	4.33E-11	0.4459459	MLM.CMS
4267546	NC_035445.2:156384352AG	NC_035445.2	156384352	4.28E-09	0.4594595	MLM.CMS
4267906	NC_035445.2:156456877AG	NC_035445.2	156456877	5.61E-09	0.4121622	MLM.CMS
4268250	NC_035445.2:156940410GA	NC_035445.2	156940410	8.60E-09	0.5	MLM.CMS
4269518	NC_035445.2:159241056GT	NC_035445.2	159241056	2.10E-09	0.4594595	MLM.CMS
4285095	NC_035446.2:1455979AC	NC_035446.2	1455979	7.34E-11	0.4594595	MLM.CMS
4376394	NC_035446.2:52076210AG	NC_035446.2	52076210	4.81E-10	0.4459459	MLM.CMS
5271507	NC_035448.2:116333387TG	NC_035448.2	116333387	5.84E-10	0.472973	MLM.CMS
5623926	NC_035449.2:102286735TC	NC_035449.2	102286735	6.68E-10	0.4459459	MLM.CMS
![image](https://github.com/user-attachments/assets/e014bc80-c6b2-49c6-94e6-948e92d845fa)


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
