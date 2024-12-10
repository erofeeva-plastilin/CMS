# Важные команды
## 1. Фильтрация VCF по хромосомам и позициям
```
tabix -p vcf {input}.vcf.gz
bcftools view -R {pos}.txt {input}.vcf.gz -o {output}.vcf
```
**Посмотреть сколько SNP в файле**
```
bcftools view -v snps Pea_short.vcf | grep -vc "^#"
```
**Посмотреть сколько образцов в файле**
```
bcftools query -l public_to_Oryza_indica.vcf | wc -l
```
**Список образцов в файле**
```
bcftools query -l input.vcf
```
## 2. Подсчет частот аллелей
```
conda activate GWAS-PIPELINE
vcftools --vcf {input}.vcf --freq --out allele_frequencies
```
## 3. Переименовывание хромосом
```
# Создаем таблицу соответствий
chromosome_map = {
    'glyma.Wm82.gnm2.Gm01': 'Chr01',
    'glyma.Wm82.gnm2.Gm02': 'Chr02',
    'glyma.Wm82.gnm2.Gm03': 'Chr03',
    'glyma.Wm82.gnm2.Gm04': 'Chr04',
    'glyma.Wm82.gnm2.Gm05': 'Chr05',
    'glyma.Wm82.gnm2.Gm06': 'Chr06',
    'glyma.Wm82.gnm2.Gm07': 'Chr07',
    'glyma.Wm82.gnm2.Gm08': 'Chr08',
    'glyma.Wm82.gnm2.Gm09': 'Chr09',
    'glyma.Wm82.gnm2.Gm10': 'Chr10',
    'glyma.Wm82.gnm2.Gm11': 'Chr11',
    'glyma.Wm82.gnm2.Gm12': 'Chr12',
    'glyma.Wm82.gnm2.Gm13': 'Chr13',
    'glyma.Wm82.gnm2.Gm14': 'Chr14',
    'glyma.Wm82.gnm2.Gm15': 'Chr15',
    'glyma.Wm82.gnm2.Gm16': 'Chr16',
    'glyma.Wm82.gnm2.Gm17': 'Chr17',
    'glyma.Wm82.gnm2.Gm18': 'Chr18',
    'glyma.Wm82.gnm2.Gm19': 'Chr19',
    'glyma.Wm82.gnm2.Gm20': 'Chr20',
}

# Открываем исходный файл для чтения и создаем новый файл для записи
with open('yield.vcf', 'r') as vcf_in, open('yield1.vcf', 'w') as vcf_out:
    for line in vcf_in:
        # Если строка начинается с #, это метаданные, их не нужно изменять
        if line.startswith('#'):
            vcf_out.write(line)
        else:
            # Разделяем строку по табуляции, чтобы получить колонки
            columns = line.strip().split('\t')
            # Заменяем идентификатор хромосомы на основе таблицы
            if columns[0] in chromosome_map:
                columns[0] = chromosome_map[columns[0]]
            # Записываем строку в выходной файл
            vcf_out.write('\t'.join(columns) + '\n')
```
## 4. Билд базы для snpEff и аннотация

1. В папке /mnt/tools/snpEff/data создать папку для референса. Например soy_old. Туда нужно положить четыре файла -- геномную фасту (sequences.fa), белковую фасту (protein.fa), кодирующие сиквенсы (cds.fa) и аннотацию (genes.gff). Названия нужны именно такие
2. Дальше в файл /usr/snpEff.config нужно вписать <название созданной папки>.genome : <название созданной папки> (где-то в районе 137 строчки файла будет раздел с нашими геномами)
3. Запустить snpeff build <название созданной папки>. Если будет ругаться можно добавить опции -noCheckCds -noCheckProtein. Через sudo
4. Аннотация ```snpeff -v <название созданной папки> <путь к нужному vcf> > <название нового vcf>```

## 5. Скачивание SRR файлов
```
/mnt/tools/sratoolkit.3.0.6-ubuntu64/bin/prefetch SRR15130914 --max-size 50G
/mnt/tools/sratoolkit.3.0.6-ubuntu64/bin/fastq-dump SRR15130914 (--split-files)
```
## 6. Создание png с распределением частот аллелей
```
import matplotlib.pyplot as plt
import pandas as pd
file_path = '22112024final.tsv'  
columns = ["CHROM", "POS", "N_ALLELES", "N_CHR", "REF_FREQ", "ALT_FREQ"]
data = pd.read_csv(file_path, sep='\t', names=columns, skiprows=1)
data = data.drop(columns=["N_ALLELES", "N_CHR"])
data[["REF", "FREQ_REF"]] = data["REF_FREQ"].str.split(":", expand=True)
data[["ALT", "FREQ_ALT"]] = data["ALT_FREQ"].str.split(":", expand=True)
data["FREQ_REF"] = pd.to_numeric(data["FREQ_REF"])
data["FREQ_ALT"] = pd.to_numeric(data["FREQ_ALT"])
data = data.drop(columns=["REF_FREQ", "ALT_FREQ"])
print(data)
plt.figure(figsize=(10, 6))
plt.hist(data["FREQ_REF"], bins=500, alpha=0.5, label='FREQ_REF')
plt.hist(data["FREQ_ALT"], bins=500, alpha=0.5, label='FREQ_ALT')
plt.xlabel('Frequency')
plt.ylabel('Count')
plt.title('Histogram of FREQ_REF and FREQ_ALT')
plt.legend()
plt.grid(axis='y')
output_image_path = 'histogram_freq_ref_alt.png'
plt.savefig(output_image_path, dpi=300)
plt.show()
```
## 7. Создание vcf с рандомными образцами в определенном количестве
```
import random
import subprocess

# Исходный файл
vcf_file = "rename_chr_our.vcf"
sample_list_file = "samples.txt"

# Список размеров выборок
sample_sizes = [10, 20, 30, 40, 50, 75, 100]
subprocess.run(f"bcftools query -l {vcf_file} > {sample_list_file}", shell=True)
with open(sample_list_file, "r") as file:
    samples = [line.strip() for line in file]
for size in sample_sizes:
    selected_samples = random.sample(samples, size)
    temp_sample_file = f"selected_samples_{size}.txt"
    with open(temp_sample_file, "w") as temp_file:
        temp_file.write("\n".join(selected_samples))
    output_vcf = f"rename_chr_our_{size}.vcf"
    subprocess.run(
        f"bcftools view -S {temp_sample_file} -o {output_vcf} -O v {vcf_file}",
        shell=True
    )
    print(f"VCF with {size} samples saved as {output_vcf}")
subprocess.run("rm selected_samples_*.txt", shell=True)
```
## 8. Добавить в гит большой файл
Скачиваем/обновляем локально наш репозиторий.
```
git lfs install # Инициализация Git LFS
git lfs track "CoordTransfer/Chain_files/Pisum_sativum_v1a.fa.to.GCF_024323335.1_CAAS_Psat_ZW6_1.0_genomic.unmasked.fna.over.chain" # Добавляем файл отслеживания
cd /mnt/users/erofeevan/Converting_Genome_Coordinates_pipeline/CoordTransfer # Перемещаемся в корневую директорию репозитория и все оставшиеся команды там делаем!
git add CoordTransfer/Chain_files/Pisum_sativum_v1a.fa.to.GCF_024323335.1_CAAS_Psat_ZW6_1.0_genomic.unmasked.fna.over.chain # Добавляем файл
git status # Проверяем статус
git commit -m "Add large chain file for Pisum sativum genome" # Делаем коммит
git push origin main # Отправляем изменения
```






