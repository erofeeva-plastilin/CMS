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
## 2. Подсчет частот аллелей
```
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
/mnt/tools/sratoolkit.3.0.6-ubuntu64/bin/fastq-dump SRR15130914
```
## 6. Создание png с распределением частот аллелей
```
import pandas as pd
import matplotlib.pyplot as plt
file_path = "our_sunflower.tsv"
columns = ["CHROM", "POS", "N_ALLELES", "N_CHR", "ALLELE_FREQ"]
allele_freq = pd.read_csv(file_path, sep="\t", names=columns, skiprows=1, on_bad_lines="skip")

freqs = []
for row in allele_freq["ALLELE_FREQ"]:
    if isinstance(row, str):  
        #row = row.strip() 
        for allele in row.split(" "):  
            try:
                freqs.append(float(allele.split(":")[1]))  
            except (IndexError, ValueError):
                continue

plt.figure(figsize=(8, 6))
plt.hist(freqs, bins=50, alpha=0.7, color="blue", edgecolor="black")
plt.title("Distribution of allele frequencies in our_sunflower VCF")
plt.xlabel("Allele frequency")
plt.ylabel("Number")
plt.grid(axis="y", linestyle="--", alpha=0.7)
output_file = "our_sunflower.png"
plt.savefig(output_file, dpi=300, bbox_inches="tight")
#plt.show()
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
