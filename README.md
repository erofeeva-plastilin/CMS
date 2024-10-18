# Важные команды
## 1. Фильтрация VCF по хромосомам и позициям
```
tabix -p vcf {input}.vcf.gz
bcftools view -R {pos}.txt {input}.vcf.gz -o {output}.vcf
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
