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
sed -i 's/^{old_name}\t/{new_name}\t/' {input}.vcf > {output}.vcf
```
## 4. Билд базы для snpEff

1. В папке /mnt/tools/snpEff/data создать папку для референса. Например soy_old. Туда нужно положить четыре файла -- геномную фасту (sequences.fa), белковую фасту (protein.fa), кодирующие сиквенсы (cds.fa) и аннотацию (genes.gff). Названия нужны именно такие
2. Дальше в файл /usr/snpEff.config нужно вписать <название созданной папки>.genome : <название созданной папки> (где-то в районе 137 строчки файла будет раздел с нашими геномами)
3. Запустить snpeff build <название созданной папки>. Если будет ругаться можно добавить опции -noCheckCds -noCheckProtein. Через sudo

