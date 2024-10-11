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


