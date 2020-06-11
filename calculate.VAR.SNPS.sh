echo "Input species list: "
read SPECIES


bcftools view -a -S $SPECIES combined.genotyped.SNP.filtered.PASS.NOMISSING.vcf.gz | \
	gzip > combined.genotyped.SNP.filtered.PASS.NOMISSING.$SPECIES.vcf.gz
echo "Number of SNPs in total vcf file: "
gzip -cd combined.genotyped.SNP.filtered.PASS.NOMISSING.$SPECIES.vcf.gz | grep -vc '^#'

bcftools view -m 2 -M 2 -i 'ALT!="."' combined.genotyped.SNP.filtered.PASS.NOMISSING.$SPECIES.vcf.gz | \
	gzip > combined.genotyped.SNP.filtered.PASS.NOMISSING.$SPECIES.VAR.vcf.gz
echo "Number of SNPs only in x species: "
gzip -cd combined.genotyped.SNP.filtered.PASS.NOMISSING.$SPECIES.VAR.vcf.gz | grep -vc '^#'


echo $SPECIES >> SNP.count.txt
gzip -cd combined.genotyped.SNP.filtered.PASS.NOMISSING.$SPECIES.VAR.vcf.gz | grep -vc '^#' >> SNP.count.txt

