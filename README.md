# panTro6 maf to vcf conversion



This repo contains the code to convert `hg19.panTro6.synNet.maf.gz` to an all sites vcf file where the only individual is panTro6 and is homozygous at every site. It should be noted that the code in this repo was inspired by Simon Martin's [ `genomics_general` repo](https://github.com/simonhmartin/genomics_general) that I repurposed for my own specific use. In this `README.md` I will outline the steps for the maf to vcf conversion.

## Step 1: Download `hg19.panTro6.synNet.maf.gz`

This can be done a lot of different ways, but here is how I did it.

```bash
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/vsPanTro6/hg19.panTro6.synNet.maf.gz
```



## Step 2: Run `chimp_maf_parser.py`

Once `hg19.panTro6.synNet.maf.gz` has been downloaded the next step is to extract the allele calls for panTro6, which can be accomplished with the following code.

```bash
for CHR in {1..22} X; do
module load tabix python/3.7.4; python3 maf_parser.py hg19.panTro6.synNet.maf.gz ${CHR} | bgzip > hg19_panTro6_chr${CHR}.txt.gz
done
```



## Step 3: Run `chimp_maf_output_to_vcf.py`

Now that all of the genotype information has been extracted for panTro6 the next step is to convert that information into a vcf file, which can be accomplished with the following code and the `generic_vcf_header.txt`.

```bash
for CHR in {1..22} X; do
module load tabix python/3.7.4; python3 chimp_maf_output_to_vcf.py hg19_panTro6_chr${CHR}.txt.gz | bgzip > hg19_panTro6_chr${CHR}.vcf.gz
done
```



## Step 4: Index the vcf files

Now that all panTro6 data has been converted to a vcf the last step is to create index files for those vcf files, which can be accomplished with the following code.

```bash
for VCF in *.vcf.gz; do
module load tabix; tabix -p vcf ${VCF}
done
```

