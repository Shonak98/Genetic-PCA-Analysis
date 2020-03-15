plink2 \
  --vcf $1 \
  --make-bed \
  --snps-only \
  --maf $5 \
  --geno $6 \
  --mind $7 \
  --recode ${4:- }\
  --out $2/$3
