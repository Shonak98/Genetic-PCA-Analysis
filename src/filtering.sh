plink2 \
  --vcf $1 \
  --make-bed \
  --snps-only \
  --maf $4 \
  --geno $5 \
  --mind $6 \
  --recode ${7:- }\
  --out $2/$3
