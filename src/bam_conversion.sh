samtools index $1

gatk HaplotypeCaller  \
   -R $2 \
   -I $1 \
   -O $3.gz