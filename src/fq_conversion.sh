bwa mem -R ${4:-"@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1"} \
    -p $2 $1 > $3.sam
    

gatk SortSam \
   --INPUT $3.sam \
   --OUTPUT $3.bam \
   --SO coordinate