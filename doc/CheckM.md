
+ only parse quality file

      BioSAK CheckM -i checkm_op.txt -o MAG_qualities.txt
      BioSAK CheckM -i checkm_op.txt -o MAG_qualities.txt -r

+ get the quality of qualified bins 

      BioSAK CheckM -i checkm_op.txt -cpl 99 -o MAG_qualities_cpl99.txt
      BioSAK CheckM -i checkm_op.txt -cpl 99 -ctm 5 -o MAG_qualities_cpl99_ctm5.txt

+ get the quality of qualified MAG and copy them into a separate folder

      BioSAK CheckM -i checkm_op.txt -g MAG_folder -x fa -cpl 99 -o MAG_qualities_cpl99.txt
      BioSAK CheckM -i checkm_op.txt -g MAG_folder -x fa -cpl 99 -ctm 5 -o MAG_qualities_cpl99_ctm5.txt

+ Note
adjusted_contamination = contamination x (100 - heterogeneity) / 100
