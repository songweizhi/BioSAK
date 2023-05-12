+ Example commands

      BioSAK sampling_GTDB_gnms -p Firm_o -meta bac120_metadata_r207.tsv -taxon p__Firmicutes -r o
      BioSAK sampling_GTDB_gnms -p Firm_f -meta bac120_metadata_r207.tsv -taxon p__Firmicutes -r f -cpl 85 -ctm 5 -rs
      BioSAK sampling_GTDB_gnms -p Firm_g -meta bac120_metadata_r207.tsv -taxon p__Firmicutes -r g -cpl 85 -ctm 5 -ts

+ -meta: bac120_metadata_r207.tsv or ar122_metadata_r207.tsv
+ Sampling rank need to be lower than the rank of specified taxon (-taxon).
+ Genome(s) with highest quality score(s)(defined as completeness-5*contamination) will be selected.
