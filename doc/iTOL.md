
### Dataset types

+ ColorStrip https://itol.embl.de/help.cgi#strip

      BioSAK iTOL -ColorStrip -lg MagTaxon.txt -lt Phylum -o ColorStrip_taxon.txt

+ ColorRange https://itol.embl.de/help.cgi#dsRanges

      BioSAK iTOL -ColorRange -lg MagTaxon.txt -lt Phylum -o ColorRange_taxon.txt
      BioSAK iTOL -ColorRange -taxon Taxonomy.txt -rank f -lt Family -o ColorRange_taxon.txt

+ SimpleBar https://itol.embl.de/help.cgi#bar

      BioSAK iTOL -SimpleBar -lv MagSize.txt -scale 0-3-6-9 -lt Size -o SimpleBar_size.txt

+ Heatmap https://itol.embl.de/help.cgi#heatmap

      BioSAK iTOL -Heatmap -lm MagAbundance.txt -lt Abundance -o Heatmap_abundance.txt

+ Binary https://itol.embl.de/help.cgi#binary

      BioSAK iTOL -Binary -lm Binary_matrix.txt -lt Presence_Absence -gc lightblue -o PA_iTOL.txt
      BioSAK iTOL -Binary -lm Binary_matrix.txt -lt Presence_Absence -gc "#85C1E9" -o PA_iTOL.txt

+ ExternalShape https://itol.embl.de/help.cgi#shapes

      BioSAK iTOL -ExternalShape -lm identity_matrix.txt -lt Identity -scale 25-50-75-100 -o ExternalShape_identity.txt


### Input file format

+ Leaf-to-Group (`-lg`, tab separated, no header)

      genome_1  Bacteria
      genome_2  Archaea

+ Taxonomy (`-taxon`, tab separated, GTDB-format taxononomy string)

      genome_1	d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Dongiales;f__Dongiaceae;g__Dongia;s__Dongia mobilis
      genome_2	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Arenicellales;f__LS-SOB;g__VYGS01;s__

+ Group-to-Color (`-gc`) and Column-to-Color (`-cc`) (tab separated, no header)

      Bacteria    #CCCC00
      Archaea #9999FF
      Virus   orange
      Eukaryote   lightblue

  -Please note that only one color can be specified for Binary data, provide with `-gc lightblue` or `-gc "#85C1E9"`


+ Leaf-to-Value file format (`-lv`, tab separated, no header)

      genome_1	6.15
      genome_2	6.63

+ Leaf-to-Matrix file format (`-lm`, tab separated, header required!!!)

      Genome_id Sample_A   Sample_B   Sample_C
      genome_1	6.15    2.23    1.56
      genome_2	6.63    1.72    2.55


## Here is a tutorial

1. Here, I have a phylogenetic tree for 37 MAGs derived from six microbial communities (either surface-associated or planktonic).
I have the taxonomy info of these MAGs at the class level, their sizes and their relative abundance across the six samples.
This short note shows how to visualize all these info in one plot.

1. Download files from [demo_data/iTOL](../demo_data/iTOL)

1. File description (Please see below for how to prepare these files.):

    + A phylogenetic tree in [Newick](http://evolution.genetics.washington.edu/phylip/newicktree.html) format: [NorthSea_0_Tree.newick](files_needed/NorthSea_0_Tree.newick)
    + Taxonomy info: [NorthSea_1_Taxon_ColorStrip.txt](files_needed/NorthSea_1_Taxon_ColorStrip.txt), [NorthSea_1_Taxon_Range.txt](files_needed/NorthSea_1_Taxon_Range.txt)
    + Life-style info: [NorthSea_2_LifeStyle.txt](files_needed/NorthSea_2_LifeStyle.txt)
    + Abundance across samples: [NorthSea_3_Abundance.txt](files_needed/NorthSea_3_Abundance.txt)
    + MAG size info: [NorthSea_4_MAG_Size.txt](files_needed/NorthSea_4_MAG_Size.txt)

1. Upload **NorthSea_0_Tree.newick** to iTOL via [https://itol.embl.de/upload.cgi](https://itol.embl.de/upload.cgi).

1. Once you have the tree uploaded, you'll see the skeleton of the tree without any decoration. 
You can now play around with the control panel on the right side (e.g. change tree layout to circular).
![Step_1](figures/Step_1.jpg)

1. We are going to add the taxonomy info of our MAG to the tree now, which is really easy to do in iTOL. 
You just need to drag and drop **NorthSea_1_Taxon_ColorStrip.txt**  to the **tree area**.

1. Do the same thing to **NorthSea_2_LifeStyle.txt**, **NorthSea_3_Abundance.txt** and **NorthSea_4_MAG_Size.txt** to add life-style, abundance and size info, 
you'll see trees like this:
![Tree_1](figures/Tree_1.jpg)

1. To get a tree with a circular layout and MAG classes coloured as in the figure below. 
You need to use **NorthSea_1_Taxon_Range.txt** instead of **NorthSea_1_Taxon_ColorStrip.txt**,
choose "**Circular**" mode in the control panel, click "**At tips**" and then turn it **off**.
![Tree_2](figures/Tree_2.jpg)

1. Go to the **Export** panel, choose desired file format and export your tree to file. 
Remember to turn on **Colored ranges legend**, if you are using **NorthSea_1_Taxon_Range.txt** to color MAG classes.
![Step_2](figures/Step_2.jpg)


## How to prepare iTOL-recognizable files

1. Some tools suggested by iTOL: https://itol.embl.de/help.cgi#external

2. You can also use [BioSAK](https://github.com/songweizhi/BioSAK)'s iTOL module to generate iTOL-recognizable files. Please refers to the "Installation" section [here](https://github.com/songweizhi/BioSAK) for its installation.
    
  The purpose for developing this module is to generate iTOL-recognizable file for your dataset, 
  parameters (e.g., colour, font size...) provided in the generated file might need 
  further adjustment.
  
  Input files for the following commands can be found in `files_needed/raw_data`.
     
      BioSAK iTOL -ColorRange -lg raw_MAG_taxon.txt -lt Class -out NorthSea_1_Taxon_Range.txt
      BioSAK iTOL -ColorStrip -lg raw_MAG_taxon.txt -lt Class -out NorthSea_1_Taxon_ColorStrip.txt
      BioSAK iTOL -ColorStrip -lg raw_MAG_LifeStyle.txt -lt LifeStyle -out NorthSea_2_LifeStyle.txt
      BioSAK iTOL -Heatmap -lm raw_MAG_abundance.txt -lt Abundance -out NorthSea_3_Abundance.txt
      BioSAK iTOL -SimpleBar -lv raw_MAG_size.txt -scale 0-3-6-9 -lt MAG_Size -out NorthSea_4_MAG_Size.txt
      
      # for help
      BioSAK iTOL -h


## Help information

1. More examples: [https://itol.embl.de/help.cgi](https://itol.embl.de/help.cgi)
1. The Newick tree format: [http://evolution.genetics.washington.edu/phylip/newicktree.html](http://evolution.genetics.washington.edu/phylip/newicktree.html)
1. Hex Color Codes: [https://htmlcolorcodes.com](https://htmlcolorcodes.com) and [https://www.color-hex.com](https://www.color-hex.com)
