Assess the quality of marker genes for phylogeny inference with split score
---


## Reference

Please consider cite the following publication if you found this approach is helpful

+ Dombrowski, N., Williams, T.A., Sun, J. et al. Undinarchaeota illuminate DPANN phylogeny and the impact of gene transfer on archaeal evolution. Nat Commun 11, 3939 (2020). [https://doi.org/10.1038/s41467-020-17408-w](https://doi.org/10.1038/s41467-020-17408-w)

        
## Step one: infer gene tree

    TreeSAK SplitScore1 -i OrthologousGroups.txt -s OrthologousGroupsFasta -o step1_op_dir -t 6 -f
    TreeSAK SplitScore1 -i OrthologousGroups.txt -s OrthologousGroupsFasta -o step1_op_dir -t 6 -f -u interested_gnm.txt
    

## Step two: calculate split score

    # Please ensure that all the commands produced in step one have been executed before proceeding to step two.
    TreeSAK SplitScore2 -i step1_op_dir -g gnm_cluster.tsv -k gnm_taxon.txt -f -t 10 -o step_2_op_dir

