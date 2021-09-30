from typing import List

from ncbi.datasets.openapi import ApiClient as DatasetsApiClient
from ncbi.datasets.openapi import ApiException as DatasetsApiException
from ncbi.datasets.openapi import GenomeApi as DatasetsGenomeApi

# 	SB0662_bin_1	Chloroflexi bacterium
accessions = ['GCA_009840555.1', 'GCA_009840575.1', 'GCA_009840525.1']
zipfile = '/Users/songweizhi/Desktop/SB0662_bin_1.zip'




def genome_download_by_asm_accessions(genome_assembly_accessions: List[str],
                                      zipfile_name: str):
    """Downloads an NCBI Datasets Genome Data Package

    Args:
        genome_assembly_accessions: A list of NCBI assembly accessions
        zipfile_name: output file name
    """
    with DatasetsApiClient() as api_client:
        genome_api = DatasetsGenomeApi(api_client)
        try:
            print('Begin download of genome data package ...')
            genome_ds_download = genome_api.download_assembly_package(
                genome_assembly_accessions,
                include_annotation_type=['RNA_FASTA', 'PROT_FASTA'],
                _preload_content=False)

            with open(zipfile_name, 'wb') as f:
                f.write(genome_ds_download.data)
            print(f'Download completed -- see {zipfile_name}')
        except DatasetsApiException as e:
            print(f'Exception when calling download_assembly_package: {e}\n')

#
# genome_download_by_asm_accessions(accessions, zipfile)



needed_id = set()
for each in open('/Users/songweizhi/Desktop/64genomes.txt'):
    needed_id.add(each.strip())
print(needed_id)

for each in open('/Users/songweizhi/Desktop/prokaryotes-2.csv'):

    assembly_id = each.strip().split(',')[5][1:-1]
    if assembly_id in needed_id:
        print(each.strip())
