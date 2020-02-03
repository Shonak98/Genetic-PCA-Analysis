"""
Shonak Shah
DSC 180A Genetics
"""
#!C:/Users/Shonak/capstone_180a/Scripts/python

def data_pull(chromosomes, outpath):
    """
    Method goes through the genetics-config json file to determine
    which chromosomes to pull and the locations that I want to place them
    into. The files are downloaded to the same format as they were in the
    server and then unzipped to be usable for analysis.

    """
    import urllib.request as ureq
    import os

    if not os.path.exists(outpath):
        os.mkdir(outpath)


    #### shared by all the files
    directory = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    link_end = ".20130502.genotypes.vcf.gz"

    #### broken up based on which chromosome that is wanted
    numbered_chrom = ".phase3_shapeit2_mvncall_integrated_v5a"
    x_chrom = ".phase3_shapeit2_mvncall_integrated_v1b"
    y_chrom = ".phase3_integrated_v2a"

    for chrom in chromosomes:
        if chrom.lower() == "x":
            chrom_name = x_chrom
        elif chrom.lower() == "y":
            chrom_name = y_chrom
        else:
            chrom_name = numbered_chrom

        ureq.urlretrieve(
                         directory + f'ALL.chr{chrom}' + chrom_name + link_end,
                         f'{outpath}/Chr{chrom}.vcf.gz'
                         )
        # print("unzipping time")
        # unzip(f'{outpath}/Chr{chrom}.vcf.gz')


def unzip(file):
    """
    Method unzips the GZ files that are downloaded and creates an
    unzipped file in the same directory

    """
    import time
    import gzip
    import shutil

    with gzip.open(file, 'rb') as f_in:
        # name excludes last 3 letters (.gz) so it is now an unzipped version
        with open(file[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            # time.sleep(0.001)
    return
