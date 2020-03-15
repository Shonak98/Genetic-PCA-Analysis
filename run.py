import sys
sys.path.insert(0, f"{sys.path[0]}/src")

import json
from etl import *

def load(params):
    links = params['chromosome_links']
    outpath = params['outpath']

    for link in links:
        file_name = link[-link[::-1].find('/'):]
        data_pull(link, file_name, outpath)
        unzip(f"{outpath}/{file_name}")


def project(params):

    ######## create reference files
    create_references(params["references"], params["new_ref_location"])

    ######## conversions
    fasta = params["new_ref_location"] + '/' + get_fasta_location(params["new_ref_location"])
    fastq_to_bam(params['test_fastq'], fasta, params['conversion_dir'], result_name = params['bam_name'])
    bam_to_vcf(params['bam_name'], fasta, params['conversion_dir'], output_name = params['vcf_name'])

    ######## vcf analysis
    chrom_lst = get_chrom_list(params['chromosomes'])
    directory = params['directory']
    extension = params['extension']
    vcf_dir = params['vcf_dir']

    vcf_lst = [f'{directory}{chrom}{extension}' for chrom in chrom_lst]

    for ind, vcf in enumerate(vcf_lst):
        filtering_vcf(vcf, chrom_lst[ind], vcf_dir, params['maf'],
                   params['geno'], params['mind'], param['fil_output'],
                   keep_files = params['keep_filtered'])

    zipped_lst = [f"{chrom}.vcf.gz" for chrom in chrom_lst]
    vcf_concat(zipped_lst, params['merged_vcf'], vcf_dir)

    pca(params['merged_vcf'], param['fil_output'], params['eigen_file'], vcf_dir)

    eigen_file = vcf_dir + '/' + params['eigen_file'] + '.eigenvec'
    plot_clusters(eigen_file, chrom_lst, picture_name = params['picture_name'])


if __name__ ==  '__main__':
    arguments = sys.argv
    params = json.load(open("config/params.json"))

    if 'load' in arguments:
        load(params)

    if 'test-project' in arguments:
        project(params)
