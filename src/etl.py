
#### imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json
import os
import urllib.request as ureq
import gzip
import shutil
import time


############### param config functions

def get_chrom_list(chr_lst):
    chromes = []

    for ind, curr in enumerate(chr_lst):
        if curr == '-':
            chromes += list(range(chr_lst[ind-1]+1, chr_lst[ind+1]))
        else:
            chromes.append(curr)
    arr = np.array(chromes)

    return np.unique(arr[(arr>0) & (arr<23)])



############## directory set up functions

def get_files_in_directory(directory):
    files = []
    locations_file = "locations.txt"
    if os.path.isdir(directory):
        os.system(f"cd {directory}; dir>{locations_file}")

        with open(directory + f'/{locations_file}') as fh:
            for i in fh:
                files += i.split()

        os.system(f"rm -rf {directory}/{locations_file}")
    return files


def get_fasta_location(directory):
    files = get_files_in_directory(directory)
    location = ''
    for j in files:
        if j[-6:] == '.fasta':
            location = j
    return location


def validate_dict_file(directory):
    files = get_files_in_directory(directory)
    fasta = get_fasta_location(directory)
    if fasta != '':
        file_name = fasta[:-6]
        if file_name +'.dict' in files:
            return True
    return False


def create_references(location, references):
    if location[-1] != '/':
        location += '/'
    if os.path.isdir(references) == False:
        os.system(f'mkdir {references}')

    fasta = get_fasta_location(references)

    if (fasta == '') and (validate_dict_file(references) == False):
        os.system(f"cp -r {location}* {references}")


    fasta = get_fasta_location(references)
    os.system(f"samtools faidx {references}/{fasta}")

    if validate_dict_file(references) == False:
        dict_file = fasta[:-6] + ".dict"
        os.system(f"gatk CreateSequenceDictionary --REFERENCE '{references}/{fasta}' \
                   --OUTPUT '{references}/{dict_file}'")



################## data pulling

def data_pull(chromosomes_paths, output_name, output_dir = None):
    """
    Method goes through the genetics-config json file to determine
    which chromosomes to pull and the locations that I want to place them
    into. The files are downloaded to the same format as they were in the
    server and then unzipped to be usable for analysis.
    """

    if output_dir != None:
        if os.path.exists(output_dir) == False:
            os.mkdir(output_dir)

    for chrom in chromosomes_paths:
        if output_dir != None:
            ureq.urlretrieve(chrom, f'{output_dir}/{output_name}')
        else:
            ureq.urlretrieve(chrom, f'{output_name}')


def unzip(file):
    """
    Method unzips the GZ files that are downloaded and creates an
    unzipped file in the same directory
    """

    with gzip.open(file, 'rb') as f_in:
        # name excludes last 3 letters (.gz) so it is now an unzipped version
        with open(file[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            time.sleep(0.01)
    return



################## data conversion functions

def fastq_to_bam(fastq, fasta, directory, read_group = "default", result_name = None, keep_sam = True):
    if os.path.isdir(directory) == False:
        os.system(f"mkdir {directory}")

    no_ext = fastq[-fastq[::-1].find('/'):-3]
    if result_name != None:
        no_ext = str(result_name)

    if no_ext[-4:] == ".bam":
        no_ext = no_ext[:-4]

    os.system("chmod 700 src/fq_conversion.sh")
    if read_group == "default":
        os.system(f"./src/fq_conversion.sh {fastq} {fasta} {directory}/{no_ext}")
    else:
        os.system(f"./src/fq_conversion.sh {fastq} {fasta} {directory}/{no_ext} {read_group}")

    if not keep_sam:
        os.system(f"rm -rf {no_ext}.sam")


def bam_to_vcf(bam, fasta, directory, output_name = None):
    if os.path.isdir(directory) == False:
        os.system(f"mkdir {directory}")

    if '/' in bam:
        no_ext = bam[-bam[::-1].find('/'):-4]
    else:
        no_ext = bam[:-4]

    if output_name != None:
        no_ext = str(output_name)

    os.system("chmod 700 src/bam_conversion.sh")
    os.system(f"./src/bam_conversion.sh {directory}/{bam} {fasta} {directory}/{no_ext}")



################## data manipulations

def filtering_vcf(vcf, output_name, directory, maf, geno, mind, output_type = '', keep_files = True):
    if os.path.isdir(directory) == False:
        os.system(f"mkdir {directory}")

    os.system("chmod 700 src/filtering.sh")

    if output_type != '':
        os.system(f"./src/filtering.sh {vcf} {directory} {output_name} {maf} {geno} {mind} {output_type}")
    else:
        os.system(f"./src/filtering.sh {vcf} {directory} {output_name} {maf} {geno} {mind}")

    if "vcf" in output_type:
        os.system(f"bgzip -f {directory}/{output_name}.vcf > {directory}/{output_name}.vcf.gz")

    if keep_files == False:
        os.system(f'rm -rf {directory}/{output_name}.bed')
        os.system(f'rm -rf {directory}/{output_name}.bim')
        os.system(f'rm -rf {directory}/{output_name}.bam')
        os.system(f'rm -rf {directory}/{output_name}.nosex')


def vcf_concat(vcf_lst, output, directory):
    string_lst = "bcftools concat "
    for vcf in vcf_lst:
        string_lst += f"{directory}/{vcf} "

    string_lst += f"-o {directory}/{output}"
    os.system(string_lst)


def pca(file, file_type, output, directory, pca_num = 2):
    if file_type == 'vcf':
        os.system(f"plink2 --vcf {directory}/{file} --pca {pca_num} --out {directory}/{output}")
    elif file_type == 'bfile':
        os.system(f"plink2 --bfile {directory}/{file} --pca {pca_num} --out {directory}/{output}")


################## analysis functions

def plot_clusters(eig_file, chrom_lst, picture_name = None, codes_file = "data/igsr_samples.tsv",
                  with_colors = True):
    """
    Plots clusters from filtered genetic data

    @param eig_file : 4 column eigenvec file created from plink2 pca
    @param picture_name: (optional) file name that cluster plot will be saved to
    @param pop_codes : (optional) igsr_samples.tsv default (comes from 1000 Genomes page)
    @param with_color : (optional) Boolean that determines if resulting plot should have colors
    """
    eigen_df = pd.read_csv(eig_file, sep = ' ', names = ['Sample name', 'Population', 'x', 'y'])
    pop_codes = pd.read_csv(codes_file, sep = '\t')[['Sample name', 'Superpopulation code']]

    super_code = eigen_df.join(pop_codes.set_index('Sample name'), on = 'Sample name')

    if with_colors:
        sns.scatterplot(x = super_code.x, y = super_code.y, hue = super_code['Superpopulation code'])
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    else:
        sns.scatterplot(x = super_code.x, y = super_code.y, edgecolor = None, s= 5)

    plt.title(f"PCA Analysis on Chromosomes {', '.join(map(str,chrom_lst))}")
    plt.xlabel("PCA Component 1")
    plt.ylabel("PCA Component 2")

    if picture_name != None:
        plt.savefig(picture_name, bbox_inches='tight')   
