import argparse
import os
from Bio import SeqIO
from Bio import Entrez

log_file = open('miniProject.log','w')



def download_data(SRR):
    wget = 'wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/' +  SRR + '/' + SRR + '.1'
    fastq = 'fastq-dump -I --split-files' +SRR + '.1'
    os.system(wget)
    os.system(fastq)



parser = argparse.ArgumentParser(description='Process some SRRs and split paired reads.')
parser.add_argument('SRR', metavar='N', type=str, nargs='+',
                    help='Compare SRR files')
parser.add_argument('--download_SRR')

args = parser.parse_args()

for i in args.SRR:
    download_data(i)