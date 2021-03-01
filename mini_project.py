import os
from Bio import SeqIO
from Bio import Entrez

log_file = open('miniProject.log','w')

def arg_parser():
    #sets up the argument parser and returns it so it can run forward"
    parser = argparse.ArgumentParser(description = "Input sra files")
    parser.add_argument(parser.add_argument('download_SRR_files', metavar = 'N/A', type=str, nargs = '*', help = 'SRR download')
    return parser.parse_args()

if args.donload_SRR_files != "N/A":
    for i in args.do:

    #finish this and change the parser

def dowload_data(SRR):
	 os.system('wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/' +  SRR + '/' + SRR + '.1')
         os.system('fastq-dump -I --split-files' + SRR + '.1')



