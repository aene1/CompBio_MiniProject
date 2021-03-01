import argparse
import os
from Bio import SeqIO
from Bio import Entrez

log_file = open('miniProject.log','w')



#download the SRR and fastq them
def download_data(SRR):
    wget = 'wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/' +  SRR + '/' + SRR + '.1'
    fastq = 'fastq-dump -I --split-files ' +SRR + '.1'
    rename = 'mv ' + SRR + '.1' +' ' + SRR
    os.system(wget)
    os.system(rename)
    os.system(fastq)


#extract the CDS in order to construct the transcriptome
def extract_CDS():
    Entrez.email = 'aene@luc.edu'
    handle = Entrez.efetch(db = 'nucleotide', id = 'EF999921', rettype='fasta')
    records = list(SeqIO.parse(handle,'fasta'))
    count = 0
    outfile = open('CDS_EF999921.fa','w') #the outfile where the CDS records will  be written out to
    for i in records.features:
        if i.type == 'CDS':
            outfile.write('>' + str(i.qualifiers['protein_id']) + '\n' + str(i.location.extract(i).seq)+'\n')
            count += 1
    log_file.write('The HCMV genome (EF999921) has ' + str(count) + '\n')








parser = argparse.ArgumentParser(description='Process some SRRs and split paired reads.')
parser.add_argument('SRR', metavar='N', type=str, nargs='+',
                    help='Compare SRR files')
# parser.add_argument('--download_SRR')
args = parser.parse_args()

for i in args.SRR:
    download_data(i)

extract_CDS()