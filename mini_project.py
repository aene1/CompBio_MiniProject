import argparse
import os
from Bio import SeqIO
from Bio import Entrez
import glob

log_file = open('miniProject.log','w+')
path = os.getcwd()

#speed up spades run

#download the SRR and fastq them
def download_data(SRR):
    wget = 'wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/' + SRR + '/' + SRR + '.1'
    rename = 'mv ' + SRR + '.1' + SRR
    os.system(wget)
    os.system(rename)



def fastq(SRR):
    fastq = 'fastq-dump -I --split-files ' +SRR + '.1'
    rename = 'mv ' + SRR + '.1' + SRR
    os.system(fastq)
    os.system(rename)
    print(fastq)


#extract the CDS in order to construct the transcriptome
def extract_CDS():
    Entrez.email = 'aene@luc.edu'
    handle = Entrez.efetch(db = 'nucleotide', id = 'EF999921', rettype='gb',retmode = 'text')
    records = SeqIO.parse(handle,'genbank')
    count = 0
    outfile = open('CDS_EF999921.fa','w') #the outfile where the CDS records will  be written out to
    for i in records:
        for j in i.features:
            if j.type == 'CDS':
                outfile.write('>' + str(j.qualifiers['protein_id']).replace("'",'').replace('[','').replace(']','') + '\n' + str(j.location.extract(i).seq)+'\n')
                count += 1
    log_file.write('The HCMV genome (EF999921) has ' + str(count) + '\n')


#runn kallisto based on SRR numbers and the SRR files dowloaded before
def kallisto(SRR):
    kallisto_index = 'time kallisto index -i HCMV_index.idx CDS_EF999921.fa'
    os.system(kallisto_index)
    run_kallisto = 'time kallisto quant -i HCMV_index.idx -o' + path + '/results_' + SRR +' -b 30 -t 4 '+ SRR + '.1_1.fastq '+ SRR+ '.1_2.fastq'
    os.system(run_kallisto)



def SleuthInput(SRR):
    #output file that goes in R
    output = open('input_sleuth.txt', 'w+')
    # initial line in file
    output.write('sample' + '\t' + 'condition' + '\t' + 'path' + '\n')
    # based on SRR number, write condition and path to output file
    for i in SRR:
        path = i
        if int(i[3:]) % 2 == 0:  #if it is even then it is condition 1 as in 2dpi, if it is not then it is condition 2 as in 6dpi
            output.write(str(i) + '\t' + '2dpi' + '\t' +  path + '\n')
        else:
            output.write(str(i) + '\t' + '6dpi' + '\t' + path + '\n')
    output.close()

def Sleuth():
    Sleuth_command = 'Rscript Sleuth.R'
    os.system(Sleuth_command)
    read_sleuth = open('R_sleuth_output.txt','r').readlines()
    for i in read_sleuth:
        log_file.write(i + 'n')


#build bowtie2 and run it to get SAM files

def index_bowtie():
    handle = Entrez.efetch(db='nucleotide', id='EF999921', rettype='gb', retmode='text')
    records = list(SeqIO.parse(handle, 'fasta'))
    for i in records:
        print(i)
        SeqIO.write(i, 'EF999921.fa', "fasta")


def bowtie2(SRR):
    build_bowtie2 = 'bowtie2-build EF999921.fa EF99992_1'
    os.system(build_bowtie2)
    bowtie2 = 'bowtie2 --quiet --no-unal --al-conc EF999921_' + SRR + '.fastq -x EF999921_1 -1 '+ SRR+ '.1_1.fastq -2 ' + SRR+ '.1_2.fastq -S EF999921_' + SRR+ '.sam'
    os.system(bowtie2)



parser = argparse.ArgumentParser(description='Process some SRRs and split paired reads.')
parser.add_argument('SRR', metavar='N', type=str, nargs='+',
                    help='SRR files we want for the comparasion')
parser.add_argument('--download_files', metavar='N', type=str, nargs='+',
                    help='Download the SRR, instead of using the server ones')
args = parser.parse_args()


in_path = os.getcwd()
files = glob.glob(("**/*"), recursive=True)
files = [f for f in files if os.path.isfile(f)]

for i in args.SRR:
    if i not in files:
        download_data(i)

extract_CDS()
index_bowtie()
for i in args.SRR:
    fastq(i)
    kallisto(i)
    bowtie2(i)

SleuthInput(args.SRR)
Sleuth()




