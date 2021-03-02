import argparse
import os
from Bio import SeqIO
from Bio import Entrez

log_file = open('miniProject.log','w+')


#download the SRR and fastq them
def download_data(SRR):
    wget = 'wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/' + SRR + '/' + SRR + '.1'
    fastq = 'fastq-dump -I --split-files ' +SRR + '.1'
    rename = 'mv ' + SRR + '.1' +' ' + SRR
    os.system(wget)
    os.system(rename)
    os.system(fastq)


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
    kallisto_index = 'time kallisto index -i HCMVindex.idx CDS_EF999921.fa'
    os.system(kallisto_index)
    run_kallisto = 'time kallisto quant -i HCMVindex.idx -o ./' + SRR +' -b 30 -t 4 '+ SRR + '_1.fastq '+ SRR+ '_2.fastq'
    os.system(run_kallisto)



def SleuthInput(SRRs):
    #output file that goes in R
    output = open('input_sleut.txt', 'w')
    # initial line in file
    output.write('sample' + '\t' + 'condition' + '\t' + 'path' + '\n')
    # based on SRR number, write condition and path to output file
    for i in SRRs:
        if int(i[3:]) % 2 == 0:  #if it is even then it is condition 1 as in 2dpi, if it is not then it is condition 2 as in 6dpi
            output.writeln(str(i) + '\t' + "2dpi" + '\t')
        else:
            output.writeln(str(i) + '\t' + "6dpi" + '\t')
    output.close()

def Sleuth():
    Sleuth_command = 'R Sleuth.R'
    os.system(Sleuth_command)
    sleuth_output = 'R_sleuth_output.txt'
    read_sleuth = open(sleuth_output).readlines()
    for i in read_sleuth:
        log_file.write(i + 'n')





parser = argparse.ArgumentParser(description='Process some SRRs and split paired reads.')
parser.add_argument('SRR', metavar='N', type=str, nargs='+',
                    help='SRR files we want for the comparasion')
parser.add_argument('--download_files', metavar='N', type=str, nargs='+',
                    help='Download the SRR, instead of using the server ones')
args = parser.parse_args()
extract_CDS()


if args.download_files != 'N':
    for i in args.SRR:
        download_data(i)

extract_CDS()
for i in args.SRR:
    kallisto(i)
    SleuthInput()
    Sleuth()




