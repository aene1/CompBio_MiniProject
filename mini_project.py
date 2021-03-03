import argparse
import os
from Bio import SeqIO
from Bio import Entrez
import glob

log_file = open('miniProject.log', 'w+')
path = os.getcwd()


# speed up spades run

# download the SRR and fastq them
def download_data(SRR):
    wget = 'wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/' + SRR + '/' + SRR + '.1'
    rename = 'mv ' + SRR + '.1 ' + SRR
    os.system(wget)
    os.system(rename)


def fastq(SRR):
    fastq = 'fastq-dump -I --split-files ' + SRR + '.1'
    os.system(fastq)
    print(fastq)


# extract the CDS in order to construct the transcriptome
def extract_CDS():
    Entrez.email = 'aene@luc.edu'
    handle = Entrez.efetch(db='nucleotide', id='EF999921', rettype='gb', retmode='text')
    records = SeqIO.parse(handle, 'genbank')
    count = 0
    outfile = open('CDS_EF999921.fa', 'w')  # the outfile where the CDS records will  be written out to
    for i in records:
        for j in i.features:
            if j.type == 'CDS':
                outfile.write('>' + str(j.qualifiers['protein_id']).replace("'", '').replace('[', '').replace(']',
                                                                                                              '') + '\n' + str(
                    j.location.extract(i).seq) + '\n')
                count += 1
    log_file.write('The HCMV genome (EF999921) has ' + str(count) + '\n')


# runn kallisto based on SRR numbers and the SRR files dowloaded before
def kallisto(SRR):
    kallisto_index = 'time kallisto index -i HCMV_index.idx CDS_EF999921.fa'
    os.system(kallisto_index)
    run_kallisto = 'time kallisto quant -i HCMV_index.idx -o' + path + '/results_' + SRR + ' -b 30 -t 4 ' + SRR + '.1_1.fastq ' + SRR + '.1_2.fastq'
    os.system(run_kallisto)


def SleuthInput(SRR):
    # output file that goes in R
    output = open('input_sleuth.txt', 'w')
    # initial line in file
    output.write('sample' + '\t' + 'condition' + '\t' + 'path' + '\n')
    # based on SRR number, write condition and path to output file
    path = os.getcwd()
    for i in SRR:
        path1 = path + '/' + 'results_' + i
        print(path1)
        if int(i[
               3:]) % 2 == 0:  # if it is even then it is condition 1 as in 2dpi, if it is not then it is condition 2 as in 6dpi
            output.write(str(i) + '\t' + '2dpi' + '\t' + path1 + '\n')
        else:
            output.write(str(i) + '\t' + '6dpi' + '\t' + path1 + '\n')
    output.close()


def Sleuth():
    Sleuth_command = 'Rscript sleuth.R'
    os.system(Sleuth_command)
    read_sleuth = open('R_sleuth_output.txt', 'r').readlines()
    for i in read_sleuth:
        log_file.write(i + 'n')


# build bowtie2 and run it to get SAM files

def index_bowtie():
    handle = Entrez.efetch(db='nucleotide', id='EF999921', rettype='fasta', retmode='text')
    records = list(SeqIO.parse(handle, 'fasta'))
    for i in records:
        print(i)
        SeqIO.write(i, 'EF999921.fa', "fasta")


def bowtie2(SRR):
    build_bowtie2 = 'bowtie2-build EF999921.fa EF999921_1'
    os.system(build_bowtie2)
    bowtie2 = 'bowtie2 --quiet --no-unal --al-conc bowtie2_' + SRR + '.fastq -x EF999921_1 -1 ' + SRR + '.1_1.fastq -2 ' + SRR + '.1_2.fastq -S EF999921_'
    os.system(bowtie2)


def Count_bowtie(SRR):
    bowtie_SRR1 = open('bowtie2_' + SRR + '.1.fastq').readlines()
    bowtie_SRR2 = open('bowtie2_' + SRR + '.2.fastq').readlines()
    donor = ''
    if SRR == "SRR5660030":
        donor += 'Donor 1 (2dpi)'
    elif SRR == 'SRR5660033':
        donor += 'Donor 1 (6dpi)'
    elif SRR == 'SRR5660044':
        donor += 'Donor 3 (2dpi)'
    elif SRR == 'SRR5660045	':
        donor += 'Donor 3 (6dpi)'
    len_bowtie = ((len(bowtie_SRR1) + len(bowtie_SRR2)) / 8)
    original1 = open(SRR + '.1_1.fastq').readlines()
    original2 = open(SRR + '.1_2.fastq').readlines()
    original = (len(original1) + len(original2)) / 8
    # write out to the log file
    log_file.write(donor + " had " + str(original) + ' read pairs before Bowtie2 filtering and ' + str(len_bowtie) + ' read pairs after \n')


def run_spades(SRR1, SRR2, SRR3, SRR4):  # use only assembler to make Spades a lot shorter
    path = os.getcwd()
    SRR1 = path + "/" + 'bowtie2_' + SRR1
    SRR2 = path + "/" + 'bowtie2_' + SRR2
    SRR3 = path + "/" + 'bowtie2_' + SRR3
    SRR4 = path + "/" + 'bowtie2_' + SRR4
    spades_command = 'spades -k 55,77,99,127 --only-assembler -t 2 --pe1-1 ' + SRR1 + '.1.fastq --pe1-2 ' + SRR1 + '.2.fastq --pe2-1 '+ SRR2 + '.1.fastq --pe2-2  ' + SRR2 + '.2.fastq --pe3-1 '+ SRR3 + '.1.fastq --pe3-2 '+ SRR3 + '.2.fastq --pe4-1 ' + SRR4 + '.1.fastq --pe4-2 ' + SRR4 + '.2.fastq -o '+ path + '/Spades/'
    # run SPades and print command to log file
    print(spades_command)
    log_file.write(spades_command + '\n')


# def contig_calc():
#     os.chdir(path + '/spades_assembly')  # navigate to spades_assembly folder
#     record = SeqIO.parse('contigs.fasta', 'fasta')
#     count = 0
#     seq_len = 0
#     longest
#     for i in record:
#         if len(i.seq) > 1000:
#             count += 1
#     log_file.write('There are ' + str(count) + ' contigs > 1000 bp in the assembly.')
#     log_file.write('\n')


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
print("INDEX BOWTIE BUILT")

# for i in args.SRR:
#     fastq(i)

for i in args.SRR:
    kallisto(i)
    print("BUILD INDEX")
    bowtie2(i)
    print("RAN BOWTIE")
#
SleuthInput(args.SRR)
print('SLEUTHinput WORKED')
Sleuth()
#
#
# for i in args.SRR:
#     Count_bowtie(i)
#
# run_spades(args.SRR[0], args.SRR[1], args.SRR[2], args.SRR[3])
