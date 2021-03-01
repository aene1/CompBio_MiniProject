import argparse
import os
from Bio import SeqIO
from Bio import Entrez

log_file = open('miniProject.log','w')



#download the SRR and fastq them
# def download_data(SRR):
#     wget = 'wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/' +  SRR + '/' + SRR + '.1'
#     fastq = 'fastq-dump -I --split-files ' +SRR + '.1'
#     rename = 'mv ' + SRR + '.1' +' ' + SRR
#     os.system(wget)
#     os.system(rename)
#     os.system(fastq)


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







# parser = argparse.ArgumentParser(description='Process some SRRs and split paired reads.')
# parser.add_argument('SRR', metavar='N', type=str, nargs='+',
#                     help='Compare SRR files')
# # parser.add_argument('--download_SRR')
# args = parser.parse_args()

# for i in args.SRR:
#     download_data(i)

extract_CDS()