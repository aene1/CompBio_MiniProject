# CompBio_MiniProject
Author : Adriana Ene

***Software Requirments:***
* Linux/Unix
* Python3
* Biopython
* Kallisto
* Bowtie2
* SPAdes

Python Script

There are 2 options of running the script. If the SRR files are already downloaded in the directory those are the files
that will be used to run the script. If they are not, they will automatically get downloaded from the NCBI SRR database. 

#### INSTALL

    git clone https://github.com/aene1/CompBio_MiniProject.git
    
#### RUN
    python3 mini_project.py <SRR1> <SRR2> <SRR3> <SRR4>
    
#### EXAMPLE
    python3 mini_project.py SRR5660030 SRR5660033 SRR5660044 SRR5660045
    
## TO RUN TEST DATA ##
Go to the TEST directory and run the following command:

    python3 mini_project.py SRR5660030 SRR5660033 SRR5660044 SRR5660045 --Test

The test data contains the first 10000 lines of the paired fastq files. Which makes it significantly faster to run.
***Files included in Repo***

* mini_project.py
   >  the whole pipeline which is composed of multiple functions that give each outputs
                    
* sleuth.R
    > R script which outputs a text while which identifies the diffrence between 2 timepoints in the expressed genes 

* sequences.fa
    > nucleotide database out of all RefSeq sequences of the Betaherpesvirinae found on 3/4/2021
                 
#### Output files (the most important ones)
* miniProject.log
    >a log file that has information from the data ran with the script
*  HCMV_index.idx
    > index created to be used by kallisto
* CDS_EF999921.fa 
    > fasta file with CDS of EF999921
