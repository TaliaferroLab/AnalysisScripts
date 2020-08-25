#!/usr/bin/env python3

import argparse
import logging
import subprocess
import os

def downloadfastqs(srr, ssh_key, output_directory):
        # Get the textual representation of the run. We specifically need the fastq_ftp bit
    logging.info("Querying ENA for FTP paths ..")
    text = subprocess.check_output("curl --silent 'https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={}&result=read_run&fields=fastq_ftp&download=txt'".format(
        srr),shell=True)

    ftp_urls = []
    header=True
    for line in text.decode('utf8').split('\n'):
        if header:
            header=False
        else:
            for url in line.split(';'):
                if url.strip() != '': ftp_urls.append(url.strip())
    logging.info("Found {} FTP URLs for download e.g. {}".format(len(ftp_urls), ftp_urls[0]))

    '''
    aspera_commands = []
    for url in ftp_urls:
        cmd = "ascp -QT -l 300m -P33001 -i {} era-fasp@fasp.sra.ebi.ac.uk:{} {}".format(
            ssh_key_file,
            url.replace('ftp.sra.ebi.ac.uk',''), output_directory)
        logging.info("Running command: {}".format(cmd))
        subprocess.check_call(cmd,shell=True)
    '''

    for url in ftp_urls:
        cmd = 'wget {0}'.format(url)
        subprocess.check_call(cmd, shell = True)

def renamefiles(srr, samplename, output_directory):
    allfastqs = os.listdir(output_directory)
    srrfastqs = [f for f in allfastqs if srr in f and f.endswith('.fastq.gz')]
    
    #singleend
    if len(srrfastqs) == 1:
        oldfilename = os.path.join(output_directory, srr + '.fastq.gz')
        newfilename = os.path.join(output_directory, samplename + '.fastq.gz')
        os.rename(oldfilename, newfilename)

    #pairedend
    elif len(srrfastqs) >= 2:
        oldfilename1 = os.path.join(output_directory, srr + '_1.fastq.gz')
        newfilename1 = os.path.join(output_directory, samplename + '_1.fastq.gz')
        os.rename(oldfilename1, newfilename1)

        oldfilename2 = os.path.join(output_directory, srr + '_2.fastq.gz')
        newfilename2 = os.path.join(output_directory, samplename + '_2.fastq.gz')  
        os.rename(oldfilename2, newfilename2)      




if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Quickly download FASTQ files from the European Nucleotide Archive (ENA) using aspera.\n\n'
        'Requires curl and ascp (i.e. aspera, see https://www.biostars.org/p/325010/#389254) to be in the $PATH.')
    parser.add_argument('--downloadlist',help='Tab separated list of files to downloaded. Column 1 is SRR id and Column 2 is sample name.')
    parser.add_argument('--output_directory',help='Output files to this directory [default: \'.\']',default='.')
    parser.add_argument('--ssh_key',help='\'linux\' or \'osx\' for default paths used in each OS respectively, \
    otherwise a path to the openssh key to used for aspera (i.e. the -i flag of ascp) [default: \'linux\']',
                        default='linux')
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # Use the example set out at this very helpful post:
    # https://www.biostars.org/p/325010

    if args.ssh_key == 'linux':
        ssh_key_file = '$HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh'
    elif args.ssh_key == 'osx':
        ssh_key_file = '$HOME/Applications/Aspera\ Connect.app/Contents/Resources/asperaweb_id_dsa.openssh'
    else:
        ssh_key_file = args.ssh_key
    logging.info("Using aspera ssh key file: {}".format(ssh_key_file))

    output_directory = args.output_directory

    srrs = []
    samplenames = []

    with open(args.downloadlist, 'r') as infh:
        for line in infh:
            line = line.strip().split('\t')
            srrs.append(line[0])
            samplenames.append(line[1])

    for idx, srr in enumerate(srrs):
        samplename = samplenames[idx]
        downloadfastqs(srr, ssh_key_file, output_directory)
        renamefiles(srr, samplename, output_directory)


logging.info("All done.")
