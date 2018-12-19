import os
import zipfile
import csv
import re
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--cutoff', default=1,
                    help='Shows the cut-off for the FDR value. Metabolites with FDR >= cut-off are not shown. Optional. Float in [0:1].')
parser.add_argument('--property', default=str(),
                    help='Shows what property of metabolite to highlight. Optional. Can be "fdr", "msm", "intensity".')

args = parser.parse_args()
property = args.property
cutoff = args.cutoff


def hmdb_kegg():# parsing hmdb, extracting corresponding kegg_ids, writing to 'HMDB_KEGG.csv'
    print('Creating "HMDB_KEGG.csv". It will take some time.')
    f1 = open("HMDB_KEGG.csv", "w")
    with zipfile.ZipFile("hmdb_metabolites.zip") as z:
        with z.open("hmdb_metabolites.xml") as f:
            hmdb_accession, kegg_id = str(), str()
            for line in f:
                line = line.decode('utf-8')
                if line.startswith('  <accession>'):
                    hmdb_accession = re.search('<accession>(.*)</accession>', line).group(1)
                if '<kegg_id>' in line:
                    kegg_id = re.search('<kegg_id>(.*)</kegg_id>', line).group(1)
                if line.startswith('</metabolite>'):
                    if kegg_id != str():
                        f1.write(hmdb_accession + '\t' + kegg_id + '\n')
                    hmdb_accession, kegg_id = str(), str()
    f1.close()
    print('File "HMDB_KEGG.csv" is created.')


def kegg_coordinates(): #parsing 'KEGG PATHWAY: Metabolic pathways - Reference pathway.html' and creating "kegg_coordinates.csv"
    print('Creating "kegg_coordinates.csv"')
    f1 = open('kegg_coordinates.csv', 'w')
    f2 = open('KEGG PATHWAY: Metabolic pathways - Reference pathway.html', 'r')
    for line in f2.readlines():
        if line.startswith('<area shape="circle"'):
            coordinates = re.search('coords="(.*),2"', line).group(1)  # we need only two first coordinates
            title = re.search('title="(.*)"', line).group(1)[0:6]  # kegg_ids are 6-letters
            f1.write(title + '\t' + coordinates + "\n")
    f1.close()
    f2.close()
    print('File "kegg_coordinates.csv" is created.')


def check_files_existance(): #check if all necessary files are in wd and creates 'HMDB_KEGG.csv' if needed
    print('Checking files in the directory')
    files = os.listdir(os.getcwd())
    ready_to_run = True
    if 'metaspace_annotations.csv' not in files:
        print("No input from METASPACE. File 'metaspace_annotations.csv' is expected. The program will stop.")
        ready_to_run = False
    if ready_to_run == True:
        if 'HMDB_KEGG.csv' not in files:
            if 'hmdb_metabolites.zip' in files:
                hmdb_kegg()
            else:
                ready_to_run = False
                print("File 'hmdb_metabolites.zip' is needed. http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip")
        if 'kegg_coordinates.csv' not in files:
            if 'KEGG PATHWAY: Metabolic pathways - Reference pathway.html' in files:
                kegg_coordinates()
            else:
                ready_to_run = False
                print("File 'KEGG PATHWAY: Metabolic pathways - Reference pathway.html' is needed.")
    return ready_to_run


def two_columns_tsv_to_dict(tsv_file): #reads tsv_file with two columns into dict: dict[row[0]] = row[1]
    d = dict()
    with open(tsv_file) as csvfile:
        table = csv.reader(csvfile, delimiter='\t')
        for row in table:
            d[row[0]] = row[1]
    return d


def annotated_metabolites(): #should be changed to interaction with the python-client for metaspace
    with open('metaspace_annotations.csv') as csvfile:
        table = csv.reader(csvfile, delimiter=',')
        accessions = dict()
        for row in table:
            if len(row) == 13 and row[0] != 'group':
                if float(row[7]) <= cutoff:
                    for accession in row[12].split(','):
                        color = 1
                        if property == 'msm':
                            color = row[6]
                        if property == 'fdr':
                            color = row[7]
                        accessions[accession.replace(' ', '')] = color
        return accessions


def main():
    if check_files_existance():
        print('All necessary files exist. Mapping starts.')
        hmdb_kegg = two_columns_tsv_to_dict('HMDB_KEGG.csv')
        kegg_coordinates = two_columns_tsv_to_dict('kegg_coordinates.csv') #here second row contains x,y coordinates
        f = open("data_for_ili.csv", "w")
        f.write('kegg_id,X,Y,Z,radius,dummy' + '\n')
        annotation = annotated_metabolites()
        for accession in annotation:
            if accession in hmdb_kegg:
                if hmdb_kegg[accession] in kegg_coordinates:
                    f.write(hmdb_kegg[accession] + "," + kegg_coordinates[hmdb_kegg[accession]] + ',,' + '10,' + str(annotation[accession]) + "\n")
        f.close()
        print('Ready. Files "data_for_ili.csv" and KEGG_EC_metabolitemap_bw.png (1320x790pi) can be drag-and-dropped to https://ili.embl.de')


main()
