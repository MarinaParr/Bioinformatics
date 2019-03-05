import os
import zipfile
import csv
import re
import pandas as pd
from sm_annotation_utils import SMInstance


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
            coordinates = re.search('coords="(.*),7"', line).group(1)  # we need only two first coordinates
            title = re.search('title="(.*)"', line).group(1)
            kegg_id = title[0:6]  # kegg_ids are 6-letters
            name = title[7:].split('"')[0].replace('(', '').replace(")", '')
            f1.write(kegg_id + '\t' + coordinates + "\t" +name + "\n")
    f1.close()
    f2.close()
    print('File "kegg_coordinates.csv" is created.')


def check_files_existance(): #check if all necessary files are in wd and creates 'HMDB_KEGG.csv' if needed
    print('Checking files in the directory')
    files = os.listdir(os.getcwd())
    ready_to_run = True
    #This part is for when you manually download data from METASPACE
    #if 'metaspace_annotations.csv' not in files:
    #    print("No input from METASPACE. File 'metaspace_annotations.csv' is expected. The program will stop.")
    #    ready_to_run = False
    if ready_to_run:
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
            d[row[0]] = row[1:]
    return d


#This function for csv files manually downloaded from METASPACE
'''
def annotated_metabolites(): #should be changed to interaction with the python-client for metaspace
    with open('metaspace_annotations.csv') as csvfile:
        table = csv.reader(csvfile, delimiter=',')
        accessions = dict()
        isomer_group = 0
        for row in table:
            isomer_group += 1
            if len(row) == 13 and row[0] != 'group':
                if float(row[7]) <= cutoff:
                    for accession in row[12].split(','):
                        property_float = 1
                        if property_for_radius == 'msm': # [0,1]
                            property_float = float(row[6])
                        if property_for_radius == 'fdr': # usually [0,0.1], so *10 for displaying as radius
                            property_float = float(row[7])*10
                        accessions[accession.replace(' ', '')] = [property_float, isomer_group]
        return accessions, isomer_group #isomer_group now means the number of isomer groups
'''


#This function uses python-client for METASPACE
def annotated_metabolites(df, cutoff, property_for_radius):
    accessions = dict()
    for index, row in df.iterrows():
        if row["fdr"] <= cutoff:
            formula = list(index)[0]
            if formula not in accessions:
                accessions[formula] = [row["fdr"], row["msm"], row["intensity"], row["moleculeIds"]]
            else:
                if accessions[formula][0] > row["fdr"] and accessions[formula][1] < row["msm"]:
                    accessions[formula] = [row["fdr"], row["msm"], row["intensity"], row["moleculeIds"]]
    metabolites = dict()
    isomer_group = 1
    for formula in accessions:
        for hmdb_id in accessions[formula][3]:
            property_float = 1
            if property_for_radius == 'msm':  # [0,1]
                property_float = accessions[formula][1]
            if property_for_radius == 'fdr':  # usually [0,0.1], so *10 for displaying as radius
                property_float = accessions[formula][0] * 10
            if property_for_radius == 'intensity':
                property_float = accessions[formula][2] / 10000 # normalising to 0:1
            metabolites[hmdb_id] = [property_float, isomer_group, accessions[formula][0], accessions[formula][1], formula]
        isomer_group += 1
    return metabolites, isomer_group  # isomer_group here means the number of isomer groups + 1


def main(email, password, names, cutoff, property_for_radius):
    sm = SMInstance()
    sm.login(email, password)
    ds_names = names.split(',')
    if check_files_existance():
        print('All necessary files exist. Mapping starts.')
        hmdb_kegg_d = two_columns_tsv_to_dict('HMDB_KEGG.csv')
        kegg_coordinates_d = two_columns_tsv_to_dict('kegg_coordinates.csv')
        for ds_name in ds_names:
            f = open(ds_name + "_data_for_ili.csv", "w")
            f.write('kegg_name,X,Y,Z,radius,group_of_isomers' + '\n')
            df = pd.DataFrame(sm.dataset(name=ds_name).results(database="HMDB-v4"),columns=['msm', 'moc', 'rhoSpatial', 'rhoSpectral', 'fdr', 'mz', 'moleculeNames','moleculeIds','intensity'])
            annotation = annotated_metabolites(df,cutoff, property_for_radius)
            values = annotation[0]
            isomer_groups = annotation[1] - 1
            for accession in values:
                if accession in hmdb_kegg_d:
                    if hmdb_kegg_d[accession][0] in kegg_coordinates_d:
                        kegg_id = hmdb_kegg_d[accession][0]
                        x_y = kegg_coordinates_d[kegg_id][0]
                        name = kegg_coordinates_d[kegg_id][1]
                        name = name.replace(',', '_')
                        radius = values[accession][0]*10 #reflects the property
                        color = values[accession][1]/isomer_groups
                        #formula = values[accession][4]
                        f.write(name + "," + x_y + ',,' + str(radius) + ',' + str(color) + "\n")
            f.close()
        print('Ready. Files _data_for_ili.csv and KEGG_EC_metabolitemap_bw.png can be drag-and-dropped to ili.embl.de')
