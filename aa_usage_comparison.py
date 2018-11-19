#This program compares the occurence of amino acids in two sets of proteomes
from Bio import SeqIO
import os
import scipy.stats as ss
import matplotlib.pyplot as plt

#input names of groups and paths to folders
name_set_1 = str(input("Print name of 1st group: ")) #example: name_1
path_to_set_1 = str(input("Print path to 1st group: ")) #example: folder_1/
name_set_2 = str(input("Print name of 2nd group: ")) #example: name_2
path_to_set_2 = str(input("Print path to 2nd group: ")) #example: folder_2/
path_to_results = str(input("Print path to folder with results: ")) #example: results_folder/

if not os.path.exists(path_to_results): #creating folder for result figures
    os.mkdir(path_to_results)

amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

def creating_dict(type_of_elements): #creates dict d: d[amino_acid] == 0 or d[amino_acid] == list()
    d = dict()
    if type_of_elements == 'ints':
        for amino_acid in amino_acids:
            d[amino_acid] = 0
    if type_of_elements == 'lists':
        for amino_acid in amino_acids:
            d[amino_acid] = list()
    return d


def aa_usage_in_proteome(proteome_name, path_to_set): #calculates usage of amino acid in one particular proteome
    proteins_counter = 0 #calculates number of proteins
    d = creating_dict('ints')
    for record in SeqIO.parse(path_to_set + proteome_name, "fasta"): #parsing .fasta file
        sequence = str(record.seq)
        if len(sequence) != 0:
            for amino_acid in amino_acids:
                d[amino_acid] += sequence.count(amino_acid)/len(sequence) #counts aa usage in one sequence and sum it to total score
            proteins_counter += 1
    for amino_acid in d:
        d[amino_acid] = round(d[amino_acid]/proteins_counter, 3) #divide the total score on the number of proteins
    return d


def creating_distribution_in_set(path_to_set): #for each aa creates distribution of frequencies in groups of proteomes
    d = creating_dict('lists')
    for proteome in os.listdir(path_to_set):
        results_for_proteome = aa_usage_in_proteome(proteome, path_to_set)
        for amino_acid in amino_acids:
            d[amino_acid].append(results_for_proteome[amino_acid])
    return d


def plotting(l1, l2, text, amino_acid, name_1, name_2): #plots boxplots of distributions
    fig = plt.figure(1, figsize=(9, 6))
    ax = fig.add_subplot(111)
    ax.boxplot([l1, l2])
    ax.set_xticklabels([name_1, name_2])
    ax.set_title(amino_acid + ":" + text)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.xlabel("Groups of proteomes")
    plt.ylabel("Occurrence of amino acid " + amino_acid)
    fig.savefig(path_to_results + "/" + amino_acid + '.png', bbox_inches='tight')
    plt.clf()


def main(): #calculates, compares, plots
    dists_1 = creating_distribution_in_set(path_to_set_1)
    dists_2 = creating_distribution_in_set(path_to_set_2)
    for amino_acid in amino_acids:
        dist_set1 = dists_1[amino_acid]
        dist_set2 = dists_2[amino_acid]
        p_value = ss.mannwhitneyu(dist_set1, dist_set2)[1]#Mann-Whitney test is used for comparison
        if p_value < 0.05:
            text = ' signigicant difference, p-value = '+ str(round(p_value, 7))
        else:
            text = ' no significant difference, p-value = ' + str(round(p_value, 7))
        plotting(dist_set1, dist_set2, text, amino_acid, name_set_1, name_set_2)

main()
