"""Ce code permet de créer des familles de séquences à partir de fichiers
fasta.

Usage:
------
    python3 create_family.py arguments

    arguments: le ou les fichier.s fasta à analyser
"""


############ Modules à importer ############

import sys
from tqdm import tqdm
from Levenshtein import *
import operator
import copy

############################################


def arguments():
    """Vérifier le format et le nombre d'arguments renseigné.

    Returns
    -------
    fichiers: liste de tous les fichiers donnés en argument.
    """

    fichiers = []
    
    if len(sys.argv) < 2:
        sys.exit("Veuillez renseigner au moins un fichier fasta à lire")
    
    for index in range(1, len(sys.argv)):
        if not sys.argv[index].endswith(".fas"):
            sys.exit("Les fichiers renseignés doivent être au format fasta")
        fichiers.append(str(sys.argv[index]))
    
    return fichiers


def save_data(fichier):
    """Lit un fichier fasta.

    Parameters
    ----------
    fichier : string
        fichier fasta à lire

    Returns
    -------
    extracted_data: list
        liste de tuple contenant le nom des séquences et les séquences.
    """
    
    name_comments = []
    sequence = []
    extracted_data = []

    with open(fichier, "r") as filin:
        lines = filin.readlines()
        for line in lines:
            if line.startswith(">"):
                name_comments.append(line.strip())
            else:
                sequence.append(line.strip())
        extracted_data = [tuple(name_comments), tuple(sequence)]
    
    return extracted_data


def kept_all(extracted_all_data):
    """lit une liste.

    Parameters
    ----------
    fichier : list
        liste de toutes les séquences de tous les fichiers fastas

    Returns
    -------
    wanted_seq: dictionnary
        dictionnaire contenant toutes les séquences présentes à 
        1000 occurrences ou plus.
    
    compte: dictionnary
        dictionnaire de comptage contenant le nombre de fois où chaque
        séquence apparait dans tous les fichiers fastas.
    """    

    compte = {}.fromkeys(set(extracted_all_data), 0)
    
    for seq in tqdm(extracted_all_data):
        compte[seq] += 1

    wanted_seq = {}
    
    for cle, valeur in tqdm(compte.items()):
        if valeur >= 1000:
            wanted_seq[cle] = valeur
    
    return wanted_seq, compte


def kept_data(extracted_data):
    """lit une liste.

    Parameters
    ----------
    fichier : list
        liste de toutes les séquences d'un fichier fasta

    Returns
    -------
    wanted_seq: dictionnary
        dictionnaire contenant toutes les séquences présentes à 
        1000 occurrences ou plus.
    """

    compte = {}.fromkeys(set(extracted_data[1]), 0)
    
    for seq in tqdm(extracted_data[1]):
        compte[seq] += 1

    wanted_seq = {}
    
    for cle, valeur in tqdm(compte.items()):
        if valeur >= 1000:
            wanted_seq[cle] = valeur
    
    return wanted_seq


def save_data_in_txt_file(data, fichier):
    """sauvegarde un dictionnaire dans un fichier texte.

    Parameters
    ----------
    data : dictionnary
        dictionnaire contenant les valeurs à sauvegarder

    fichier: string
        le nom du fichier dans lequel les sauvegarder
    """

    with open(fichier, "w") as filout:
        for cle, valeur in data.items():
            filout.write(f"{cle} {valeur}\n")


def save_one_key_dict_in_txt_file(data, fichier):
    """sauvegarde les données d'une clé d'un dictionnaire dans 
    un fichier texte.

    Parameters
    ----------
    data : dictionnary values
        valeurs contenues dans une clé d'un dictionnaire

    fichier: string
        le nom du fichier dans lequel les sauvegarder
    """
    with open(fichier, "w") as filout:
        for seq in data:
            filout.write(f"{seq}\n")


def extract_max_in_dict(dictionnaire):
    """donne le maximum d'un dictionnaire.

    Parameters
    ----------
    dictionnaire : dictionnary
        un dictionnaire contenant des valeurs numériques

    returns
    -------

    number_seq_max : int
        la valeur maximale du dictionnaire
    
    seq_max : string
        la séquence associée à cette valeur maximale
    """
    number_seq_max = max(dictionnaire.items(), key = operator.itemgetter(1))[1]
    seq_max = max(dictionnaire.items(), key = operator.itemgetter(1))[0]
    
    return number_seq_max, seq_max


def create_families(compte_all, nombre_famille, dist_max):
    """créer des familles de séquences.

    Parameters
    ----------
    compte_all : dictionnary
        un dictionnaire contenant le nombre d'occurences de
        chaque séquences au sein de tous les fichiers fastas

    nombre_famille: int
        le nombre de famille voulu

    dist_max: int
        la distance de Levenshtein maximum entre deux séquences
        pour qu'elles soient considéré comme faisant parti de la
        même famille

    returns
    -------

    fam_seq : dictionnary
        dictionnaire contenant les familles de séquences,
        avec en clé le numéro de famille et en valeur associée
        la liste des séquences différentes présentes dans cette famille.

    seq_fam: dictionnary
        dictionnaire contenant les familles de séquences,
        avec en clé la séquence et en valeur associée
        le numéro de famille.
    
    fam_seq_complete: dictionnary
        dictionnaire contenant les familles de séquences,
        avec en clé le numéro de famille et en valeur associée
        la liste des séquences totales présentes dans cette famille.

    seq_ref: dictionnary
        dictionnaire contenant la séquence de référence pour chaque famille
        ainsi que le nombre de fois où elle apparait dans la famille.
    """    
    fam_seq = {}
    seq_fam = {}
    fam_seq_complete = {}
    seq_ref = {}
    num_famille = 0
    dict_miroir = copy.deepcopy(dictionnaire)
    
    while dict_miroir:

        if len(fam_seq) == nombre_famille:
            break
        
        number_seq_max, seq_max = extract_max_in_dict(dictionnaire)
        num_famille += 1
        fam_seq[num_famille] = []
        fam_seq_complete[num_famille] = []
        fam_seq[num_famille].append(seq_max)
        tmp = seq_max.split()
        fam_seq_complete[num_famille] += tmp * number_seq_max
        seq_fam[seq_max] = num_famille
        seq_ref[num_famille] = [(seq_max), (number_seq_max)]
        dictionnaire.pop(seq_max)
        dict_miroir.pop(seq_max)
        print(num_famille)

        for cle, valeur in tqdm(dictionnaire.items()):
            if (cle, valeur) in dict_miroir.items():
                diff = distance(cle, seq_max)
                if diff <= dist_max:
                    fam_seq[num_famille].append(cle)
                    seq_fam[cle] = num_famille
                    tmp = []
                    tmp.append(cle)
                    fam_seq_complete[num_famille] += tmp * valeur
                    dict_miroir.pop(cle)
        
        for seq in tqdm(fam_seq[num_famille]):
            if seq in dictionnaire:
                dictionnaire.pop(seq)
    
    return fam_seq, seq_fam, fam_seq_complete, seq_ref

############################################

#def common_string_two_list(a, b):
#
#    return set(a) & set(b)

############################################
############################################

#def diff_string_two_list(a, b):
#
#    return set(a) ^ set(b)

############################################
############################################

#def count_diff_and_common_seq_in_files(dictionnaire):
    
#    nbr_seq_eq = {}
#    nbr_seq_diff = {}
#    
#    for cle, valeur in tqdm(dictionnaire.items()):
#        for cle2, valeur2 in dictionnaire.items():
#            if cle == cle2:
#                continue
#            else:
#                nbr_seq_eq[int(cle.strip('R')), int(cle2.strip('R'))] = len(common_string_two_list(dictionnaire[cle], dictionnaire[cle2]))
#                nbr_seq_diff[int(cle.strip('R')), int(cle2.strip('R'))] = len(diff_string_two_list(dictionnaire[cle], dictionnaire[cle2]))
#
#
#    
#    return nbr_seq_eq, nbr_seq_diff

############################################
############################################

#def normalize_values(dictionnaire):
#
#    nbr_seq_diff_in_rounds = [1586796, 894642, 1132902, 978828, 928735,
#    657549, 677082, 304929, 193240, 225455, 154717, 239797, 359010, 105078,
#    141490, 123612, 131028, 138229]
#    normalized_values = {}
#    keys = []
#
#    for cle in dictionnaire.keys():
#        if cle[0] not in keys:
#            keys.append(cle[0])
#    
#    diff_seq_in_Round = dict(zip(keys, nbr_seq_diff_in_rounds))
#
#    for cle, valeur in dictionnaire.items():
#        for k, v in diff_seq_in_Round.items():
#            for i in keys:
#                if cle[0] == i:
#                    normalized_values[cle] = 0
#                    normalized_values[cle] = valeur / min(diff_seq_in_Round[cle[0]], diff_seq_in_Round[cle[1]])
#
#    return normalized_values

############################################

def main():
    """Le main du programme."""

    dist_max = int(input("Quelle est la distance de Levenshtein maximum entre deux séquences pour former une famille ? "))
    nombre_famille = int(input("Quel est le nombre de famille souhaité ? "))
    
    extracted_all_data_dict = {}
    extracted_all_data_list = []
    fichiers = arguments()
    
    for fichier in fichiers:
        extracted_data = save_data(fichier)
        wanted_seq = kept_data(extracted_data)
        extracted_all_data_dict[fichier[0:3]] = list(extracted_data[1])
        extracted_all_data_list = extracted_all_data_list + list(extracted_data[1])
        seq_kept = fichier.strip(".fastq_result.fas") + "_kept_data.txt"
        save_data_in_txt_file(wanted_seq, seq_kept)
    
    wanted_seq_all, compte_all = kept_all(extracted_all_data_list)
    fam_seq, seq_fam, fam_seq_complete, seq_ref = create_families(compte_all, nombre_famille, dist_max)

    save_seq_ref = "seq_de_reference_pour_familles.txt"
    save_data_in_txt_file(seq_ref, save_seq_ref)

    for cle, valeur in fam_seq.items():
        save_by_family_diff_seq = f"family_{cle}_diff_seq.txt"
        save_one_key_dict_in_txt_file(fam_seq[cle], save_by_family_diff_seq)
    
    for cle, valeur in fam_seq_complete.items():
        save_by_family_all_seq = f"family_{cle}_all_seq.txt"
        save_one_key_dict_in_txt_file(fam_seq_complete[cle], save_by_family_all_seq)


    seq_kept_all = "seq_sup_1000_occ.txt"
    save_data_in_txt_file(wanted_seq_all, seq_kept_all)


#    nbr_seq_eq, nbr_seq_diff = count_diff_and_common_seq_in_files(extracted_all_data_dict)
#    normalized_values = normalize_values(nbr_seq_eq)
#    normalized_values_diff = normalize_values(nbr_seq_diff)
#    seq_diff_in_files_norm = "seq_diff_in_files_norm.txt"
#    common_seq_in_files_norm = "common_seq_in_files_norm.txt"
#    save_data_in_txt_file(normalized_values_diff, seq_diff_in_files_norm)
#    save_data_in_txt_file(normalized_values, common_seq_in_files_norm)
#    common_seq_in_files = "common_seq_in_files.txt"
#    seq_diff_in_files = "seq_diff_in_files.txt"
#    save_data_in_txt_file(nbr_seq_eq, common_seq_in_files)
#    save_data_in_txt_file(nbr_seq_diff, seq_diff_in_files)



if __name__ == "__main__":
    main()