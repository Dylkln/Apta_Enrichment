"""Ce code permet de l'entropie de shannon de chaque famille créée.

Usage:
------
    python3 entropy.py arguments

    arguments: le ou les fichier.s texte à analyser
"""


############ Modules à importer ############


import sys
from tqdm import tqdm
import math


############################################

def arguments():
    """Vérifier le format et le nombre d'arguments renseigné

    Returns
    -------
    fichiers: liste de tous les fichiers donnés en argument
    """

    fichiers = []
    
    if len(sys.argv) < 2:
        sys.exit("Veuillez renseigner au moins un fichier txt à lire")
    
    for index in range(1, len(sys.argv)):
        if not sys.argv[index].endswith(".txt"):
            sys.exit("Les fichiers renseignés doivent être au format txt")
        fichiers.append(str(sys.argv[index]))
    
    return fichiers


def save_data(fichier):
    """Lit un fichier.

    Parameters
    ----------
    fichier : string
        fichier texte à lire

    Returns
    -------
    sequence: list
        liste de des séquences présentes dans cette famille
    taille_famille : int
        la taille de la famille
    """ 
    sequence = []
    taille_famille = 0
    
    with open(fichier, "r") as filin:
        lines = filin.readlines()
        taille_famille = len(lines)
        for line in lines:
            sequence.append(line.strip())
    
    return sequence, taille_famille


def compte_seq(sequence):
    """compte le nombre d'occurrences de chaque séquence
    dans la liste de séquences.

    Parameters
    ----------
    sequence: list
        liste des séquences de la famille

    Returns
    -------
    compte: dictionnary
        dictionnaire contenant le nombre d'occurrences de chaque séquence
        de la liste renseignée en argument
    """ 

    compte = {}.fromkeys(set(sequence), 0)

    for seq in tqdm(sequence):
        compte[seq] += 1

    return compte


def calc_freq(compte, taille_famille):
    """calcule la fréquence de chaque séquence de la famille.

    Parameters
    ----------
    compte: dictionnary
        dictionnaire contenant le nombre d'occurrences de chaque séquence
        de la liste renseignée en argument

    taille_famille : int
        la taille de la famille        

    Returns
    -------
    frequence: dictionnary
        dictionnaire contenant la fréquence de chaque séquence de la famille
    """ 

    frequence = {}.fromkeys(set(compte.keys()), 0)
   
    for cle, valeur in compte.items():
        frequence[cle] += valeur / taille_famille

    return frequence


def calc_shannon_entropy(frequence, fichier, taille_seq_uniques):
    """calcule l'entropie de Shannon de la famille.

    Parameters
    ----------
    frequence: dictionnary
        dictionnaire contenant la fréquence de chaque séquence de la famille

    fichier: string
        le nom du fichier renseigné en argument
    
    taille_seq_uniques : dictionnary
        dictionnaire contenant la taille de chaque famille créée        

    Returns
    -------
    entropy: float
        entropy de Shannon pour la famille testée
    """ 

    cle = int(fichier.strip("family__all_seq.txt"))
    base = taille_seq_uniques[cle]
    entropy = 0
    
    for valeur in frequence.values():
        entropy += - valeur * math.log(valeur, base)
    
    return entropy


def save_dict(entropy, fichier):
    """sauvegarde un dictionnaire dans un fichier texte.

    Parameters
    ----------
    entropy : dictionnary
        dictionnaire contenant les valeurs à sauvegarder

    fichier: string
        le nom du fichier dans lequel les sauvegarder
    """
    with open(fichier, "w") as filout:
        for cle, valeur in entropy.items():
            filout.write(f"{cle} {valeur}\n")


def create_dict_taille_seq_uniques():
    """créer un dictionnaire contenant la taille de chaque famille.
    
    returns
    -------

    taille_seq_uniques: dictionnary
        dictionnaire contenant la taille de chaque famille créée  
    """    
    taille_seq_uniques = {}
    
    with open("taille_des_familles.txt", "r") as filin:
        lines = filin.readlines()
        for line in lines:
            taille_seq_uniques[int(line[0:3].strip(" \n "))] = 0
            taille_seq_uniques[int(line[0:3].strip(" \n "))] += int(line[4:13].strip(" \n "))
    
    return taille_seq_uniques

def main():
    """Le main du programme."""   

    taille_seq_uniques = create_dict_taille_seq_uniques()
    entropy_dict = {}
    fichiers = arguments()
    

    for fichier in fichiers:
        sequence, taille_famille = save_data(fichier)
        compte = compte_seq(sequence)
        frequence = calc_freq(compte, taille_famille)
        entropy = calc_shannon_entropy(frequence, fichier, taille_seq_uniques)
        entropy_dict[fichier.strip("family__all_seq.txt")] = 0
        entropy_dict[fichier.strip("family__all_seq.txt")] += entropy
        freq = "frequence_par_seq_famille_{}_all_seq.txt".format(fichier.strip("family__all_seq.txt"))
        frequence = dict(sorted(frequence.items(), key = lambda t: t[1], reverse = True))
        save_dict(frequence, freq)
        compte_par_famille = "nbr_occ_seq_famille_{}.txt".format(fichier.strip("family__all_seq.txt"))
        compte = dict(sorted(compte.items(), key = lambda t: t[1], reverse = True))
        save_dict(compte, compte_par_famille)
    
    entropy_file = "shannon_entropy.txt"
    save_dict(entropy_dict, entropy_file)


if __name__ == "__main__":
    main()