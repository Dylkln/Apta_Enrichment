"""Ce code permet de calculer l'entropie de Shannon à chaque position
des séquences d'une famille à chaque round.

Usage:
------
    python3 mutation.py arguments

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
    data: dictionnary
        dictionnaire contenant le profil de la famille pour le round donné
    """ 

    data = {}
    with open(fichier, "r") as filin:
        lines = filin.readlines()
        lines = filter(str.strip, lines)
        for line in lines:
            if line[0:2] not in data.keys():
                data[line[0:2]] = {}
                data[line[0:2]][line[3]] = int(line[4:20].strip())
            if line[0:2] in data.keys():
                data[line[0:2]][line[3]] = int(line[4:20].strip())

    return data

def mutation(data, fichier, longueur_seq, seq_by_round):
    """Calcule la différence moyenne entre les séquences.

    Parameters
    ----------
    data: dictionnary
        dictionnaire contenant le profil de la famille pour le round donné

    fichier : string
        fichier texte à lire

    longueur_seq: int
        la longueur des séquences du fichier

    seq_by_round: dictionnary
        dictionnaire contenant le nombre de séquences présentes par round
        pour la famille testée

    Returns
    -------
    diff_moy_seq: float
        différence moyenne entre deux séquences pour la famille testée
    """
    somme_base = 0
    diff_moy_seq = 0

    taille_fam = seq_by_round[fichier[16:19]]

    for cle, valeur in data.items():
        for cle2, valeur2 in valeur.items():
            tmp = valeur2 - 1
            somme_base += valeur2 * tmp

    tmp = taille_fam - 1
    taille_fam = taille_fam * tmp
    tmp = somme_base / taille_fam
    diff_moy_seq = longueur_seq - tmp


    return diff_moy_seq

    

def frequence_par_round(data, seq_by_round):
    """calcule la fréquence de base azotée par round.

    Parameters
    ----------
    data: dictionnary
        dictionnaire contenant le profil de la famille pour le round donné

    seq_by_round: dictionnary
        dictionnaire contenant le nombre de séquences présentes par round
        pour la famille testée       

    Returns
    -------
    freq_par_round: dictionnary
        dictionnaire contenant la fréquence de chaque base azotée par round
    """     
    freq_par_round = {}

    bases = ["A", "T", "C", "G"]

    for cle in seq_by_round.keys():
        freq_par_round[cle] = {}

    for cle, valeur in data.items():
        for cle2, valeur2 in valeur.items():
            freq_par_round[cle][cle2] = {}
            for cle3, valeur3 in valeur2.items():
                taille = seq_by_round[cle]
                freq_par_round[cle][cle2][cle3] = valeur3 / taille
                for base in bases:
                    if base not in valeur2.keys():
                        freq_par_round[cle][cle2][base] = 0




    return freq_par_round


def calc_shannon_entropy(freq_par_round):
    """calcule l'entropie de Shannon pour chaque position.

    Parameters
    ----------
    freq_par_round: dictionnary
        dictionnaire contenant la fréquence de chaque base azotée par round

    Returns
    -------
    entropy_dict: dictionnary
        dictionnaire contenant l'entropie de Shannon pour chaque position
    """ 
    entropy_dict = {}
    base = 4
    entropy = 0

    for cle, valeur in frequence.items():
        entropy_dict[cle] = {}
        for cle2, valeur2 in valeur.items():
            entropy_dict[cle][cle2] = {}
            for cle3, valeur3 in valeur2.items():
                if valeur3 == 0:
                    continue
                else:
                    entropy += - valeur3 * math.log(valeur3, base)
            entropy_dict[cle][cle2] = entropy
            entropy = 0


    return entropy_dict

def save_data_nested_dic(fichier, dictionnaire):
    """sauvegarde un dictionnaire imbriqué dans un fichier texte.

    Parameters
    ----------
    fichier: string
        le nom du fichier dans lequel les sauvegarder
    
    dictionnaire : dictionnary
        dictionnaire contenant les valeurs à sauvegarder
    """
    with open(fichier, "w") as filout:
        for cle, valeur in dictionnaire.items():
            if isinstance(valeur, dict):
                for cle2, valeur2 in valeur.items():
                        filout.write("{} {} {}\n".format(cle, cle2, valeur2))
            if not isinstance(valeur, dict):
                filout.write("{} {}\n".format(cle, valeur))



def main():
    """Le main du programme."""    
    
    fam = input("Quelle famille voulez vous traiter ? ")
    longueur_seq = input("Quelle est la longueur des séquences ? ")
    rounds = []
    fichiers = arguments()

    seq_by_round = {"R00" : 0, "R01" : 0,
    "R02" : 0, "R03" : 0, "R04" : 0, "R05" : 0,
    "R06" : 0, "R07" : 0, "R08" : 0, "R09" : 0, "R10" : 0,
    "R11" : 0, "R12" : 0, "R13" : 0, "R14" : 0, "R15" : 0,
    "R17" : 0, "R18" : 0}

    for fichier in fichiers:
        rounds.append(fichier[16:19])
   
    profil_par_round = {}.fromkeys(set(rounds), {})
    diff_moy_seq = {}.fromkeys(set(rounds), {})


    for fichier in fichiers:
        data = save_data(fichier)
        for cle, valeur in data["24"].items():
            seq_by_round[fichier[16:19]] += valeur
        profil_par_round[fichier[16:19]] = data
        diff_moy_seq[fichier[16:19]] = mutation(data, fichier, longueur_seq, seq_by_round)


    diff_moy_seq = dict(sorted(diff_moy_seq.items(), key = lambda t: t[0]))
    diff_moy = f"difference_moy_entre_deux_seq_fam_{fam}.txt"
    save_data_nested_dic(diff_moy, diff_moy_seq)


    freq_par_round = frequence_par_round(profil_par_round, seq_by_round)
    entropy_dict = calc_shannon_entropy(freq_par_round)
    
    for cle in entropy_dict.keys():    
        entropy_file = f"shannon_entropy_fam_{fam}_{cle}.txt"
        save_data_nested_dic(entropy_file, entropy_dict[cle])
    
    freq_par_round = dict(sorted(freq_par_round.items(), key = lambda t: t[0]))
    
    for k, v in freq_par_round.items():
        for k2, v2 in v.items():
            freq_par_round[k][k2] = dict(sorted(freq_par_round[k][k2].items(), key = lambda t: t[0]))
    
    for cle in freq_par_round.keys():
        freq = f"frequence_par_pos_{cle}_fam_{fam}.txt"
        save_data_nested_dic(freq, freq_par_round[cle])


if __name__ == "__main__":
    main()