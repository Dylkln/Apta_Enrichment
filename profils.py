"""Ce code permet de créer le profil de chaque famille indépendamment du round.

Usage:
------
    python3 profils.py arguments

    arguments: le ou les fichier.s texte à analyser
"""


############ Modules à importer ############


import sys
from tqdm import tqdm
import operator


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


def read_txt_files(fichier):
    """Lit un fichier.

    Parameters
    ----------
    fichier : string
        fichier texte à lire

    Returns
    -------
    data: list
        liste des séquences de la famille
    """ 

    data = []
    
    with open(fichier, "r") as filin:
        lines = filin.readlines()
        for line in lines:
            data.append(line.strip())

    return data


def compte_seq_len(data):
    """dénombre la taille des séquences
    
    Parameters
    ----------
    data: list
        liste des séquences de la famille

    Returns
    -------
    compte : dictionnary
        dictionnaire contenant les différentes tailles des séquences
        et leur nombre d'occurrences au sein de la famille

    """

    compte = {}
    for seq in data:
        cle = len(seq)
        if cle in compte.keys():
            compte[cle] += 1
        else:
            compte[cle] = 0
            compte[cle] += 1

    return compte


def extract_max_from_dict(dictionnaire):
    """donne le maximum d'un dictionnaire.

    Parameters
    ----------
    dictionnaire : dictionnary
        un dictionnaire contenant des valeurs numériques

    returns
    -------    
    seq_len_max : int
        la longueur de la séquence apparaissant le plus de fois dans
        le dictionnaire fournit en argument.
    """

    seq_len_max = max(dictionnaire.items(), key = operator.itemgetter(1))[0]

    return seq_len_max


def create_profils(compte, seq_len_max, data):
    """créer le profil de chaque famille.

    Parameters
    ----------
    compte : dictionnary
        dictionnaire contenant chaque séquence d'une famille à un round donné

    seq_len_max : int
        la longueur de la séquence apparaissant le plus de fois dans
        le dictionnaire fournit en argument.

    data: list
        liste des séquences de la famille

    Returns
    -------
    profils: dictionnary
        dictionnaire contenant le profil de la famille
    """

    profils = {}
    wanted_seq = [data[i] for i in range(len(data)) if len(data[i]) == seq_len_max]
    
    for seq in tqdm(wanted_seq):
        for j, base in enumerate(seq):
            
            if (j + 24) not in profils.keys():
                profils[j + 24] = {}
            
            if base not in profils[j + 24].keys():
                profils[j + 24][base] = 1
            
            else:
                profils[j + 24][base] += 1


    return profils


def save_in_text_file(data, fichier):
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


def save_profil(data, fichier):
    """sauvegarde un dictionnaire imbriqué dans un fichier texte.

    Parameters
    ----------
    data : dictionnary
        dictionnaire contenant les valeurs à sauvegarder
    fichier: string
        le nom du fichier dans lequel les sauvegarder
    """
   
    with open(fichier, "w") as filout:
        for pos, valeur in data.items():
            for base, nombre in valeur.items():
                filout.write(f"{pos} {base} {nombre}\n")


def main():
    """Le main du programme.""" 
    
    fichiers = arguments()
    
    for fichier in fichiers:
        data = read_txt_files(fichier)
        compte = compte_seq_len(data)
        compte = dict(sorted(compte.items(), key = lambda t: t[0]))
        compte_data = "nbr_seq_len_family_{}.txt".format(fichier.strip("family__all_seq.txt"))
        save_in_text_file(compte, compte_data)
        seq_len_max = extract_max_from_dict(compte)
        profils = create_profils(compte, seq_len_max, data)
        profil_file = "profil_famille_{}.txt".format(fichier.strip("family__all_seq.txt"))
        save_profil(profils, profil_file)




if __name__ == "__main__":
    main()