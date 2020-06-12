"""Ce code permet d'extraire la taille des familles créées.

Usage:
------
    python3 nbr_seq_in_files.py arguments

    arguments: fichier texte contenant des séquences nucléiques d'une famille

"""


############ Modules à importer ############


import sys


############################################


def arguments():
    """Vérifier le format et le nombre d'arguments renseigné

    Returns
    -------
    fichiers: liste de tous les fichiers donnés en argument
    """  
    fichiers = []
    
    if len(sys.argv) < 2:
        sys.exit("Veuillez renseigner au moins un fichier texte à lire")
    
    for index in range(1, len(sys.argv)):
        if not sys.argv[index].endswith(".txt"):
            sys.exit("Les fichiers renseignés doivent être au format texte")
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
    nbr_ligne: int
        nombre de séquences dans le fichier
    """    
    nbr_ligne = 0
    
    with open(fichier, "r") as filin:
        lines = filin.readlines()
        for line in lines:
            nbr_ligne += 1
    
    return nbr_ligne

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


def main():
    """Le main du programme."""
    nbr_seq_in_files = {}
    fichiers = arguments()

#    for fichier in fichiers:
#        if len(sys.argv) > 2:
#            text_file = "nbr_seq_par_round.txt"
#            nbr_seq = save_data(fichier)
#            nbr_seq_in_files[fichier[1:3]] = nbr_seq
#        else:
#            text_file = "nbr_seq_all_round.txt"
#            nbr_seq = save_data(fichier)
#            nbr_seq_in_files["nbr_seq_supp_to_1000_all_rounds"] = nbr_seq

    for fichier in fichiers:
        text_file = "taille_des_familles.txt"
        nbr_seq = save_data(fichier)
        if fichier.endswith("_diff_seq.txt"):
            nbr_seq_in_files[fichier.strip("family__diff_seq.txt")] = nbr_seq
        if fichier.endswith("_all_seq.txt"):
            nbr_seq_in_files[fichier.strip("family__all_seq.txt")] = nbr_seq

    save_in_text_file(nbr_seq_in_files, text_file)


if __name__ == "__main__":
    main()