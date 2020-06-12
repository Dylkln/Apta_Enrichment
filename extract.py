"""Ce code permet d'extraire une famille d'un fichier texte contenant
toutes les familles.

Usage:
------
    python3 extract.py arguments

    arguments: fichier texte contenant des valeurs de toutes les familles

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
        sys.exit("Veuillez renseigner au moins un fichier txt à lire")
    
    for index in range(1, len(sys.argv)):
        if not sys.argv[index].endswith(".txt"):
            sys.exit("Les fichiers renseignés doivent être au format txt")
        fichiers.append(str(sys.argv[index]))
    
    return fichiers


def save_data(fichier, answer):
    """Lit un fichier.

    Parameters
    ----------
    fichier : string
        fichier texte à lire
    
    answer: string
        réponse à l'input lors du lancement du programme

    Returns
    -------
    data: dictionnary
        dictionnaire contenant les valeurs de la famille voulue
    """ 
    data = {}
    
    with open(fichier, "r") as filin:
        lines = filin.readlines()
        lines = filter(str.strip, lines)
        for line in lines:
            if line[0] == answer and line[1] == " ":
                line.strip(line[0:2])
                if line[2:4] not in data.keys():
                    data[line[2:4]] = {}
                    data[line[2:4]][line[5]] = int(line[7:20].strip())
                if line[2:4] in data.keys():
                    data[line[2:4]][line[5]] = int(line[7:20].strip())
            else:
                continue

    return data

def save_dict_in_txt_file(fichier, dictionnaire):
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
            for cle2, valeur2 in valeur.items():
                filout.write("{} {} {}\n".format(cle, cle2, valeur2))


def save_data_in_txt_file(fichier, liste):
    """sauvegarde un dictionnaire dans un fichier texte.

    Parameters
    ----------
    fichier: string
        le nom du fichier dans lequel les sauvegarder
    liste: list
        liste contenant les fichiers créés par ce programme
    """
    with open(fichier, "w") as filout:
        for file in liste:
            filout.write("{}\n".format(file))


def main():

    answer = input("Quelle famille voulez vous extraire ? ")
    fichiers_saved = []
    fichiers = arguments()
    for fichier in fichiers:
        data = save_data(fichier, answer)
        profil_file = f"profil_fam_{answer}_in_{fichier[23:26]}.txt"
        fichiers_saved.append(profil_file)
        save_dict_in_txt_file(profil_file, data)

    liste_file_saved = f"liste_fichiers_profil_fam_{answer}.txt"
    save_data_in_txt_file(liste_file_saved, fichiers_saved)

if __name__ == "__main__":
    main()
