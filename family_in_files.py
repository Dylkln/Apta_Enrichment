"""Ce code permet d'effectuer plusieurs comparaison entre
les familles de séquences représentées dans des fichiers texte, et
les fichiers fasta contenant les séquences.

Usage:
------
    python3 family_in_files.py arguments arguments2

    arguments: fichier texte contenant des séquences nucléiques

    arguments2: fichier fasta contenant des séquences nucléiques
"""


############ Modules à importer ############


import sys
from tqdm import tqdm
import operator


############################################


def arguments(answer):
    """Vérifier le format et le nombre d'arguments renseigné.

    Returns
    -------
    fichiers: liste de tous les fichiers donnés en argument.
    """

    fichiers = []
    if answer != "oui":
        sys.exit("veuillez renseigner des fichiers texte et fasta à lire")
    
    if len(sys.argv) < 3:
        sys.exit("Veuillez renseigner au moins un fichier texte et un fichier fasta à lire")
    
    for i in range(1, len(sys.argv)):
        fichiers.append(str(sys.argv[i]))
    
    return fichiers


def save_data(fichiers):
    """Lit la liste de tous les fichiers renseignés en argument.

    Parameters
    ----------
    fichiers : list
        liste des fichiers à lire

    Returns
    -------
    fichiers_txt: list
        liste des fichiers textes renseignés.

    fichiers_fasta: list
        liste des fichiers fastas renseignés.
    """

    fichiers_txt = []
    fichiers_fasta = []

    for fichier in fichiers:
        if fichier.endswith(".txt"):
            fichiers_txt.append(fichier)
        if fichier.endswith(".fas"):
            fichiers_fasta.append(fichier)

    return fichiers_txt, fichiers_fasta


def read_txt_files(fichier_txt):
    """Lit un fichier texte.

    Parameters
    ----------
    fichiers_txt : list
        liste des fichiers à lire

    Returns
    -------
    data: list
        liste contenant les séquences présentes dans le fichier texte.
    """

    data = []
    
    with open(fichier_txt, "r") as filin:
        lines = filin.readlines()
        for line in lines:
            data.append(line.strip())

    return data


def read_fasta_files(fichier_fasta):
    """Lit un fichier fasta.

    Parameters
    ----------
    fichiers_txt : list
        liste des fichiers à lire

    Returns
    -------
    data: list
        liste contenant les séquences présentes dans le fichier fasta.
    """

    data = []

    with open(fichier_fasta, "r") as filin:
        lines = filin.readlines()
        for line in lines:
            if line.startswith(">"):
                continue
            else:
                data.append(line.strip())

    return data


def extract_common_seq_len(common_seq, taille_round):
    """donne le nombre de séquences communes entre une
    famille et le round, et la fréquence associée.
    
    Parameters
    ----------
    common_seq: list
        liste des séquences communes entre une famille et le round donné

    taille_round: int
        le nombre de séquences présent à ce round
    """

    common_seq_len = 0
    freq_common_seq_len = 0


    for seq in common_seq:
        common_seq_len += 1

    
    freq_common_seq_len = len(common_seq) / taille_round

    
    return common_seq_len, freq_common_seq_len


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


def create_profils(compte, dic_b, common_seq):
    """créer le profil de chaque famille.

    Parameters
    ----------
    compte : dictionnary
        dictionnaire contenant chaque séquence d'une famille à un round donné

    dic_b: dictionnary
        dictionnaire contenant les séquences d'un round donné et leur
        nombre d'occurrences.

    common_seq: list
        list contenant toutes les séquences communes à une famille donnée
        à un round donné

    Returns
    -------
    profils: dictionnary
        dictionnaire contenant le profil de la famille donné à un round donné
    """
    profils = {}

    seq_len_max = extract_max_from_dict(compte)

    wanted_seq = [common_seq[i] for i in range(len(common_seq)) if len(common_seq[i]) == seq_len_max]

    for seq in wanted_seq:
        for j, base in enumerate(seq):
            
            if (j + 24) not in profils.keys():
                profils[j + 24] = {}
            
            if base not in profils[j + 24].keys():
                profils[j + 24][base] = 1
            
            else:
                profils[j + 24][base] += 1
    
    return profils


def create_count_dict(a, b):
    """créer un dictionnaire.

    Parameters
    ----------
    a : list
        liste de séquences

    b: list
        liste de séquences

    Returns
    -------
    dic_b: dictionnary
        dictionnaire contenant les séquences de la liste b et leur
        nombre d'occurrences.

    common_seq: list
        list contenant toutes les séquences communes aux deux listes a et b

    compte: dictionnary
        dictionnaire contenant le le nombre de fois qu'apparait une séquence
        d'une taille donnée
    """  
    
    dic_b = {}.fromkeys(set(b), 0)
    compte = {}
    common_seq = list(set(a) & set(b))

    for seq in b:
        dic_b[seq] += 1

    for seq, nbr in tqdm(dic_b.items()):
        if seq in common_seq:
            tmp = seq.split()
            multiplicateur = nbr - 1
            if multiplicateur != 0:
                common_seq += tmp * multiplicateur

    for seq in common_seq:
        cle = len(seq)
        if cle in compte.keys():
            compte[cle] += 1
        else:
            compte[cle] = 0
            compte[cle] += 1

    return compte, common_seq


def find_common_seq(a, b):
    """trouve les séquences communes entre deux liste de séquences.

    Parameters
    ----------
    a : list
        liste de séquences

    b: list
        liste de séquences

    Returns
    -------
    dic_b: dictionnary
        dictionnaire contenant les séquences de la liste b et leur
        nombre d'occurrences.

    common_seq: list
        list contenant toutes les séquences communes aux deux listes a et b

    taille_round: int
        le nombre de séquences présentes à ce round
    """    

    dic_b = {}.fromkeys(set(b), 0)
    common_seq = list(set(a) & set(b))

    taille_round = len(b)
    for seq in b:
        dic_b[seq] += 1

    for seq, nbr in tqdm(dic_b.items()):
        if seq in common_seq:
            tmp = seq.split()
            multiplicateur = nbr - 1
            if multiplicateur != 0:
                common_seq += tmp * multiplicateur

    return dic_b, common_seq, taille_round


def count_file_seq_in_family(family_dict, round_dict):
    """compte le nombre de séquences présentes à chaque round pour
    chaque famille.

    Parameters
    ----------
    family_dict : dictionnary
        dictionnaire contenant la famille à tester

    round_dict: dictionnary
        dictionnaire contenant les séquences présentes à chaque round

    returns
    -------

    nbr_seq_in_families : dictionnary
        dictionnaire contenant le nombre de séquences présentes par
        famille par round.

    freq_seq_in_families: dictionnary
        dictionnaire contenant la fréquence de chaque famille par round
    
    profils: dictionnary
        dictionnaire contenant les profils de chaque famille par round

    compte: dictionnary
        dictionnaire contenant le nombre d'occurrences de chaque séquence
        de chaque famille par round
    """

    nbr_seq_in_families = {}
    freq_seq_in_families = {}
    profils = {}
    compte = {}

    for cle in round_dict.keys():
        nbr_seq_in_families[cle] = {}
        freq_seq_in_families[cle] = {}
        profils[cle] = {}
        compte[cle] = {}

    for num_fam, seq_list in family_dict.items():
        for Round, seq_list2 in round_dict.items():
            
            compte[Round][num_fam], common_seq = create_count_dict(seq_list, seq_list2)
            dic_b, common_seq, taille_round = find_common_seq(seq_list, seq_list2)
            nbr_seq_in_families[Round][num_fam], freq_seq_in_families[Round][num_fam] = extract_common_seq_len(common_seq, taille_round)
            if common_seq:
                profils[Round][num_fam] = create_profils(compte[Round][num_fam], dic_b, common_seq)
            if not common_seq:
                profils[Round][num_fam] = "Famille absente de ce Round"
    
    return nbr_seq_in_families, freq_seq_in_families, profils, compte


def diversification_mesure(family_dict, round_dict):
    """Mesure la diversification d'une famille au cours des différents rounds.

    Parameters
    ----------
    family_dict : dictionnary
        dictionnaire contenant la famille à tester

    round_dict: dictionnary
        dictionnaire contenant les séquences présentes à chaque round

    Returns
    -------
    div_dic: dictionnary
        dictionnaire contenant toutes les séquences apparaissant pour la
        première fois, à chaque round, plus le nombre de fois ou ces séquences
        sont apparues à chaque round.
    """
    common_seq = {}
    div_dict = {}
    s = set()


    for num_fam, seq_list in tqdm(family_dict.items()):
        for Round, seq_list2 in round_dict.items():
            common_seq[Round] = []
            dic_b, common_seq[Round], taille_round = find_common_seq(seq_list, seq_list2)


    rounds = ("R00", "R01", "R02", "R03", "R04", "R05", "R06", "R07", "R08",
    "R09", "R10", "R11", "R12", "R13", "R14", "R15", "R17", "R18")
    
    for i in rounds:
        div_dict[i] = {}    

    sequences = set()
    
    for r in rounds:
        for seq in common_seq[r]:
            if seq not in sequences:
                if seq not in div_dict[r].keys():
                    div_dict[r][seq] = 0
                div_dict[r][seq] += 1
        for cle, seq in div_dict[r].items():
            sequences.add(seq)


    return div_dict


def save_data_in_txt_file(fichier, fichier2, data):
    """sauvegarde un dictionnaire dans un fichier texte.

    Parameters
    ----------
    fichier: string
        le nom du fichier dans lequel les sauvegarder

    fichiers2: string
        le nom du fichier dans lequel les sauvegarder

    data : dictionnary
        dictionnaire contenant les valeurs à sauvegarder
    """
    somme = {}
    for cle, valeur in data.items():
        somme[cle] = 0    
        for cle2, valeur2 in valeur.items():
            somme[cle] += valeur2
    
    with open(fichier, "w") as filout:
        for cle, valeur in data.items():
            for cle2, valeur2 in valeur.items():    
                filout.write("{} {} {}\n".format(cle, cle2, valeur2))
    
    with open(fichier2, "w") as filout:
        for cle, valeur in data.items():
            filout.write("{} {} {}\n".format(cle, len(valeur), somme[cle]))


def save_data_of_nested_dic(fichier, dictionnaire):
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
                    for cle3, valeur3 in valeur2.items():
                        filout.write("{} {} {} {}\n".format(cle, cle2, cle3, valeur3))
                    filout.write("\n")
            if isinstance(valeur, str):
                filout.write(f"{valeur}\n\n")


def save_data_in_text_file(fichier, dictionnaire):
    """sauvegarde un dictionnaire dans un fichier texte.

    Parameters
    ----------
    fichier: string
        le nom du fichier dans lequel les sauvegarder

    fichiers2: string
        le nom du fichier dans lequel les sauvegarder

    data : dictionnary
        dictionnaire contenant les valeurs à sauvegarder
    """
    
    with open(fichier, "w") as filout:
        for cle, valeur in dictionnaire.items():
            if isinstance(valeur, dict):
                for cle2, valeur2 in valeur.items():
                    filout.write(f"{cle} {cle2} {valeur2}\n")
            
            elif not valeur:
                filout.write(f"la famille {cle} n'est pas présente dans ce Round\n")
            
            elif valeur:
                filout.write(f"{cle} {valeur}\n")


def main():
    """Le main du programme."""

    answer = input("Avez vous renseigné des fichiers texte PUIS des fichiers fasta ? (oui/non) ")
    answer2 = input("Quelle famille voulez vous traiter ? ")

    fichiers = arguments(answer)
    fichiers_txt, fichiers_fasta = save_data(fichiers)

    keys = []
    
    for fichier in tqdm(fichiers_txt):
        if fichier.endswith("_diff_seq.txt"):
            keys.append(int(fichier.strip("family__diff_seq.txt")))
        if fichier.endswith("_all_seq.txt"):
            keys.append(int(fichier.strip("family__all_seq.txt")))
    
    family_dict = {}.fromkeys(set(keys), [])
    
    for fichier in tqdm(fichiers_txt):
        if fichier.endswith("_diff_seq.txt"):
            data = read_txt_files(fichier)
            family_dict[int(fichier.strip("family__diff_seq.txt"))] = list(data)
        if fichier.endswith("_all_seq.txt"):
            data = read_txt_files(fichier)
            family_dict[int(fichier.strip("family__all_seq.txt"))] = list(data)
    
    keys = []
    
    for fichier in tqdm(fichiers_fasta):
        keys.append(fichier[0:3])
    
    round_dict = {}.fromkeys(set(keys), [])
    
    for fichier in tqdm(fichiers_fasta):
        data = read_fasta_files(fichier)
        round_dict[fichier[0:3]] = list(data)


#    diversification_dict = diversification_mesure(family_dict, round_dict)
#    div_file = f"diversification_famille_{answer2}.txt"
#    div_file2 = f"diversification_famille_{answer2}_count.txt"
#    save_data_in_txt_file(div_file, div_file2,  diversification_dict)

    nbr_seq_in_families, freq_seq_in_families, profils, compte = count_file_seq_in_family(family_dict, round_dict)
    
    nbr_seq_in_families = dict(sorted(nbr_seq_in_families.items(), key = lambda t: t[0][0]))
    freq_seq_in_families = dict(sorted(freq_seq_in_families.items(), key = lambda t: t[0][0]))
    compte = dict(sorted(compte.items(), key = lambda t: t[0][0][0]))
    profils = dict(sorted(profils.items(), key = lambda t: t[0][0], reverse = True))

    for cle in compte.keys():
        if fichiers_txt[0].endswith("_diff_seq.txt"):
            fichier_compte = f"seq_count_by_fam_in_round_{cle}_diff_seq.txt"
            save_data_in_text_file(fichier_compte, compte[cle])
        if fichiers_txt[0].endswith("_all_seq.txt"):
            fichier_compte = f"seq_count_by_fam_in_round_{cle}_all_seq.txt"
            save_data_in_text_file(fichier_compte, compte[cle])

    for cle in profils.keys():
        if fichiers_txt[0].endswith("_diff_seq.txt"):
            profil_seq = f"profil_by_fam_in_round_{cle}_diff_seq.txt"
            save_data_of_nested_dic(profil_seq, profils[cle])
        if fichiers_txt[0].endswith("_all_seq.txt"):
            profil_seq = f"profil_by_fam_in_round_{cle}_all_seq.txt"
            save_data_of_nested_dic(profil_seq, profils[cle])

    for cle in nbr_seq_in_families.keys():
        if fichiers_txt[0].endswith("_diff_seq.txt"):
            nbr_seq = f"seq_by_family_in_round_{cle}_diff_seq.txt"
            save_data_in_text_file(nbr_seq, nbr_seq_in_families[cle])
        if fichiers_txt[0].endswith("_all_seq.txt"):
            nbr_seq = f"seq_by_family_in_round_{cle}_all_seq.txt"
            save_data_in_text_file(nbr_seq, nbr_seq_in_families[cle])

    for cle in freq_seq_in_families.keys():
        if fichiers_txt[0].endswith("_diff_seq.txt"):
            freq_seq = f"freq_by_family_in_round_{cle}_diff_seq.txt"
            save_data_in_text_file(freq_seq, freq_seq_in_families[cle])
        if fichiers_txt[0].endswith("_all_seq.txt"):
            freq_seq = f"freq_by_family_in_round_{cle}_all_seq.txt"
            save_data_in_text_file(freq_seq, freq_seq_in_families[cle])

if __name__ == "__main__":
    main()