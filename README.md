Projet de Stage : Etudes de familles de séquences
==================================================

## Description ##

Ces scripts permettent d’utiliser les données d’un SELEX déjà réalisé et
ayant conduit à l’identification de familles d’aptamères et de parasites,
pour essayer d’identifier si elles ont des profils d’évolution différents.


## Installation ##

Pour ces scripts, les modules suivants déjà installés dans la version standard de Python sont nécessaires :

- os
- sys
- math
- operator
- copy

Il est également nécessaire d'avoir installé les modules pythons suivants :

- python-Levenshtein (module permettant entre autre de calculer la distance de Levenshtein
- tqdm (ajoute une barre de progression du programme)

Aucun autre module n'est nécessaire pour faire fonctionner les scripts.

Les scripts utilisent Python 3.7 qui est la version de Python disponible par défaut avec "Anaconda3" ou "Miniconda3"

## Utilisation ##

Le script python "create_family.py" est utilisé pour créer un nombre de famille défini à partir 
des séquences présentes dans des fichiers fasta.

Le script python "nbr_seq_in_files.py" est utilisé pour extraire la taille des familles créées.

Le script python "family_in_files.py" est utilisé pour effectuer des comparaison entre les séquences des familles
et les séquences présentes dans les fichiers fasta.

Le script python "entropy.py" est utilisé pour calculer l'entropie de Shannon.

Le script python "profils.py" est utilisé pour créer les profils des différentes familles créées.

Le script python "mutation.py" permet de calculer l'entropie de Shannon à chaque position
des séquences d'une famille à chaque round.

Le script python "extract.py" permet d'extraire une famille d'un fichier texte contenant
toutes les familles.



## Support ##

En cas de problème, merci de me contacter au mail suivant:

    klein.dylan@outlook.com

