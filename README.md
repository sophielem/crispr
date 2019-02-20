# Pré-calcul des noeuds : intersection et union

Les pré-calculs sont stockés sous format *pickle* dans un dossier *nodes* contenant les sous-dossiers : *inter* et *union*.
Le dictionnaire sauvé contient deux clefs principales: **metadata** et **data**.
Data contient un dictionnaire contenant les séquences sgRNA comme clef, suivi de l'organisme, de la référence du génome afin de stocker une liste de coordonnées.
Metadata contient une liste de fichiers, qui recensent les fils n+1 du noeud étudié. Ainsi, il est possible de reconstruire le noeud directement à partir des noeuds et feuilles fils de ce noeud. 

## Intersection  
Si le noeud étudié ne contient que des feuilles, un simple lancement du script *allgenomes.py* est effectué. De ce fait, pour le génome le plus petit, son fichier pickle va être chargé, contenant le dictionnaire contenant toutes les séquences sgRNA, et des *bowtie* vont être effectués afin de trouver les séquences sgRNA communes à tous les autres génomes.

Si le noeud contient aussi des noeuds, alors la première chose effectuée est le calcul de l'intersection des différents noeuds. Pour cela, les fichiers *pickle* contenant les séquences sgRNA, stockées sous la forme d'un dictionnaire en tant que clef, qui sont communes aux noeuds sont chargées et seules les séquences sgRNA, donc les clefs, communes aux différents dictionnaires sont gardées.
Si des feuilles sont ajoutées, alors le dictionnaire de séquence obtenue lors de l'intersection des noeuds sert de référence pour la recherche des séquences à l'aide de *bowtie*.

## Union
Que ce soit pour les noeuds ou les feuilles, le même princpe est appliqué. Le fichier *pickle* contenant le dictionnaire contenant les séquences sgRNA de la feuille ou de l'union s'il s'agit d'un noeud est chargé. De là, une recherche pour avoir une liste de clefs unique, donc de séquences sgRNA unique, est effectuée. Ainsi, la redondance des clefs est évitée.
