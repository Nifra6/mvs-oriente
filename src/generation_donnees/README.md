# Générations des données

Ce dossier contient l'ensemble des outils pour préparer les données pour les tests de reconstruction. Il comporte deux types de codes : la mise en forme et le bruitage.

## Mise en forme des données

Les codes **conversion_simulateur_orthographique** et **conversion_simulateur_perspectif** sont les premiers codes à exécuter en cas d'obtentions de nouvelles données depuis le simulateur. Ces scripts vont convertir le fichier **simulateur.mat** contenu dans le dossier **data/[orthographique|perspectif]** en un fichier **simulateur_formate.mat** qui est utilisable par les restes des codes.
Il y a deux options :
- filtrage\_points est un booléen permettant de retirer les pixels qui ne sont pas présents dans l'ensemble des images du jeu de test, afin de renforcer la confiance dans les statistiques générées (les occultations ne sont pas considérés dans ces codes). En cas de crash pour mémoire insuffisante, y compris lorsque c'est le seul processus actif, vous pouvez augmentez la valeur de grille\_pixel, un paramètre limitant les pixels utilisés. Ne pas utiliser si la calotte, ou tout autre scène non continue est utilisée.
- masque_calotte est un booléen permettant de travailler avec les images qui comportent une scène limitée. Dans le cadre d'une calotte de sphère, seuls les pixels appartenant à la calotte seront étudiés.

La section *Mise en forme des données* de ces codes vous renseigneront sur les variables contenus dans les fichiers formatés.

## Bruitage

Deux types de bruitages additifs sont possibles : bruit blanc uniforme ou bruit poivre et sel. La force des bruits peut être réglée et ils sont appliqués sur le fichier **simulateur_formate.mat** courant.

## Utilisation

Petit récapitulatif de l'utilisation standard avec des conseils de nomenclatures :
1. Sauvegarder votre jeu de donnée dans le répertoire data correspondant à son type de caméra (orthographique ou perspectif), sous le nom **simulateur_votreNom.mat**.
2. Le copier sous le nom **simulateur.mat**.
3. Utiliser le code conversion\_simulateur correspondant au type de votre caméra.
4. Optionnel : appliquer le bruit de votre choix en utilisant le code de bruitage voulu.
5. Copier le fichier **simulateur_formate.mat** sous le nom **simulateur_votreNom_formate.mat**. Le nom du jeu de donnée *votreNom* sera la chaîne de caractère à utiliser comme paramètres pour ce jeu de donnée.
