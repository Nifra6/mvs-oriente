# Lancement de la stéréoscopie multi-vues

Ce dossier contient trois fichier principaux et a pour but de lancer une estimation des cartes de profondeurs ou des cartes de normales. Le seul fichier à utiliser est **lancement_test.m**. Le coeur de l'implémentation se trouve dans les fichiers **mvs.m** et **mvs_modifie.m**. Ces deux fichiers font le travail effectif de recherche exhaustive de profondeur, d'estimation des normales, des reprojections.

## Utilisation

Pour réaliser un test ou une série de test, il faut modifier le fichier **lancement_test.m** pour utiliser les paramètres voulus. Les variables nommées *liste_quelquechose* admettent plusieurs valeurs pour lancer une série de test.
- **liste_surface** doit contenir les noms des jeux de données à utiliser. Ils doivent être enregistré dans le répertoire **data** correspondant au type de la caméra sous le nom **simulateur_nom_formate.mat**.
- **liste_rayon_voisinage** doit contenir les moitiés de tailles de patchs souhaitées. Ainsi, pour utiliser des vosinages 9x9, il suffit d'entrer la valeur 4. Des patchs plus grands ralentissent grandement le processus.
- **liste_ecart_type** I ou grad pemettent de régler l'écart type du filtre utilisés sur les images (pour le critère photométrique) ou les gradients (pour filtrer le gradient utilisé à l'estimation des normales). Si aucun filtrage n'est appliqué, il vaut mieux laisser ces valeurs à -1.
- **liste_nombre_vues** indique combien d'images doit être utilisées pour la reconstrution. Attention, cela compte l'image de référence. Pour utiliser 4 images témoins, il faut renseigner la valeur 5. Un grand nombre d'image ralentit le processus.
- **liste_nombre_profondeur_iteration** permet de chosir le nombre d'échantillon prit sur le rayon pour chque pixel. Augmenter ce nombre augmente linéairement le temps d'exécution, c'est la boucle principale des algorithmes.
- **valeur_bruitage** permet de spécifier la valeur de bruitage utilisée sur les données préalablement bruitées.
- **filtrage** permet d'indiquer si un filtrage doit être appliqué aux images et/ou gradients. Actuellement il semble que son usage soit supplanté par la valeur de **liste_ecart_type_[I|grad]**, cela doit être réparé.
- **grille_pixels** indiquent la distance spatiale entre 2 pixels étudiés. Une valeur à 2 (minimale sinon le code plante) permet donc de considérer un pixel sur 4 (un pixel toutes les 2 lignes et 2 colonnes).
- **utilisation_profondeur_GT** permet de fixer la profondeurs des pixels à la vérité terrain enregistrées. Ce paramètre permet d'étudier l'estimation des normales. Les tests ne comporteront plus qu'une seule itération (pas de recherche exhaustive).
- **utilisation_normale_GT** permet d'utiliser les normales vérité terrains pour la reprojections des patchs de voisinages. Ce paramètre permet d'étudier le gain engendré par la bonne orientations des patchs selon la normale à la surface.

Une fois le paramétrage fait, il suffit d'exécuter le fichier **lancement_test.m** avant d'aller se faire un café en salle de pose, lire un article, aller s'acheter le repas de midi et manger, dormir, lancer pédantix/cémantix, continuer de coder, lire un autre article.
