# Analyse des résultats et tracé de figures

Ce dossier contient les différents scripts permettant de tracer les figures, graphes et représentations 3D nécessaires à l'analyse des résultats obtenus avec la stéréoscopie multi-vues orientée. La pluspart des script s'utilisent de la même manière mais produisent des résultats différents.

## Utilisations

Les scripts contiennent une partie *Paramètres* qui permettent de sélectionner les données à utiliser. Les noms de paramètres et leur signification sont identiques aux paramètres de lancement de tests, disponible dans les dossiers **orthographique** et **perspectif**. Si le test n'a pas été lancé, une erreur sera affichée, indiquant le nom du fichier qui aurait du être ouvert. Parfois le fichier que vous souhaitez étudié est bien présent mais il suffit qu'un paramètre soit mal renseigné pour que le nom de fichier soit différent et corrspondent donc à un test non lancé. Pensez donc à vérifier dans vos dossiers result/tests/[perspectif|orthographique] quels fichiers existent.

## Les différents graphes

Voici la liste des scripts usuels et les figures qu'ils génèrent :
- **graphe_histogramme_angle.m** est le script le plus fourni et qui utilise le plus de données. Son but premier est d'afficher les histogrammes d'erreurs de profondeurs en fonction de l'angle entre l'axe optique et la normale à la surface. Des représentations 3D des reconstructions sont également affichées, indiquant les zones angulaires correspondant à l'histogramme (par tranche de 10°). Deux représentations 3D sont également fournies, qui utilisent les profondeurs estimées par la stéréoscopie multi-vues standard et notre propositions, colorisées en fonction de l'erreur de profondeur commmise.
Les erreurs d'estimations de normales sont également fournies, sous forme d'erreur d'angle entre la normale vérité terrain et la normale estimée.  Le script affiche également sous représentation 3D les différences angulaires entre normales vérité terrains et axe optique, ainsi qu'entre normales estimées et axe optique.
- **graphe_estimation_normale.m** permet d'étudier l'estimation des normales lorsque la profondeur vérité terrain est utilisée. Il affiche une représentation 3D avec les zones angulaires ainsi que la carte des erreurs angulaires d'estimations des normales. Les représentations 3D des erreurs entre normales vérité terrain et axe optique et entre normales estimées et axe optique sont également affichées.
- **graphe_angle_filtre_I.m** permettent d'étudier l'évolution de l'erreur de profondeur en fonction des paramètres des filtres utilisés pour débruiter les images, avec utilisation des normales vérité terrains.
- **graphe_angle_filtre_grad.m** permettent d'étudier l'évolution de l'erreur angulaire d'estimation des normales en fonction des paramètres des filtres utilisés pour débruiter les gradients des images, avec utilisation des profondeurs vérité terrains.
- **graphe_boxlot_angle.m** affiche sous forme de boîte à moustache l'évolution des erreurs de profondeurs en fonction des zones angulaires.
- **graphe_boxlot_vues.m** affiche sous forme de boîte à moustache l'évolution des erreurs de profondeurs en fonction du nombre de vue utilisées.
- **graphe_representation_normale.m** permet d'afficher en 3D la reconstruction avec notre stéréoscopie multi-vues orientées, accompagnée des normales vérité terrains et normales estimées.
