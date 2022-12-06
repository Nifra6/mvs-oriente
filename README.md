# Amélioration de la stéréoscopie multi-vues par photoclinométrie

Ce repo contient les différents codes MATLAB élaborés au cours de mon stage de fin d'études.

## Structure

Le dossier **src** contient les différents codes utilisables, regroupés par catégorie :
1. Le dossier **calotte** contient des codes permettant de simuler une calotte de sphère et de faire du rendu orthographique
2. Le dossier **generation_donnees** permet de convertir des données issues du simulateur lambertien disponible [ici](https://github.com/bbrument/lambertianRendering_v1). Il contient également des codes permettant de bruiter les images.
3. Les dossiers **orthographique** et **perspectif** permettent de lancer une série de reconstruction 3D en utilisant la stéréoscopie multi-vues, dans le cadre de projection orthographique ou perspective. De multiples paramètres sont disponibles et la reconstruction 3D s'effectue à la fois en homographie fronto-parrallèle et en homographie orientée par la normale à la surface.
4. Les dossiers **graphes** et **graphes_perspectif** permettent d'anaylser les résultats obtenus par recosntruction grâce à différents types de graphes.
5. Le dossier **toolbox** contient quelques fonctions utiles pour les autres codes.
6. Le dossier **old** contient une majorité de codes qui ne sont plus d'actualités mais qui peuvent être utiles.

## Utilisation

Pour estimer les profondeurs et analyser les résultats, il suffit de se référer aux README.md dans les dossiers suivants :
1. Le dossier **generation_donnees**
2. Le dossier **perspectif** 
3. Le dossier **graphes_perspectifs**
