%% Trucs de Matlab
% Clear
clear;
close all;
% Paramètres d'affichage
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);
% Imports de fonctions utiles
addpath(genpath('../toolbox/'));

%% Paramètres
valeur_bruitage = 2;
liste_surface = ["gaussienne_1_bruitee_" + int2str(valeur_bruitage), "gaussienne_1_pepper_bruitee_" ...
	+ int2str(valeur_bruitage), "gaussienne_2_bruitee_" + int2str(valeur_bruitage), "sinc_1_bruitee_" ...
	+ int2str(valeur_bruitage)];
nombre_vues = 2;
rayon_voisinage = 4;
nombre_iteration = 1;
liste_ecart_type_I = [0:0.5:3];
nombre_profondeur_iteration = 5000;
ecart_type_grad_lock = 0;
utilisation_profondeur_GT = 0;
utilisation_normale_GT = 1;

%% Variables
nb_surface = size(liste_surface,2);
nb_ecart_type_I = size(liste_ecart_type_I,2);
taille_patch = 2*rayon_voisinage+1;
erreurs_z_moy = zeros(nb_ecart_type_I,1);
erreurs_z_med = zeros(nb_ecart_type_I,1);

if (utilisation_profondeur_GT)
	fichier_profondeur_GT = "__profondeurs_GT";
	liste_ecart_type_I = 0;
	nb_ecart_type_I = 1;
else
	fichier_profondeur_GT = "";
end
if (utilisation_normale_GT)
	fichier_normale_GT = "__normales_GT";
	ecart_type_grad_lock = 0;
else
	fichier_normale_GT = "";
end

%% Algorithme
for i_surface = 1:nb_surface
	surface = liste_surface(i_surface);
	erreurs_z_moy = zeros(nb_ecart_type_I,1);
	erreurs_z_med = zeros(nb_ecart_type_I,1);
	for i_ecart_type_I = 1:nb_ecart_type_I
		ecart_type_I = liste_ecart_type_I(i_ecart_type_I);
		% Chargement des résultats
		nom_fichier = "Surface_" + surface + "__nb_vues_" + int2str(nombre_vues) + "__patch_" ...
			+ int2str(taille_patch) + "x" + int2str(taille_patch) + "__nb_profondeur_" ...
			+ int2str(nombre_profondeur_iteration) + "__bruite_" + int2str(valeur_bruitage) ...
			+ "__filtre_I_" + int2str(ecart_type_I) + "__filtre_grad_" ...
			+ int2str(ecart_type_grad_lock) + fichier_profondeur_GT + fichier_normale_GT + ".mat";
		path = "../../result/tests/";
		load(path+nom_fichier);
		% Extraction des données intéressantes
		erreurs_z_moy(i_ecart_type_I) = mean(erreur_z_mvsm);
		erreurs_z_med(i_ecart_type_I) = median(erreur_z_mvsm);
	end
	figure;
	plot(liste_ecart_type_I,erreurs_z_moy,'bv-','LineWidth',1.5);
	hold on;
	plot(liste_ecart_type_I,erreurs_z_med,'mv-','LineWidth',1.5);
	hold off;
	grid on;
	if (utilisation_normale_GT)
		sous_titre = "Avec les normales vérité terrain"; 
	else
		sous_titre = "Avec sigma = " + int2str(ecart_type_grad_lock) + " fixé pour les gradients"; 
	end
	titre = ["Erreurs de profondeurs sur la surface " + surface ; sous_titre];
	title(titre);
	xlabel('Écart type')
	ylabel('Erreurs de profondeurs')
	legend('Erreurs moyennes','Erreurs médianes','Location','best');
	%ylim([2e-3,0.1]);
end

