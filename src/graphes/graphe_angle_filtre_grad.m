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
valeur_bruitage = 4;
liste_surface = ["gaussienne_1_bruitee_" + int2str(valeur_bruitage), "gaussienne_1_pepper_bruitee_" ...
	+ int2str(valeur_bruitage), "gaussienne_2_bruitee_" + int2str(valeur_bruitage), "sinc_1_bruitee_" ...
	+ int2str(valeur_bruitage)];
%liste_surface = ["gaussienne_1_bruitee_" + int2str(valeur_bruitage)];
nombre_vues = 5;
rayon_voisinage = 4;
nombre_iteration = 1;
liste_ecart_type_grad = [0:10];
nombre_profondeur_iteration = 5000;
ecart_type_I_lock = 1;
utilisation_profondeur_GT = 1;
utilisation_normale_GT = 0;
utilisation_mediane_normale = 1;

%% Variables
nb_surface = size(liste_surface,2);
nb_ecart_type_grad = size(liste_ecart_type_grad,2);
taille_patch = 2*rayon_voisinage+1;
erreurs_angulaires_moyennes = zeros(nb_ecart_type_grad,1);
erreurs_angulaires_medianes = zeros(nb_ecart_type_grad,1);

if (utilisation_profondeur_GT)
	fichier_profondeur_GT = "__profondeurs_GT";
	ecart_type_I_lock = 0;
else
	fichier_profondeur_GT = "";
end
if (utilisation_normale_GT)
	fichier_normale_GT = "__normales_GT";
	liste_ecart_type_grad = 0;
	nb_ecart_type_grad = 1;
else
	fichier_normale_GT = "";
end
if (utilisation_mediane_normale)
	fichier_mediane = "__normales_medianes";
else
	fichier_mediane = "";
end

%% Algorithme
for i_surface = 1:nb_surface
	surface = liste_surface(i_surface);
	for i_ecart_type_grad = 1:nb_ecart_type_grad
		ecart_type_grad = liste_ecart_type_grad(i_ecart_type_grad);
		% Chargement des résultats
		nom_fichier = "Surface_" + surface + "__nb_vues_" + int2str(nombre_vues) + "__patch_" ...
			+ int2str(taille_patch) + "x" + int2str(taille_patch) + "__nb_profondeur_" ...
			+ int2str(nombre_profondeur_iteration) + "__bruite_" + int2str(valeur_bruitage) ...
			+ "__filtre_I_" + num2str(ecart_type_I_lock) + "__filtre_grad_" ...
			+ num2str(ecart_type_grad) + fichier_profondeur_GT + fichier_normale_GT ...
			+ fichier_mediane + ".mat";
		path = "../../result/tests/orthographique/";
		load(path+nom_fichier);
		% Extraction des données intéressantes
		erreurs_angulaires_moyennes(i_ecart_type_grad) = erreur_angle_moy;
		erreurs_angulaires_medianes(i_ecart_type_grad) = erreur_angle_med;
		erreurs_z_moy(i_ecart_type_grad) = mean(erreur_z_mvsm);
		erreurs_z_med(i_ecart_type_grad) = median(erreur_z_mvsm);
	end
	% Affichage du graphe
	figure;
	plot(liste_ecart_type_grad,erreurs_angulaires_moyennes,'bs-','LineWidth',1.5);
	hold on;
	plot(liste_ecart_type_grad,erreurs_angulaires_medianes,'gs-','LineWidth',1.5);
	hold off;
	grid on;
	titre = "Erreurs angulaires sur la surface " + surface;
	title(titre);
	xlabel('Écart type')
	ylabel('Erreurs angulaires')
	legend('Erreurs moyennes','Erreurs médianes','Location','best');
	%ylim([4.9,45]);
end
