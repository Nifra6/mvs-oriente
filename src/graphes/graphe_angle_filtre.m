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
liste_surface = ["gaussienne_1_bruitee_2", "gaussienne_1_pepper_bruitee_2", "gaussienne_2_bruitee_2", "sinc_1_bruitee_2"];
nombre_vues = 2;
rayon_voisinage = 4;
nombre_iteration = 1;
liste_ecart_type_grad = [0:10];
liste_ecart_type_I = [0:6];
nombre_profondeur_iteration = 5000;
ecart_type_grad_lock = 0;
ecart_type_I_lock = 6;
utilisation_profondeur_GT = 0;
utilisation_normale_GT = 1;

%% Variables
nb_surface = size(liste_surface,2);
nb_ecart_type_grad = size(liste_ecart_type_grad,2);
nb_ecart_type_I = size(liste_ecart_type_I,2);
taille_patch = 2*rayon_voisinage+1;
erreurs_angulaires_moyennes = zeros(nb_ecart_type_grad,1);
erreurs_angulaires_medianes = zeros(nb_ecart_type_grad,1);
erreurs_z_moy = zeros(nb_ecart_type_grad,1);
erreurs_z_med = zeros(nb_ecart_type_grad,1);

%% Algorithme
for i_surface = 1:nb_surface
	surface = liste_surface(i_surface);
	for i_ecart_type_grad = 1:nb_ecart_type_grad
		ecart_type_grad = liste_ecart_type_grad(i_ecart_type_grad);
		% Chargement des résultats
		if (utilisation_profondeur_GT)
			fichier_profondeur_GT = "__profondeurs_GT";
		else
			fichier_profondeur_GT = "";
		end
		if (utilisation_normale_GT)
			fichier_normale_GT = "__normales_GT";
		else
			fichier_normale_GT = "";
		end
		if (ecart_type_I_lock > 1)
			nom_fichier = "Surface_" + surface + "__nb_vues_" + int2str(nombre_vues) + "__patch_" ...
				+ int2str(taille_patch) + "x" + int2str(taille_patch) + "__nb_profondeur_" ...
				+ int2str(nombre_profondeur_iteration) + "__bruite_" + int2str(valeur_bruitage) ...
				+ "__filtre_I_" + int2str(ecart_type_I_lock) + "__filtre_grad_" ...
				+ int2str(ecart_type_grad) + fichier_profondeur_GT + fichier_normale_GT + ".mat";
		else
			nom_fichier = "Surface_" + surface + "__nb_vues_" + int2str(nombre_vues) + "__patch_" ...
				+ int2str(taille_patch) + "x" + int2str(taille_patch) + "__nb_profondeur_" ...
				+ int2str(nombre_profondeur_iteration) + "__bruite_" + int2str(valeur_bruitage) ...
				+ "__filtre_" + int2str(ecart_type_grad) + fichier_profondeur_GT ...
				+ fichier_normale_GT + ".mat";
		end
		path = "../../result/tests/";
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
	ylim([4.9,45]);

	if (ecart_type_grad_lock >= 0)
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
		titre = ["Erreurs de profondeurs sur la surface " + surface ; "Avec sigma = " + int2str(ecart_type_grad_lock) + " fixé pour les gradients"];
		title(titre);
		xlabel('Écart type')
		ylabel('Erreurs de profondeurs')
		legend('Erreurs moyennes','Erreurs médianes','Location','best');
		ylim([2e-3,0.1]);
	else	
		figure;
		plot(liste_ecart_type_grad,erreurs_z_moy,'bv-','LineWidth',1.5);
		hold on;
		plot(liste_ecart_type_grad,erreurs_z_med,'mv-','LineWidth',1.5);
		hold off;
		grid on;
		titre = "Erreurs de profondeurs sur la surface " + surface;
		title(titre);
		xlabel('Écart type')
		ylabel('Erreurs de profondeurs')
		legend('Erreurs moyennes','Erreurs médianes','Location','best');
		ylim([2e-3,0.1]);
	end
end

