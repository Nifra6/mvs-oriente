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
liste_surface = ["gaussienne_1", "gaussienne_1_pepper", "gaussienne_2", "sinc_1"];
liste_nombre_vues = [2:9];
liste_rayon_voisinage = 4;
nombre_iteration = 1;
ecart_type = -1;
nombre_profondeur_iteration = 5000;
liste_modifiers = ["o:", "v--", "d-.", "s-"];
reconstruction = "MVSm";
utilisation_profondeur_GT = 0;
utilisation_normale_GT = 1;
utilisation_mediane_normale = 1;

%reconstruction = "all";
%liste_rayon_voisinage = 4;
%liste_modifiers = ["s--","d-"];

%% Modifications noms de fichiers
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
if (utilisation_mediane_normale)
	fichier_mediane = "__normales_medianes";
else
	fichier_mediane = "";
end

%% Variables
nb_rayon_voisinage = size(liste_rayon_voisinage,2);
nb_nombre_vues = size(liste_nombre_vues,2);
nb_surface = size(liste_surface,2);


%% Algorithme
for i_surface = 1:nb_surface
	surface = liste_surface(i_surface);
	erreurs_mvs = [];
	erreurs_mvsm = [];
	% Création de la figure
	figure;
	hold on;
	for i_rayon_voisinage = 1:nb_rayon_voisinage
		rayon_voisinage = liste_rayon_voisinage(i_rayon_voisinage);
		taille_patch = 2*rayon_voisinage + 1;
		for i_nombre_vues = 1:nb_nombre_vues
			nombre_vues = liste_nombre_vues(i_nombre_vues);
			% Chargement des résultats
			nom_fichier = "Surface_" + surface + "__nb_vues_" + int2str(nombre_vues) + "__patch_" ...
				+ int2str(taille_patch) + "x" + int2str(taille_patch) + "__nb_profondeur_" ...
				+ int2str(nombre_profondeur_iteration) + fichier_profondeur_GT ...
				+ fichier_normale_GT + fichier_mediane + ".mat";
			path = "../../result/tests/";
			load(path+nom_fichier);
			% Extraction des données intéressantes
			switch (reconstruction)
				case "MVS"
					erreur_z = erreur_z_mvs;
					erreurs_mvs = [erreurs_mvs , erreur_z];
				case "MVSm"
					erreur_z = erreur_z_mvsm;
					erreurs_mvsm = [erreurs_mvsm , erreur_z];
			end
		end
		% Tracé
		boxplot(erreurs_mvsm);
		%plot(liste_nombre_vues,erreurs_moyennes,'b'+liste_modifiers(i_rayon_voisinage),'LineWidth',1.5);
		%plot(liste_nombre_vues,erreurs_medianes,'g'+liste_modifiers(i_rayon_voisinage),'LineWidth',1.5);
		%legende = [legende , reconstruction + " moyenne en " + int2str(taille_patch) + "x" + int2str(taille_patch) , reconstruction + " médiane en " + int2str(taille_patch) + "x" + int2str(taille_patch) ];
	end
	% Affichage du graphe
	hold off;
	grid on;
	titre = "Erreurs sur la surface " + surface + ", " + int2str(nombre_profondeur_iteration) + " profondeurs testées";
	complement_titre = "Voisinages " + int2str(taille_patch) + "x" + int2str(taille_patch);
	title([titre,complement_titre]);
	xlabel("Nombre de vues utilisées");
	ylabel("Erreurs de profondeur");
	%legend(legende,'Location','best');
	ylim([0,3e-3]);
end
