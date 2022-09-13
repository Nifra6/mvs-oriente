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

%reconstruction = "all";
%liste_rayon_voisinage = 4;
%liste_modifiers = ["s--","d-"];

%% Variables
nb_rayon_voisinage = size(liste_rayon_voisinage,2);
nb_nombre_vues = size(liste_nombre_vues,2);
nb_surface = size(liste_surface,2);

erreurs_moyennes = zeros(1,nb_nombre_vues);
erreurs_medianes = zeros(1,nb_nombre_vues);
erreurs_moyennes_mvs = zeros(1,nb_nombre_vues);
erreurs_medianes_mvs = zeros(1,nb_nombre_vues);
erreurs_moyennes_mvsm = zeros(1,nb_nombre_vues);
erreurs_medianes_mvsm = zeros(1,nb_nombre_vues);
legende = [];

%% Algorithme
for i_surface = 1:nb_surface
	surface = liste_surface(i_surface);
	% Création de la figure
	figure;
	hold on;
	for i_rayon_voisinage = 1:nb_rayon_voisinage
		rayon_voisinage = liste_rayon_voisinage(i_rayon_voisinage);
		taille_patch = 2*rayon_voisinage + 1;
		for i_nombre_vues = 1:nb_nombre_vues
			nombre_vues = liste_nombre_vues(i_nombre_vues);
			% Chargement des résultats
			nom_fichier = "Surface_" + surface + "__nb_vues_" + int2str(nombre_vues) + "__patch_" + int2str(taille_patch) + "x" ...
				+ int2str(taille_patch) + "__nb_profondeur_" + int2str(nombre_profondeur_iteration) + ".mat";
			path = "../../result/tests/orthographique/";
			load(path+nom_fichier);
			% Extraction des données intéressantes
			switch (reconstruction)
				case "MVS"
					erreur_z = erreur_z_mvs;
					erreurs_moyennes(i_nombre_vues) = mean(erreur_z);
					erreurs_medianes(i_nombre_vues) = median(erreur_z);
				case "MVSm"
					erreur_z = erreur_z_mvsm;
					erreurs_moyennes(i_nombre_vues) = mean(erreur_z);
					erreurs_medianes(i_nombre_vues) = median(erreur_z);
				case "all"
					erreurs_moyennes_mvs(i_nombre_vues) = mean(erreur_z_mvs);
					erreurs_medianes_mvs(i_nombre_vues) = median(erreur_z_mvs);
					erreurs_moyennes_mvsm(i_nombre_vues) = mean(erreur_z_mvsm);
					erreurs_medianes_mvsm(i_nombre_vues) = median(erreur_z_mvsm);
			end
		end
		% Tracé
		switch (reconstruction)
			case 'all'
				plot(liste_nombre_vues,erreurs_moyennes_mvs,'b'+liste_modifiers(1),'LineWidth',1.5);
				plot(liste_nombre_vues,erreurs_medianes_mvs,'b'+liste_modifiers(2),'LineWidth',1.5);
				plot(liste_nombre_vues,erreurs_moyennes_mvsm,'m'+liste_modifiers(1),'LineWidth',1.5);
				plot(liste_nombre_vues,erreurs_medianes_mvsm,'m'+liste_modifiers(2),'LineWidth',1.5);
				legende = [legende , "MVS moyenne en " + int2str(taille_patch) + "x" + int2str(taille_patch) , "MVS médiane en " + int2str(taille_patch) + "x" + int2str(taille_patch) , "MVSm moyenne en " + int2str(taille_patch) + "x" + int2str(taille_patch) , "MVSm médiane en " + int2str(taille_patch) + "x" + int2str(taille_patch)];
			otherwise
				plot(liste_nombre_vues,erreurs_moyennes,'b'+liste_modifiers(i_rayon_voisinage),'LineWidth',1.5);
				plot(liste_nombre_vues,erreurs_medianes,'g'+liste_modifiers(i_rayon_voisinage),'LineWidth',1.5);
				legende = [legende , reconstruction + " moyenne en " + int2str(taille_patch) + "x" + int2str(taille_patch) , reconstruction + " médiane en " + int2str(taille_patch) + "x" + int2str(taille_patch) ];
		end
	end
	% Affichage du graphe
	hold off;
	grid on;
	titre = "Erreurs sur la surface " + surface + ", " + int2str(nombre_profondeur_iteration) + " profondeurs testées";
	title(titre);
	xlabel("Nombre de vues utilisées");
	ylabel("Erreurs de profondeur");
	legend(legende,'Location','best');
	ylim([0,3e-3]);
end
