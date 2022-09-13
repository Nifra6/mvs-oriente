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
surface = "gaussienne_1";
nombre_vues = 5;
liste_rayon_voisinage = 4;
nombre_iteration = 1;
ecart_type = -1;
liste_nombre_profondeur_iteration = [10, 50, 100, 500, 1000, 5000, 10000];
utilisation_profondeur_GT = 0;
utilisation_normale_GT = 0;
utilisation_mediane_normale = 1;

%% Variables
nb_rayon_voisinage = size(liste_rayon_voisinage,2);
nb_nombre_profondeur = size(liste_nombre_profondeur_iteration,2);
erreurs_moyennes_mvs = zeros(1,nb_nombre_profondeur);
erreurs_moyennes_mvsm = zeros(1,nb_nombre_profondeur);
erreurs_medianes_mvs = zeros(1,nb_nombre_profondeur);
erreurs_medianes_mvsm = zeros(1,nb_nombre_profondeur);

if (utilisation_profondeur_GT)
	fichier_profondeur_GT = "__profondeurs_GT";
	ecart_type_I = 0;
else
	fichier_profondeur_GT = "";
end
if (utilisation_normale_GT)
	fichier_normale_GT = "__normales_GT";
	ecart_type_grad = 0;
else
	fichier_normale_GT = "";
end
if (utilisation_mediane_normale)
	fichier_mediane = "__normales_medianes";
else
	fichier_mediane = "";
end


%% Algorithme
for i_rayon_voisinage = 1:nb_rayon_voisinage
	rayon_voisinage = liste_rayon_voisinage(i_rayon_voisinage);
	taille_patch = 2*rayon_voisinage + 1;
	for i_nombre_profondeur = 1:nb_nombre_profondeur
		nombre_profondeur_iteration = liste_nombre_profondeur_iteration(i_nombre_profondeur);
		% Chargement des résultats
		nom_fichier = "Surface_" + surface + "__nb_vues_" + int2str(nombre_vues) + "__patch_" ...
			+ int2str(taille_patch) + "x" + int2str(taille_patch) + "__nb_profondeur_" ...
			+ int2str(nombre_profondeur_iteration) + fichier_profondeur_GT + fichier_normale_GT ...
			+ fichier_mediane + ".mat";
		path = "../../result/tests/orthographique/";
		load(path+nom_fichier);
		% Extraction des données intéressantes
		erreurs_moyennes_mvs(i_nombre_profondeur) = mean(erreur_z_mvs);
		erreurs_moyennes_mvsm(i_nombre_profondeur) = mean(erreur_z_mvsm);
		erreurs_medianes_mvs(i_nombre_profondeur) = median(erreur_z_mvs);
		erreurs_medianes_mvsm(i_nombre_profondeur) = median(erreur_z_mvsm);
	end
	% Affichage du graphe
	figure;
	loglog(liste_nombre_profondeur_iteration,erreurs_moyennes_mvs,'bs-','LineWidth',1.5);
	hold on;
	loglog(liste_nombre_profondeur_iteration,erreurs_moyennes_mvsm,'mo-','LineWidth',1.5);
	loglog(liste_nombre_profondeur_iteration,erreurs_medianes_mvs,'gs-','LineWidth',1.5);
	loglog(liste_nombre_profondeur_iteration,erreurs_medianes_mvsm,'co-','LineWidth',1.5);
	hold off;
	grid on;
	titre = "Erreurs sur la surface " + surface;
	complement_titre = int2str(nombre_vues) + " vues, patchs " + int2str(taille_patch) + "x" + int2str(taille_patch);
	title([titre; complement_titre],'interpreter','none');
	xlabel("Nombre d'échantillons de profondeur");
	ylabel("Erreurs de profondeur");
	legend("MVS moyenne", "MVSM moyenne", "MVS médiane", "MVSM médiane",'Location','best');
	%ylim([7e-4,3e-2]);
end
