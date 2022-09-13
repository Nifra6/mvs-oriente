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
	+ int2str(valeur_bruitage)]; set_surface = "";
liste_surface = ["gaussienne_1","gaussienne_1_pepper","gaussienne_2","sinc_1"]; set_surface = "";
%liste_surface = ["gaussienne_1_bis","gaussienne_1_pepper_bis","gaussienne_2_bis","sinc_1_bis"]; set_surface = " bis";
nombre_vues = 2;
rayon_voisinage = 4;
ecart_type_grad = -5;
ecart_type_I = -1; 
filtrage = 0;
nombre_profondeur_iteration = 5000;
utilisation_profondeur_GT = 0;
utilisation_normale_GT = 1;
utilisation_mediane_normale = 0;
grille_pixel = 10;
mesure = "median";
mesure = "all";

%% Variables
nb_surface = size(liste_surface,2);
taille_patch = 2*rayon_voisinage + 1;

complement_titre = nombre_vues + " vues";
if (ecart_type_grad >= 0 & filtrage)
	if (ecart_type_I >= 0)
		complement_titre = complement_titre + ", sigma_grad à " + ecart_type_grad + " et sigma_I à " + ecart_type_I;
	else
		complement_titre = complement_titre + " et sigma à " + ecart_type_grad;
	end
end
if (utilisation_profondeur_GT)
	fichier_profondeur_GT = "__profondeurs_GT";
	ecart_type_I = 0;
	complement_titre = complement_titre + ", profondeurs VT";
else
	fichier_profondeur_GT = "";
end
if (utilisation_normale_GT)
	fichier_normale_GT = "__normales_GT";
	ecart_type_grad = 0;
	complement_titre = complement_titre + ", normales VT";
else
	fichier_normale_GT = "";
end
if (utilisation_mediane_normale)
	fichier_mediane = "__normales_medianes";
else
	fichier_mediane = "";
end



%% Analyse des résultats
zones_angles = 0:10:60;
nombre_zones = size(zones_angles,2) - 1;
erreurs_mvs_moy = zeros(1,nombre_zones);
erreurs_mvs_med = zeros(1,nombre_zones);
erreurs_mvsm_moy = zeros(1,nombre_zones);
erreurs_mvsm_med = zeros(1,nombre_zones);
nombre_points_zones = zeros(1,nombre_zones);
label_zones = [];

for k = 1:nombre_zones
	erreurs_mvs = [];
	erreurs_mvsm = [];

	for i_surface = 1:nb_surface

		% Préparation de la surface
		surface = liste_surface(i_surface);
		if (ecart_type_grad >= 0 & filtrage)
			if (ecart_type_I >= 0)
				fichier_bruite = "__bruite_" + int2str(valeur_bruitage) + "__filtre_I_" ...
					+ num2str(ecart_type_I) + "__filtre_grad_" + num2str(ecart_type_grad);
			else
				fichier_bruite = "__bruite_" + num2str(valeur_bruitage) + "__filtre_" ...
					+ int2str(ecart_type_grad);
			end
		else
			fichier_bruite = "";
		end
		nom_fichier = "Surface_" + surface + "__nb_vues_" + int2str(nombre_vues) + "__patch_" ...
			+ int2str(taille_patch) + "x" + int2str(taille_patch) + "__nb_profondeur_" ...
			+ int2str(nombre_profondeur_iteration) + fichier_bruite  + fichier_profondeur_GT ...
			+ fichier_normale_GT  + fichier_mediane + ".mat";
		path = "../../result/tests/orthographique/";
		load(path+nom_fichier);

		% Préparation des normales GT
		load('../../data/orthographique/simulateur_' + surface + '_formate.mat','nombre_lignes','nombre_colonnes','facteur_k','u_0','v_0','s','N','masque');
		masque_1 = masque(:,:,1); clear masque;
		grille_pixel = 10;
		masque_1_shrink = masque_1(1:grille_pixel:end,1:grille_pixel:end);
		ind_1_shrink = find(masque_1_shrink);
		[i_k,j_k] = find(masque_1);
		ind_1 = sub2ind([nombre_lignes nombre_colonnes],i_k,j_k);
		indices_grille = (mod(i_k,grille_pixel) == 1) & (mod(j_k,grille_pixel) == 1);
		ind_1 = ind_1(find(indices_grille));
		N_1 = N(:,:,:,1); clear N;
		normales_GT = [N_1(ind_1)' ; N_1(ind_1 + nombre_lignes*nombre_colonnes)' ; N_1(ind_1 + 2*nombre_lignes*nombre_colonnes)'];

		% Récupération des angles
		normales_fronto = zeros(size(normales_mvsm));
		normales_fronto(3,:) = -1;
		angles_mvs = angle_normale(normales_fronto,normales_mvs);
		angles_mvsm = angle_normale(normales_fronto,normales_mvsm);
		angles_GT = angle_normale(normales_fronto, normales_GT);
		indices_GT = find(zones_angles(k) <= angles_GT & angles_GT < zones_angles(k+1));
		nombre_points_zones(k) = nombre_points_zones(k) + length(indices_GT);
		erreurs_mvs = [erreurs_mvs ; erreur_z_mvs(indices_GT)];
		erreurs_mvsm = [erreurs_mvsm ; erreur_z_mvsm(indices_GT)];

	end
	erreurs_mvs_moy(k) = transpose(mean(erreurs_mvs));
	erreurs_mvs_med(k) = transpose(median(erreurs_mvs));
	erreurs_mvsm_moy(k) = transpose(mean(erreurs_mvsm));
	erreurs_mvsm_med(k) = transpose(median(erreurs_mvsm));
	label_zones = [ label_zones , zones_angles(k+1) ];
end

%% Affichage
% Histogramme
figure
if (mesure == "median")
	b = bar(label_zones,[erreurs_mvs_med ; erreurs_mvsm_med]);
	legend('MVS','MVS modifié','Location','best')
	xlabel('Angles des normales avec la direction de la caméra de référence')
	ylabel('Erreurs de profondeurs médianes')
else
	b = bar(label_zones,[erreurs_mvs_moy ; erreurs_mvsm_moy ; erreurs_mvs_med ; erreurs_mvsm_med]);
	legend('Moyenne MVS','Moyenne MVS modifié','Médiane MVS','Médiane MVS modifié','Location','best')
	xlabel('Angles des normales avec la direction de la caméra de référence')
	ylabel('Erreurs de profondeurs')
end
title(["Erreurs sur toutes les surfaces" + set_surface + " avec " + int2str(nombre_profondeur_iteration) + " échantillons"; complement_titre],'interpreter','none');

% Affichage des nombres
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(nombre_points_zones);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
	'VerticalAlignment','bottom')

diff_erreurs = erreurs_mvs_med - erreurs_mvsm_med
pourcentages_diff_erreurs = 100 * (erreurs_mvs_med - erreurs_mvsm_med) ./ erreurs_mvs_med
nombre_points_zones

