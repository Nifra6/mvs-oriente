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
liste_surface = ["gaussienne_1","gaussienne_1_pepper","gaussienne_2","sinc_1"];
rayon_voisinage = 4;
ecart_type = -1;
nombre_profondeur_iteration = 5000;
nombre_vues = 4;
grille_pixel = 10;

%% Variables
zones_angles = 0:10:60;
nombre_zones = size(zones_angles,2) - 1;
erreurs_mvs_all_all = zeros(40000,nombre_zones);
erreurs_mvsm_all_all = zeros(40000,nombre_zones);
nombre_points = zeros(1,nombre_zones);

nb_surface = size(liste_surface,2);
for i_surface = 1:nb_surface
	% Sauvegarde des résultats
	nom_fichier = "donnees_histo_angles_" + surface + ".mat";
	path = "../../result/intermed/";
	load(path+nom_fichier,"erreurs_mvs_all","erreurs_mvsm_all","nombre_points_zones");

	% Rassemblement
	for k = 1:nombre_zones
		erreurs_mvs_all_all(1+nombre_points(k):nombre_points(k)+nombre_points_zones(k),k) = erreurs_mvs_all;
		erreurs_mvsm_all_all(1+nombre_points(k):nombre_points(k)+nombre_points_zones(k),k) = erreurs_mvsm_all;
		nombre_points(k) = nombre_points(k) + nombre_points_zones(k);
	end
end

erreurs_mvs = zeros(1,nombre_zones);
erreurs_mvsm = zeros(1,nombre_zones);
for k = 1:nb_surface

end









taille_patch = 2*rayon_voisinage + 1;
if (ecart_type >= 0)
	fichier_bruite = "__bruite_" + int2str(10) + "__filtre_" + int2str(ecart_type);
else
	fichier_bruite = "";
end
nom_fichier = "Surface_" + surface + "__nb_vues_" + int2str(nombre_vues) + "__patch_" + int2str(taille_patch) + "x" ...
	+ int2str(taille_patch) + "__nb_profondeur_" + int2str(nombre_profondeur_iteration) + fichier_bruite + ".mat";
path = "../../result/tests/";
load(path+nom_fichier);

% Préparation de la reconstruction
load('../../data/simulateur_' + surface + '_formate.mat','nombre_lignes','nombre_colonnes','facteur_k','u_0','v_0','s','N','masque');
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
X = 1:grille_pixel:nombre_colonnes;
X = (X - u_0) / facteur_k;
Y = 1:grille_pixel:nombre_lignes;
Y = (Y - v_0) / facteur_k;
[X,Y] = meshgrid(X,Y);


%% Analyse des résultats
normales_fronto = zeros(size(normales_mvsm));
normales_fronto(3,:) = -1;
angles_mvs = angle_normale(normales_fronto,normales_mvs);
angles_mvsm = angle_normale(normales_fronto,normales_mvsm);
angles_GT = angle_normale(normales_fronto, normales_GT);
color_map_value_GT = zeros(nombre_lignes/10,nombre_colonnes/10);
color_map_value_mvs = zeros(nombre_lignes/10,nombre_colonnes/10);
color_map_value_mvsm = zeros(nombre_lignes/10,nombre_colonnes/10);

zones_angles = 0:10:180;
zones_angles = 0:10:60;
nombre_zones = size(zones_angles,2) - 1;
erreurs_mvs = zeros(1,nombre_zones);
erreurs_mvsm = zeros(1,nombre_zones);
erreurs_mvs_all = zeros(10000,nombre_zones);
erreurs_mvsm_all = zeros(10000,nombre_zones);
nombre_points_zones = zeros(1,nombre_zones);
label_zones = [];
for k = 1:nombre_zones
	indices_GT = find(zones_angles(k) <= angles_GT & angles_GT < zones_angles(k+1));
	color_map_value_GT(ind_1_shrink(indices_GT')) = k;
	nombre_points_zones(k) = length(indices_GT);
	erreurs_mvs(k) = transpose(median(erreur_z_mvs(indices_GT)));
	erreurs_mvsm(k) = transpose(median(erreur_z_mvsm(indices_GT)));
	erreurs_mvs_all(1:nombre_points_zones(k),k) = transpose(erreur_z_mvs(indices_GT));
	erreurs_mvsm_all(1:nombre_points_zones(k),k) = transpose(erreur_z_mvsm(indices_GT));
	label_zones = [ label_zones , zones_angles(k+1) ];
end

%% Affichage
% Histogramme
figure
b = bar(label_zones,[erreurs_mvs ; erreurs_mvsm]);
legend('MVS','MVS modifié')
xlabel('Angles des normales avec la direction de la caméra de référence')
ylabel('Erreurs de profondeurs médianes')

% Affichage des nombres
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(nombre_points_zones);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
	'VerticalAlignment','bottom')

% Préparation de la reconstruction
X = 1:grille_pixel:nombre_colonnes;
X = (X - u_0) / facteur_k;
Y = 1:grille_pixel:nombre_lignes;
Y = (Y - v_0) / facteur_k;
[X,Y] = meshgrid(X,Y);

% Affichage de la reconstruction
figure('Name','Relief MVS','Position',[0,0,0.33*L,0.5*H]);
sl = surf(X,Y,-z_estime_mvs(1:grille_pixel:end,1:grille_pixel:end),color_map_value_GT);
%sl.EdgeColor = 'none';
grid off;
%colormap gray;
axis equal;
title("Relief MVS");

% Affichage de la reconstruction
figure('Name','Relief MVS modifié','Position',[0,0,0.33*L,0.5*H]);
sl = surf(X,Y,-z_estime_mvsm(1:grille_pixel:end,1:grille_pixel:end),color_map_value_GT);
%sl.EdgeColor = 'none';
grid off;
%colormap gray;
axis equal;
title("Relief MVS modifié");

erreurs_mvs - erreurs_mvsm
nombre_points_zones


% Sauvegarde des résultats
nom_fichier = "donnees_histo_angles_" + surface + ".mat";
path = "../../result/intermed/";
save(path+nom_fichier,"erreurs_mvs_all","erreurs_mvsm_all","nombre_points_zones");
