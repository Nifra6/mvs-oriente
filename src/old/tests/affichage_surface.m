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
surface = "gaussienne_1_bruitee_10";
rayon_voisinage = 4;
ecart_type = 6;
nombre_profondeur_iteration = 5000;
nombre_vues = 4;
grille_pixel = 10;

%% Variables
taille_patch = 2*rayon_voisinage + 1;
path = "../../data/";
fichier_surface = "simulateur_" + surface + "_formate.mat";
load(path+fichier_surface,'nombre_lignes','nombre_colonnes','facteur_k','u_0','v_0','s');
path = "../../result/tests/";
if (ecart_type >= 0)
	fichier_bruite = "__bruite_" + int2str(10) + "__filtre_" + int2str(ecart_type);
else
	fichier_bruite = "";
end
nom_fichier = "Surface_" + surface + "__nb_vues_" + int2str(nombre_vues) + "__patch_" ...
	+ int2str(taille_patch) + "x" + int2str(taille_patch) + "__nb_profondeur_" ...
	+ int2str(nombre_profondeur_iteration) + fichier_bruite + ".mat";
fichier_surface = "simulateur_" + surface + "_formate.mat";
load(path+nom_fichier);

%% Algorithme
% Préparation
X = 1:grille_pixel:nombre_colonnes;
X = (X - u_0) / facteur_k;
Y = 1:grille_pixel:nombre_lignes;
Y = (Y - v_0) / facteur_k;
[X,Y] = meshgrid(X,Y);
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
map_erreur_mvs = zeros(size(X,2),size(Y,2));
map_erreur_mvsm = zeros(size(X,2),size(Y,2));
map_erreur_mvs(ind_1_shrink) = erreur_z_mvs;
map_erreur_mvsm(ind_1_shrink) = erreur_z_mvsm;
min_c_map = min([min(erreur_z_mvs) min(erreur_z_mvsm)]);
max_c_map = max([max(erreur_z_mvs) max(erreur_z_mvsm)]);

complement_titre = ", " + nombre_vues + " vues";
if (ecart_type >= 0)
	complement_titre = complement_titre + " et sigma à " + ecart_type;
end


% Affichage
figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
s = surf(X,Y,-z_estime_mvs(1:grille_pixel:end,1:grille_pixel:end),map_erreur_mvs);
%s.EdgeColor = 'none';
s.CDataMapping = 'scaled';
ax = gca;
ax.CLim = [min_c_map max_c_map];
grid off;
colormap jet;
c = colorbar;
c.Label.String = 'Erreur de profondeurs (en m)';
c.Label.FontSize = 11;
axis equal;
title("Reconstruction MVS" + complement_titre);
view([-90 90]);

% Affichage
figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
s = surf(X,Y,-z_estime_mvsm(1:grille_pixel:end,1:grille_pixel:end),map_erreur_mvsm);
%s.EdgeColor = 'none';
s.CDataMapping = 'scaled';
ax = gca;
ax.CLim = [min_c_map max_c_map];
grid off;
colormap jet;
c = colorbar;
c.Label.String = 'Erreur de profondeurs (en m)';
c.Label.FontSize = 11;
axis equal;
title("Reconstruction MVSm" + complement_titre);
view([-90 90]);
