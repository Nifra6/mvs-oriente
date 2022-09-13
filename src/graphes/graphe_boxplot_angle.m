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
surface = "sinc_1_bruitee_10";
rayon_voisinage = 4;
ecart_type_grad = 6;
ecart_type_I = 1; 
nombre_profondeur_iteration = 5000;
nombre_vues = 4;
grille_pixel = 10;

%% Variables
taille_patch = 2*rayon_voisinage + 1;
if (ecart_type_grad >= 0)
	if (ecart_type_I >= 0)
		fichier_bruite = "__bruite_" + int2str(10) + "__filtre_I_" + int2str(ecart_type_I) ...
			+ "__filtre_grad_" + int2str(ecart_type_grad);
	else
		fichier_bruite = "__bruite_" + int2str(10) + "__filtre_" + int2str(ecart_type_grad);
	end
else
	fichier_bruite = "";
end
nom_fichier = "Surface_" + surface + "__nb_vues_" + int2str(nombre_vues) + "__patch_" ...
	+ int2str(taille_patch) + "x" + int2str(taille_patch) + "__nb_profondeur_" ...
	+ int2str(nombre_profondeur_iteration) + fichier_bruite + ".mat";
path = "../../result/tests/orthographique/";
load(path+nom_fichier);

% Préparation de la reconstruction
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
X = 1:grille_pixel:nombre_colonnes;
X = (X - u_0) / facteur_k;
Y = 1:grille_pixel:nombre_lignes;
Y = (Y - v_0) / facteur_k;
[X,Y] = meshgrid(X,Y);

% Préparation des erreurs
map_erreur_mvs = zeros(size(X,2),size(Y,2));
map_erreur_mvsm = zeros(size(X,2),size(Y,2));
map_erreur_mvs(ind_1_shrink) = erreur_z_mvs;
map_erreur_mvsm(ind_1_shrink) = erreur_z_mvsm;
min_c_map = min([min(erreur_z_mvs) min(erreur_z_mvsm)]);
max_c_map = max([max(erreur_z_mvs) max(erreur_z_mvsm)]);

complement_titre = ", " + nombre_vues + " vues";
if (ecart_type_grad >= 0)
	if (ecart_type_I >= 0)
		complement_titre = complement_titre + ", sigma_grad à " + ecart_type_grad + " et sigma_I à " + ecart_type_I;
		else
	complement_titre = complement_titre + " et sigma à " + ecart_type_grad;
	end
end

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
erreurs_mvs = [];
erreurs_mvsm = [];
erreurs_inter = [];
label_points_zones = [];
label_zones = [];
for k = 1:nombre_zones
	indices_GT = find(zones_angles(k) <= angles_GT & angles_GT < zones_angles(k+1));
	color_map_value_GT(ind_1_shrink(indices_GT')) = k;
	nombre_points_zones(k) = length(indices_GT);
	erreurs_mvs = [erreurs_mvs ; erreur_z_mvs(indices_GT)];
	erreurs_mvsm = [erreurs_mvsm ; erreur_z_mvsm(indices_GT)];
	erreurs_inter = [erreurs_inter ; zeros(size(indices_GT'))];
	label_points_zones = [label_points_zones ; 3 * k * ones(length(indices_GT),1)];
	label_zones = [ label_zones , zones_angles(k+1) ];
end


%% Affichage
% Histogramme
figure
b = boxplot([erreurs_mvs ; erreurs_mvsm ; erreurs_inter], [label_points_zones ; label_points_zones + 1 ; label_points_zones + 2]);
%legend('MVS','MVS modifié','Location','best')
xlabel('Angles des normales avec la direction de la caméra de référence')
ylabel('Erreurs de profondeurs médianes')
title(["Erreurs sur la surface " + surface ; "avec " + int2str(nombre_vues) + " vues, " + int2str(nombre_profondeur_iteration) + " échantillons" + complement_titre]);

% Affichage des nombres
%xtips1 = b(1).XEndPoints;
%ytips1 = b(1).YEndPoints;
%labels1 = string(nombre_points_zones);
%text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%	'VerticalAlignment','bottom')

% Préparation de la reconstruction
X = 1:grille_pixel:nombre_colonnes;
X = (X - u_0) / facteur_k;
Y = 1:grille_pixel:nombre_lignes;
Y = (Y - v_0) / facteur_k;
[X,Y] = meshgrid(X,Y);

% Affichage de la reconstruction
figure('Name','Relief MVS','Position',[0,0,0.33*L,0.5*H]);
subplot(2,1,1);
sl = surf(X,Y,-z_estime_mvs(1:grille_pixel:end,1:grille_pixel:end),color_map_value_GT);
%sl.EdgeColor = 'none';
sl.CDataMapping = 'scaled';
ax = gca;
ax.CLim = [0 5];
grid off;
colormap 'jet';
axis equal;
title("Relief MVS");
view([-90 90])

subplot(2,1,2);
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

% Affichage de la reconstruction
figure('Name','Relief MVS modifié','Position',[0,0,0.33*L,0.5*H]);
subplot(2,1,1);
sl = surf(X,Y,-z_estime_mvsm(1:grille_pixel:end,1:grille_pixel:end),color_map_value_GT);
%sl.EdgeColor = 'none';
sl.CDataMapping = 'scaled';
ax = gca;
ax.CLim = [0 5];
grid off;
colormap 'jet';
axis equal;
title("Relief MVS modifié");
view([-90 90]);

subplot(2,1,2);
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


%diff_erreurs = erreurs_mvs_med - erreurs_mvsm_med
%pourcentages_diff_erreurs = 100 * (erreurs_mvs_med - erreurs_mvsm_med) ./ erreurs_mvs_med
nombre_points_zones

