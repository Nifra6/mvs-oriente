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
surface = "gaussienne_1_bruitee_" + int2str(valeur_bruitage);
surface = "gaussienne_decentree";
nombre_vues = 9;
rayon_voisinage = 4;
ecart_type_grad = -5;
ecart_type_I = -2.5;
filtrage = 0;
nombre_profondeur_iteration = 1000;
utilisation_profondeur_GT = 0;
utilisation_normale_GT = 1;
grille_pixel = 4;
mesure = "median";
mesure = "all";
utilisation_mediane_normale = 1;

%% Variables
taille_patch = 2*rayon_voisinage + 1;
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

if (ecart_type_grad >= 0 & filtrage)
	if (ecart_type_I >= 0)
		fichier_bruite = "__bruite_" + int2str(valeur_bruitage) + "__filtre_I_" ...
			+ num2str(ecart_type_I) + "__filtre_grad_" + num2str(ecart_type_grad);
	else
		fichier_bruite = "__bruite_" + int2str(valeur_bruitage) + "__filtre_" + num2str(ecart_type_grad);
	end
else
	fichier_bruite = "";
end


nom_fichier = "Surface_" + surface + "__nb_vues_" + int2str(nombre_vues) + "__patch_" ...
	+ int2str(taille_patch) + "x" + int2str(taille_patch) + "__nb_profondeur_" ...
	+ int2str(nombre_profondeur_iteration) + fichier_bruite + fichier_profondeur_GT ...
	+ fichier_normale_GT + fichier_mediane + ".mat";
path = "../../result/tests/perspectif/";
load(path+nom_fichier);
grille_pixel = grille_pixels;

% Préparation de la reconstruction
path_data = '../../data/perspectif/simulateur_';
load(path_data + surface + '_formate.mat','nb_lignes','nb_colonnes','f','u_0','v_0','s','N','masque','R','t','K');
masque_1 = masque(:,:,1); clear masque;
masque_1_shrink = masque_1(1:grille_pixel:end,1:grille_pixel:end);
ind_1_shrink = find(masque_1_shrink);
[i_k,j_k] = find(masque_1);
ind_1 = sub2ind([nb_lignes nb_colonnes],i_k,j_k);
indices_grille = (mod(i_k,grille_pixel) == 1) & (mod(j_k,grille_pixel) == 1);
ind_1 = ind_1(find(indices_grille));
N_1 = N(:,:,:,1); clear N;
normales_GT = [N_1(ind_1)' ; N_1(ind_1 + nb_lignes*nb_colonnes)' ; N_1(ind_1 + 2*nb_lignes*nb_colonnes)'];
X_o = 1:grille_pixel:nb_colonnes;
Y_o = 1:grille_pixel:nb_lignes;
X_o = X_o - u_0;
Y_o = Y_o - v_0;
[X_o,Y_o] = meshgrid(X_o,Y_o);

% Préparation des erreurs
map_erreur_mvs = zeros(size(X_o,2),size(Y_o,2));
map_erreur_mvsm = zeros(size(X_o,2),size(Y_o,2));
size(ind_1_shrink)
size(erreur_z_mvs)
map_erreur_mvs(ind_1_shrink) = erreur_z_mvs;
size(X_o)
size(ind_1_shrink)
size(erreur_z_mvsm)
map_erreur_mvsm(ind_1_shrink) = erreur_z_mvsm;
min_c_map = min([min(erreur_z_mvs) min(erreur_z_mvsm)]);
max_c_map = max([max(erreur_z_mvs) max(erreur_z_mvsm)]);

map_erreur_angles_GT = zeros(size(X_o,2),size(Y_o,2));
erreurs_angles_GT = angle_normale(normales_GT,normales_mvsm);
map_erreur_angles_GT(ind_1_shrink) = erreurs_angles_GT;

complement_titre = ", " + nombre_vues + " vues";
if (ecart_type_grad >= 0 & filtrage)
	if (ecart_type_I >= 0)
		complement_titre = complement_titre + ", sigma_grad à " + ecart_type_grad + " et sigma_I à " + ecart_type_I;
	else
		complement_titre = complement_titre + " et sigma à " + ecart_type_grad;
	end
else
	complement_titre = complement_titre + ", sans filtrage";
end
if (utilisation_profondeur_GT)
	complement_titre = complement_titre + ", profondeurs VT";
end
if (utilisation_normale_GT)
	complement_titre = complement_titre + ", normales VT";
end
if (utilisation_mediane_normale)
	complement_titre = complement_titre + ", normales médianes";
end


%% Analyse des résultats
normales_fronto = zeros(size(normales_mvsm));
normales_fronto(3,:) = -1;
angles_mvs = angle_normale(normales_fronto,normales_mvs);
angles_mvsm = angle_normale(normales_fronto,normales_mvsm);
angles_GT = angle_normale(normales_fronto, normales_GT);
color_map_value_GT = zeros(floor(nb_lignes/grille_pixel),floor(nb_colonnes/grille_pixel));
color_map_value_mvs = zeros(floor(nb_lignes/grille_pixel),floor(nb_colonnes/grille_pixel));
color_map_value_mvsm = zeros(floor(nb_lignes/grille_pixel),floor(nb_colonnes/grille_pixel));

zones_angles = 0:10:180;
zones_angles = 0:10:50;
nombre_zones = size(zones_angles,2) - 1;
erreurs_mvs_moy = zeros(1,nombre_zones);
erreurs_mvs_med = zeros(1,nombre_zones);
erreurs_mvsm_moy = zeros(1,nombre_zones);
erreurs_mvsm_med = zeros(1,nombre_zones);
nombre_points_zones = zeros(1,nombre_zones);
label_zones = [];
for k = 1:nombre_zones
	indices_GT = find(zones_angles(k) <= angles_GT & angles_GT < zones_angles(k+1));
	color_map_value_GT(ind_1_shrink(indices_GT')) = k;
	nombre_points_zones(k) = length(indices_GT);
	erreurs_mvs_moy(k) = transpose(mean(erreur_z_mvs(indices_GT)));
	erreurs_mvs_med(k) = transpose(median(erreur_z_mvs(indices_GT)));
	erreurs_mvsm_moy(k) = transpose(mean(erreur_z_mvsm(indices_GT)));
	erreurs_mvsm_med(k) = transpose(median(erreur_z_mvsm(indices_GT)));
	label_zones = [ label_zones , zones_angles(k+1) ];
end

%% Affichage
% Histogramme
if (~utilisation_profondeur_GT)
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
	%title(["Erreurs sur la surface " + surface ; "avec " + int2str(nombre_profondeur_iteration) + " échantillons" + complement_titre],'interpreter','none');

	% Affichage des nombres
	xtips1 = b(1).XEndPoints;
	ytips1 = b(1).YEndPoints;
	labels1 = string(nombre_points_zones);
	text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
		'VerticalAlignment','bottom')
end

% Préparation de la reconstruction
Z = z_estime_mvs(1:grille_pixel:end,1:grille_pixel:end);
[nb_l,nb_c] = size(X_o);
R_inv = R(:,:,1)'
t_inv = - R_inv * t(:,1);
p = repmat(Z(:)',3,1) .* (inv(K) * [X_o(:)' ; Y_o(:)' ; ones(1,nb_c*nb_l)]);
P = R_inv * p  + t_inv;
X = P(1,:);
Y = P(2,:);
Z = P(3,:);
X = reshape(X,nb_l,nb_c);
Y = reshape(Y,nb_l,nb_c);
Z = reshape(Z,nb_l,nb_c);


% Affichage de la reconstruction
figure('Name','Relief MVS','Position',[0,0,0.33*L,0.5*H]);
subplot(2,1,1);
sl = surf(X,Y,Z,color_map_value_GT);
sl.EdgeColor = 'none';
sl.CDataMapping = 'scaled';
ax = gca;
ax.CLim = [0 5];
grid off;
colormap 'jet';
axis equal;
%title("Relief MVS",'interpreter','none');
title("Orientation des normales",'interpreter','none');
colorbar;
view([-90 90])

if (~utilisation_profondeur_GT)
	subplot(2,1,2);
	s = surf(X,Y,-z_estime_mvs(1:grille_pixel:end,1:grille_pixel:end),map_erreur_mvs);
	s.EdgeColor = 'none';
	s.CDataMapping = 'scaled';
	ax = gca;
	ax.CLim = [min_c_map max_c_map];
	grid off;
	colormap jet;
	c = colorbar;
	c.Label.String = 'Erreur de profondeurs (en m)';
	c.Label.FontSize = 11;
	c.Location = "east";
	c.AxisLocation = "out";
	axis equal;
	title("Reconstruction MVS" + complement_titre,'interpreter','none');
	view([-90 90]);
end

% Préparation de la reconstruction
Z = z_estime_mvsm(1:grille_pixel:end,1:grille_pixel:end);
[nb_l,nb_c] = size(X_o);
R_inv = R(:,:,1)'
t_inv = - R_inv * t(:,1);
p = repmat(Z(:)',3,1) .* (inv(K) * [X_o(:)' ; Y_o(:)' ; ones(1,nb_c*nb_l)]);
P = R_inv * p  + t_inv;
X = P(1,:);
Y = P(2,:);
Z = P(3,:);
X = reshape(X,nb_l,nb_c);
Y = reshape(Y,nb_l,nb_c);
Z = reshape(Z,nb_l,nb_c);


% Affichage de la reconstruction
figure('Name','Relief MVS modifié','Position',[0,0,0.33*L,0.5*H]);
subplot(2,1,1);
sl = surf(X,Y,Z,color_map_value_GT);
sl.EdgeColor = 'none';
sl.CDataMapping = 'scaled';
ax = gca;
ax.CLim = [0 5];
grid off;
colormap 'jet';
axis equal;
title("Relief MVS modifié",'interpreter','none');
%view([-90 90]);

if (~utilisation_profondeur_GT)
	subplot(2,1,2);
	s = surf(X,Y,-z_estime_mvsm(1:grille_pixel:end,1:grille_pixel:end),map_erreur_mvsm);
	s.EdgeColor = 'none';
	s.CDataMapping = 'scaled';
	ax = gca;
	ax.CLim = [min_c_map max_c_map];
	grid off;
	colormap jet;
	c = colorbar;
	c.Label.String = 'Erreur de profondeurs (en m)';
	c.Label.FontSize = 11;
	c.Location = "east";
	c.AxisLocation = "out";
	axis equal;
	title("Reconstruction MVSm" + complement_titre,'interpreter','none');
	view([-90 90]);
end


diff_erreurs = erreurs_mvs_med - erreurs_mvsm_med
pourcentages_diff_erreurs = 100 * (erreurs_mvs_med - erreurs_mvsm_med) ./ erreurs_mvs_med
nombre_points_zones




% Affichage de la reconstruction
figure('Name','Différence angulaire','Position',[0,0,0.33*L,0.5*H]);
imagesc(map_erreur_angles_GT')
%sl = surf(X,Y,-z_estime_mvsm(1:grille_pixel:end,1:grille_pixel:end),map_erreur_angles_GT);
%sl.EdgeColor = 'none';
%sl.CDataMapping = 'scaled';
%ax = gca;
%ax.CLim = [0 5];
%view([-90 90]);
grid off;
colormap 'jet';
colorbar
axis equal;
%title("Différence angulaire entre normales GT et normales estimées",'interpreter','none');


map_erreur_fronto_GT = zeros(size(X,2),size(Y,2));
erreurs_fronto_GT = angle_normale(normales_GT,normales_fronto);
map_erreur_fronto_GT(ind_1_shrink) = erreurs_fronto_GT;

map_erreur_fronto_estim = zeros(size(X,2),size(Y,2));
erreurs_fronto_estim = angle_normale(normales_mvsm,normales_fronto);
map_erreur_fronto_estim(ind_1_shrink) = erreurs_fronto_estim;

min_map = min([min(map_erreur_fronto_GT,[],'all'),min(map_erreur_fronto_estim,[],'all')]);
max_map = max([max(map_erreur_fronto_GT,[],'all'),max(map_erreur_fronto_estim,[],'all')]);

figure('Name','Différence angulaire','Position',[0,0,0.33*L,0.5*H]);
sl = surf(X,Y,-z_estime_mvsm(1:grille_pixel:end,1:grille_pixel:end),map_erreur_fronto_GT);
sl.EdgeColor = 'none';
sl.CDataMapping = 'scaled';
ax = gca;
ax.CLim = [min_map max_map];
grid off;
colormap 'jet';
colorbar
axis equal;
title("Différence angulaire entre normales GT et normales frontoparallèles",'interpreter','none');
view([-90 90]);


figure('Name','Différence angulaire','Position',[0,0,0.33*L,0.5*H]);
sl = surf(X,Y,-z_estime_mvsm(1:grille_pixel:end,1:grille_pixel:end),map_erreur_fronto_estim);
sl.EdgeColor = 'none';
sl.CDataMapping = 'scaled';
ax = gca;
ax.CLim = [min_map max_map];
grid off;
colormap 'jet';
colorbar
axis equal;
title("Différence angulaire entre normales estimées et normales frontoparallèles",'interpreter','none');
view([-90 90]);
