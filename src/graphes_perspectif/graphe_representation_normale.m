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
valeur_bruitage = 0;
surface = "gaussienne_decentree_corrige";
surface = "boite_simple_scale8_GT";
%surface = "plan_peppers_21flou_16bit_cote_corrige";
%surface = "calotte"
%surface = "sinus";
nombre_vues = 5;
rayon_voisinage = 4;
ecart_type_grad = 10;
ecart_type_I = 0;
filtrage = 0;
nombre_profondeur_iteration = 5000;
utilisation_profondeur_GT = 1;
utilisation_normale_GT = 0;
mesure = "median";
mesure = "all";

%% Variables
taille_patch = 2*rayon_voisinage + 1;
if (utilisation_profondeur_GT)
	fichier_profondeur_GT = "__profondeurs_GT";
	fichier_profondeur = "";
else
	fichier_profondeur_GT = "";
	fichier_profondeur = "__nb_profondeur_" + int2str(nombre_profondeur_iteration);
end
ecart_type_I = 0;
fichier_normale_GT = "";

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
	+ int2str(taille_patch) + "x" + int2str(taille_patch) + fichier_profondeur ...
	+ fichier_bruite + fichier_profondeur_GT + fichier_normale_GT + ".mat";
path = "../../result/tests/perspectif/";
load(path+nom_fichier);
grille_pixel = grille_pixels;

% Préparation de la reconstruction
path_data = '../../data/perspectif/simulateur_';
load(path_data + surface + '_formate.mat','nb_lignes','nb_colonnes','f','u_0','v_0','s','N','masque','R','t','K','I');


figure
I_truc = masque(:,:,1).*I(:,:,1);
imshow(I_truc);

masque_1 = masque(:,:,1); clear masque;
masque_1(1:rayon_voisinage,:) = 0;
masque_1(end-rayon_voisinage:end,:) = 0;
masque_1(:,1:rayon_voisinage) = 0;
masque_1(:,end-rayon_voisinage:end) = 0;
masque_1_shrink = masque_1(1:grille_pixel:end,1:grille_pixel:end);
ind_1_shrink = find(masque_1_shrink);
[i_k,j_k] = find(masque_1);
ind_1 = sub2ind([nb_lignes nb_colonnes],i_k,j_k);
indices_grille = (mod(i_k,grille_pixel) == 1) & (mod(j_k,grille_pixel) == 1);
ind_1 = ind_1(find(indices_grille));
N_1 = N(:,:,:,1); clear N;
normales_GT = [N_1(ind_1)' ; N_1(ind_1 + nb_lignes*nb_colonnes)' ; N_1(ind_1 + 2*nb_lignes*nb_colonnes)'];
f = K(1,1); u_0 = K(1,3); v_0 = K(2,3);
X_o = 1:grille_pixel:nb_colonnes;
Y_o = 1:grille_pixel:nb_lignes;
X_o = X_o - u_0;
Y_o = Y_o - v_0;
[X_o,Y_o] = meshgrid(X_o,Y_o);

% Préparation des erreurs
map_erreur_mvs = zeros(size(X_o,1),size(Y_o,2));
map_erreur_mvsm = zeros(size(X_o,1),size(Y_o,2));
size(ind_1_shrink)
size(erreur_z_mvs)
map_erreur_mvs(ind_1_shrink) = erreur_z_mvs;
size(X_o)
size(ind_1_shrink)
size(erreur_z_mvsm)
map_erreur_mvsm(ind_1_shrink) = erreur_z_mvsm;
min_c_map = min([min(erreur_z_mvs) min(erreur_z_mvsm)]);
max_c_map = max([max(erreur_z_mvs) max(erreur_z_mvsm)]);

map_erreur_angles_GT = zeros(size(X_o,1),size(X_o,2));
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
complement_titre = complement_titre + ", profondeurs VT";


%% Analyse des résultats
normales_fronto = zeros(size(normales_mvsm));
normales_fronto(3,:) = -1;
angles_mvs = angle_normale(normales_fronto,normales_mvs);
angles_mvsm = angle_normale(normales_fronto,normales_mvsm);
angles_GT = angle_normale(normales_fronto, normales_GT);
color_map_value_GT = zeros(size(X_o,1),size(X_o,2));
color_map_value_mvs = zeros(size(X_o,1),size(X_o,2));
color_map_value_mvsm = zeros(size(X_o,1),size(X_o,2));

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

	% Bug
	size(color_map_value_GT)
	size(ind_1_shrink)
	size(indices_GT')



	color_map_value_GT(ind_1_shrink(indices_GT')) = k;
	nombre_points_zones(k) = length(indices_GT);
	erreurs_mvs_moy(k) = transpose(mean(erreur_z_mvs(indices_GT)));
	erreurs_mvs_med(k) = transpose(median(erreur_z_mvs(indices_GT)));
	erreurs_mvsm_moy(k) = transpose(mean(erreur_z_mvsm(indices_GT)));
	erreurs_mvsm_med(k) = transpose(median(erreur_z_mvsm(indices_GT)));
	label_zones = [ label_zones , zones_angles(k+1) ];
end

%% Affichage
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

n_world = R(:,:,1)' * normales_mvsm;
n_GT = R(:,:,1)' * normales_GT;
n_GT(1,1,:)

quiver3(X(ind_1_shrink),Y(ind_1_shrink),Z(ind_1_shrink),n_world(1,:)',n_world(2,:)',n_world(3,:)');
hold on;
quiver3(X(ind_1_shrink),Y(ind_1_shrink),Z(ind_1_shrink),n_GT(1,:)',n_GT(2,:)',n_GT(3,:)','r');
axis equal
hold on
surf(X,Y,Z);
colormap gray
