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
surface = "gaussienne_1_pepper";
%surface = "gaussienne_1_bis";
nombre_vues = 2;
rayon_voisinage = 4;
ecart_type_grad = -5;
ecart_type_I = -2.5;
filtrage = 0;
nombre_profondeur_iteration = 5000;
utilisation_profondeur_GT = 0;
utilisation_normale_GT = 1;
grille_pixel = 10;
mesure = "median";
mesure = "all";
utilisation_mediane_normale = 0;

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
path = "../../result/tests/";
load(path+nom_fichier);

% Préparation de la reconstruction
load('../../data/simulateur_' + surface + '_formate.mat','nombre_lignes','nombre_colonnes','facteur_k','u_0','v_0','s','N','masque','z');
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
Z_1 = z(:,:,1);
Z_1 = Z_1(1:grille_pixel:end,1:grille_pixel:end);


figure;
quiver3(X(ind_1_shrink),Y(ind_1_shrink),-Z_1(ind_1_shrink),normales_GT(1,:)',normales_GT(2,:)',-normales_GT(3,:)');
hold on
surf(X,Y,-Z_1);
axis equal;
colormap gray;
