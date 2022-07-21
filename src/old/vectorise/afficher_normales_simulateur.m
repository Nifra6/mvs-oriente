%% Clear
clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

%% Imports
addpath(genpath('../../../Développement/Ortho/Toolbox/'));
addpath(genpath('../toolbox/'));

%% Données
load ../../data/simulateur_formate.mat;
load ../../data/normales_simulateur_estimes.mat;
% Choix des images
indice_premiere_image = 1;
% Les normales
N_1 = N(:,:,:,indice_premiere_image);
% Les images
I_1 = I(:,:,indice_premiere_image);
% Les masques des images
masque_1 = masque(:,:,indice_premiere_image);


%% Paramètres
affichage_log = 0;	% Affichage d'informations diverses

%% Variables utiles
[i_1_liste, j_1_liste]	= find(masque_1);
nb_pixels_utilises		= size(i_1_liste,1);
espacement_normales = 5;
facteur_normales 	= 2;
affichage_normales_calculees = 1;

%% Algorithme

fprintf("\n");
tic;
figure;
imshow(I_1);
hold on;
% Sélection d'un pixel
for indice_pixel = 1:nb_pixels_utilises
	i_1 		= i_1_liste(indice_pixel);
	j_1 		= j_1_liste(indice_pixel);

	if mod(i_1,espacement_normales) == 1 & mod(j_1,espacement_normales) == 1

		% ----- Affichage de la normale théorique
		normale_theorique = reshape(N_1(i_1,j_1,:),3,1);
		i_arrow_th = i_1 + facteur_normales * normale_theorique(2);
		j_arrow_th = j_1 + facteur_normales * normale_theorique(1);
		X_th = [j_1 j_arrow_th];
		Y_th = [i_1 i_arrow_th];
		plot(X_th, Y_th, 'b-', 'MarkerSize', 30, 'LineWidth', 2);
		
		% ----- Affichage de la normale estimée
		normale_estimee = n_totales(:,sub2ind([nombre_lignes, nombre_colonnes], i_1, j_1));
		i_arrow_cal = i_1 + facteur_normales * normale_estimee(2);
		j_arrow_cal = j_1 + facteur_normales * normale_estimee(1);
		X_cal = [j_1 j_arrow_cal];
		Y_cal = [i_1 i_arrow_cal];
		plot(X_cal, Y_cal, 'r-', 'MarkerSize', 30, 'LineWidth', 2);

	end

end
toc;
hold off;
