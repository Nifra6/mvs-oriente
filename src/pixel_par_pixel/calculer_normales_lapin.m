%% Clear
clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

%% Imports
addpath(genpath('../toolbox/'));
addpath(genpath('../../../Développement/Ortho/Toolbox/'));

%% Données
load ../../data/data_bunny_ortho.mat;
load ../../data/normales_veritables.mat;
% Choix des images
indice_premiere_image = 1;
indice_deuxieme_image = 2;
% Les profondeurs
Z_1 = z(:,:,indice_premiere_image);
% Les images
I_1 = Im(:,:,indice_premiere_image);
I_2 = Im(:,:,indice_deuxieme_image);
% Les masques des images
masque_1 = mask(:,:,indice_premiere_image);
masque_2 = mask(:,:,indice_deuxieme_image);
% La pose
R_2 = R(:,:,indice_deuxieme_image) * R(:,:,indice_premiere_image)';
t_2 = t(:,indice_deuxieme_image) - R_2 * t(:,indice_premiere_image);
% Le gradient de l'image 2
[dx_I_1, dy_I_1] = gradient(I_1);
[dx_I_2, dy_I_2] = gradient(I_2);
% Les normales
N_1 = n_true;
% Caractéristiques de la caméra
u_0			= size(I_1,1)/2;
v_0			= size(I_1,2)/2;
pixelSize	= 1.5/size(I_1,2);
% Coordonnées
u_1_coord	= pixelSize*(1-u_0):pixelSize:pixelSize*(size(I_1,2)-u_0);
v_1_coord	= pixelSize*(1-v_0):pixelSize:pixelSize*(size(I_1,1)-v_0);
[Y,X] 		= meshgrid(v_1_coord,u_1_coord);

%% Paramètres
affichage_log = 0;	% Affichage d'informations diverses

%% Variables utiles
[i_1_liste, j_1_liste]	= find(masque_1);
nb_pixels_utilises		= size(i_1_liste,1);
liste_p_estimes			= zeros(size(I_1));
liste_q_estimes			= zeros(size(I_1));
liste_normales			= zeros(size(I_1,1), size(I_1,2), 3);
liste_normales(:,:,3)	= -1;

%% Algorithme

fprintf("\n");
tic;
% Sélection d'un pixel
for indice_pixel = 1:nb_pixels_utilises
	i_1 		= i_1_liste(indice_pixel);
	j_1 		= j_1_liste(indice_pixel);
	grad_I_1	= [dx_I_1(i_1,j_1); dy_I_1(i_1,j_1)];

	% Affichage de la progression des calculs
	if affichage_log & mod(indice_pixel, 100) == 0
		fprintf('\r');
		fprintf("Progression : %d / %d",indice_pixel,nb_pixels_utilises);
	end

	% Sélection de la profondeur
	z = Z_1(i_1,j_1);

	% Changements de repère
	P_1	= [pixelSize*(j_1 - u_0); pixelSize*(i_1 - v_0); z];
	P_2 = R_2 * P_1 + t_2;
	i_2 = round(P_2(2)/pixelSize + v_0);
	j_2 = round(P_2(1)/pixelSize + u_0);

	% Vérification si pixel hors image
	condition_image = i_2 > 0 & i_2 <= size(masque_2,1) & j_2 > 0 & j_2 <= size(masque_2,2);

	% Si le point reprojeté tombe sur le masque de la deuxième image
	if condition_image & masque_2(i_2,j_2)

		grad_I_2 		= [dx_I_2(i_2,j_2); dy_I_2(i_2,j_2)];
		numerateur_pq 	= grad_I_1 - R_2(1:2,1:2)' * grad_I_2;
		denominateur_pq = R_2(1:2,3)' * grad_I_2;

		% Si pas de division par 0, on continue
		if abs(denominateur_pq) > 0

			% Estimation de la pente
			p_estime = numerateur_pq(1) / denominateur_pq;
			liste_p_estimes(i_1,j_1) = p_estime;
			q_estime = numerateur_pq(2) / denominateur_pq;
			liste_q_estimes(i_1,j_1) = q_estime;

			% Calcul de la normale
			normale = (1 / sqrt(p_estime^2 + q_estime^2 + 1)) * [p_estime ; q_estime ; -1];
			liste_normales(i_1,j_1,:) = normale;
		end
	end
end
toc;

%% Résultats
displayNormalMap(N_1,'Carte des normales théoriques')
displayNormalMap(liste_normales,'Carte des normales estimees aux profondeurs réelles')

%% Fonctions annexes

function [i_limite, j_limite] = limites_voisinage(i_voisinage, j_voisinage)
	i_limite = i_voisinage(1,:)';
	i_limite = [i_limite ; i_voisinage(2:end-1,end)];
	i_limite = [i_limite ; wrev(i_voisinage(end,:))'];
	i_limite = [i_limite ; wrev(i_voisinage(2:end-1,1))];
	j_limite = j_voisinage(1,:)';
	j_limite = [j_limite ; j_voisinage(2:end-1,end)];
	j_limite = [j_limite ; wrev(j_voisinage(end,:))'];
	j_limite = [j_limite ; wrev(j_voisinage(2:end-1,1))];
end
