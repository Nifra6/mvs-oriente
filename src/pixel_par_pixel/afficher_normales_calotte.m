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
load ../../data/donnees_calotte.mat;
% Choix des images
indice_premiere_image = 1;
indice_deuxieme_image = 2;
% Les images
I_1 = I(:,:,indice_premiere_image);
I_2 = I(:,:,indice_deuxieme_image);
% Les masques des images
masque_1 = masque(:,:,indice_premiere_image);
masque_2 = masque(:,:,indice_deuxieme_image);
% La pose
R_1_2 = R(:,:,indice_deuxieme_image) * R(:,:,indice_premiere_image)';
t_1_2 = t(:,indice_deuxieme_image) - R_1_2 * t(:,indice_premiere_image);

% Le gradient de l'image 2
%[G,imask] = make_gradient(masque_1);
%Dx_main = G(1:2:end-1,:);
%Dy_main = G(2:2:end,:);
%clear G;
%dx_I_1_calc = Dx_main * I_1(imask);
%dy_I_1_calc = Dy_main * I_1(imask);
%dx_I_1 = zeros(size(I_1));
%dy_I_1 = zeros(size(I_1));
%dx_I_1(imask) = dx_I_1_calc;
%dy_I_1(imask) = dy_I_1_calc;
%clear imask;

%[G,imask] = make_gradient(masque_2);
%Dx_main = G(1:2:end-1,:);
%Dy_main = G(2:2:end,:);
%clear G;
%dx_I_2_calc = Dx_main * I_2(imask);
%dy_I_2_calc = Dy_main * I_2(imask);
%dx_I_2 = zeros(size(I_2));
%dy_I_2 = zeros(size(I_2));
%dx_I_2(imask) = dx_I_2_calc;
%dy_I_2(imask) = dy_I_2_calc;
%clear imask;

pixelSize	= 1.5/size(I_1,2);
[dx_I_1, dy_I_1] = gradient(I_1);
[dx_I_2, dy_I_2] = gradient(I_2);

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
	grad_I_1	= [dx_I_1(i_1,j_1); dy_I_1(i_1,j_1)];

	if mod(i_1,espacement_normales) == 1 & mod(j_1,espacement_normales) == 1

		% ----- Affichage de la normale théorique
		normale_theorique = reshape(N_1(i_1,j_1,:),3,1);
		i_arrow_th = i_1 + facteur_normales * normale_theorique(2);
		j_arrow_th = j_1 + facteur_normales * normale_theorique(1);
		X_th = [j_1 j_arrow_th];
		Y_th = [i_1 i_arrow_th];
		plot(X_th, Y_th, 'b-', 'MarkerSize', 30, 'LineWidth', 2);

		% ------ Calcul de la normale avec p et q
		if affichage_normales_calculees

			% Récupération de la profondeur
			z = Z_1(i_1,j_1);

			% Changements de repère
			u_1 = j_1 - u_0;
			v_1 = i_1 - v_0;
			P_1 = [u_1 ; v_1 ; z];
			P_2 = R_1_2 * P_1 + t_1_2;
			u_2 = round(P_2(1));
			v_2 = round(P_2(2));
			i_2 = v_2 + v_0;
			j_2 = u_2 + u_0;

			% Vérification si pixel hors image
			condition_image = i_2 > 0 & i_2 <= size(masque_2,1) & j_2 > 0 & j_2 <= size(masque_2,2);

			% Si le point reprojeté tombe sur le masque de la deuxième image
			if condition_image && masque_2(i_2,j_2)

				grad_I_2 		= [dx_I_2(i_2,j_2); dy_I_2(i_2,j_2)];
				numerateur_pq 	= grad_I_1 - R_1_2(1:2,1:2)' * grad_I_2;
				denominateur_pq = R_1_2(1:2,3)' * grad_I_2;

				% Si pas de division par 0, on continue
				if abs(denominateur_pq) > 0.001

					% Estimation de la pente
					p_estime = numerateur_pq(1) / denominateur_pq;
					q_estime = numerateur_pq(2) / denominateur_pq;

					% Calcul du plan au pixel considéré
					normale = (1 / sqrt(p_estime^2 + q_estime^2 + 1)) * [p_estime ; q_estime ; -1];
					i_arrow_cal = i_1 + facteur_normales * normale(2);
					j_arrow_cal = j_1 + facteur_normales * normale(1);
					X_cal = [j_1 j_arrow_cal];
					Y_cal = [i_1 i_arrow_cal];
					plot(X_cal, Y_cal, 'r-', 'MarkerSize', 30, 'LineWidth', 2);

				end
			end

		end

	end

end
toc;
hold off;
