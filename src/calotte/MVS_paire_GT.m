%% Clear
clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

%% Données
load donnees_calotte;

%% Paramètres
valeurs_z = 0:1:120;	% Les valeurs de profondeurs utilisées
lambda = 0;				% Paramètre pour l'évaluation du score de la reprojection

%% Variables utiles
[i_1_liste,j_1_liste] = find(masque_1);
nb_pixels_utilises = size(i_1_liste,1);
nb_profondeurs_testees = size(valeurs_z,2);
scores = 10 * ones(nb_pixels_utilises);
liste_p_estimes = zeros(nb_pixels_utilises);
liste_q_estimes = zeros(nb_pixels_utilises);
angles_trouves = zeros(taille,taille);
[dy_I_1,dx_I_1] = gradient(I_1);
[dy_I_2,dx_I_2] = gradient(I_2);

%% Algorithme
% Sélection d'un pixel
for indice_pixel = 1:size(i_1_liste,1)
	i_1 = i_1_liste(indice_pixel);
	j_1 = j_1_liste(indice_pixel);
	grad_I_1 = [dx_I_1(i_1,j_1); dy_I_1(i_1,j_1)];

	z_tilde = Z_1(i_1,j_1);
	% Changements de repère
	P_1_tilde = [X(i_1,j_1); Y(i_1,j_1); z_tilde];
	P_2_tilde = R_2' * P_1_tilde;
	i_2_tilde = P_2_tilde(1) + C_x;
	j_2_tilde = P_2_tilde(2) + C_y;

	if (i_1 == 128 & j_1 == 128)
		figure;
		imshow(I_1);
		hold on;
		plot(j_1,i_1,'r+');
		hold off;

		figure;
		imshow(I_2);
		hold on;
		plot(j_2_tilde,i_2_tilde,'r+');
		hold off;

		[i_1,j_1]
		[i_2_tilde,j_2_tilde]
	end

	% Évalution de la position du point trouvé dans les dimensions des images
	condition_image_2 = i_2_tilde > 0.5 & i_2_tilde <= size(masque_2,1) & j_2_tilde > 0.5 & j_2_tilde <= size(masque_2,2);

	% Si le point reprojeté tombe sur le masque de la deuxième image
	if condition_image_2 & masque_2(round(i_2_tilde),round(j_2_tilde))

		grad_I_2_tilde = [interp2(dx_I_2,j_2_tilde,i_2_tilde); interp2(dy_I_2,j_2_tilde,i_2_tilde)];
		%grad_I_2_tilde = [dx_I_2(round(i_2_tilde),round(j_2_tilde)); dy_I_2(round(i_2_tilde),round(j_2_tilde))];

		second_membre = grad_I_1 - R_2(1:2,1:2) * grad_I_2_tilde;
		coefficient_p_q = R_2(1,3) * grad_I_2_tilde(1) + R_2(2,3) * grad_I_2_tilde(2);

		% Si pas de division par 0, on calcule le score associé
		if abs(coefficient_p_q) > 0
			p_estime = second_membre(1) / coefficient_p_q;
			q_estime = second_membre(2) / coefficient_p_q;


			normale_calculee = [p_estime ; q_estime ; 1];
			normale_calculee = normale_calculee / vecnorm(normale_calculee);
			normale_GT = reshape(N_1(i_1,j_1,:),3,1);
			angles_trouves(i_1,j_1) = angle_normale(normale_GT,normale_calculee);

			liste_p_estimes(indice_pixel) = p_estime;
			liste_q_estimes(indice_pixel) = q_estime;

			residu_1 = I_1(i_1,j_1) - 1 / sqrt(p_estime^2 + q_estime^2 + 1);
			residu_2 = I_1(i_1,j_1) - interp2(I_2,i_2_tilde,j_2_tilde);

			scores(indice_pixel) = residu_2^2;
		end
	end
end


%% Résultats
% Affichage des normales
figure;
imagesc(angles_trouves);
colorbar;
axis equal;

figure
imshow(I_1);
figure
imshow(I_2);
