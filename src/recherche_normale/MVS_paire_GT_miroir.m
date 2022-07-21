%% Clear
clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

%% Données
load donnees_calotte_tr;

%% Paramètres
valeurs_z = 0:1:120;	% Les valeurs de profondeurs utilisées
lambda = 0;				% Paramètre pour l'évaluation du score de la reprojection
offset = 0;

%% Variables utiles
[i_1_liste,j_1_liste] = find(masque_1);
nb_pixels_utilises = size(i_1_liste,1);
nb_profondeurs_testees = size(valeurs_z,2);
scores = 10 * ones(nb_pixels_utilises,1);
liste_p_estimes = zeros(nb_pixels_utilises,1);
liste_q_estimes = zeros(nb_pixels_utilises,1);
angles_trouves = zeros(taille,taille);
angles_trouves_fronto = zeros(taille,taille);
angles_GT_fronto = zeros(taille,taille);
%[dy_I_1,dx_I_1] = gradient(I_1);
%[dy_I_2,dx_I_2] = gradient(I_2);
[dx_I_1_tr,dy_I_1_tr] = gradient(I_1_tr);
[dx_I_2_tr,dy_I_2_tr] = gradient(I_2_tr);
%ens_condition_image_2 = zeros(taille,taille);
ens_condition_image_2_tr = zeros(taille,taille);

f1 = figure;
imshow(I_1_tr);
title("Point sélectionné sur l'image 1");
axis on;
f2 = figure;
imshow(I_2_tr);
title("Point sélectionné reprojeté sur l'image 2 transposée");
axis on;

%% Algorithme
% Sélection d'un pixel
for indice_pixel = 1:size(i_1_liste,1)
	i_1 = i_1_liste(indice_pixel);
	j_1 = j_1_liste(indice_pixel);
	%grad_I_1 = [dx_I_1(i_1,j_1); dy_I_1(i_1,j_1)];
	grad_I_1_tr = [dx_I_1_tr(i_1,j_1); dy_I_1_tr(i_1,j_1)];

	%z_tilde = Z_1(i_1,j_1);
	z_tilde_tr = Z_1_tr(i_1,j_1);
	% Changements de repère
	%P_1_tilde = [X(i_1,j_1); Y(i_1,j_1); z_tilde];
	P_1_tilde_tr = [X_tr(i_1,j_1); Y_tr(i_1,j_1); z_tilde_tr];
	P_1_tilde_tr = [j_1 - C_x_tr ; i_1 - C_y_tr ; z_tilde_tr];
	%P_2_tilde = R_2' * P_1_tilde;
	P_2_tilde_tr = R_2_tr * P_1_tilde_tr + t_2_tr;
	%i_2_tilde = P_2_tilde(1) + C_x;
	%j_2_tilde = P_2_tilde(2) + C_y;
	j_2_tilde_tr = P_2_tilde_tr(1) + C_x_tr;
	i_2_tilde_tr = P_2_tilde_tr(2) + C_y_tr;

	if (0 & mod(i_1,10) == 1 & mod(j_1,10) == 1)
		[i_1 , j_1 , i_2_tilde_tr , j_2_tilde_tr]
	end

	if (i_1 == taille/2+delta_y & j_1 == taille/2+delta_x | i_1 == taille/2-10+delta_y & j_1 == taille/2+delta_x | i_1 == taille/2+10+delta_y & j_1 == taille/2+delta_x | i_1 == taille/2+delta_y & j_1 == taille/2-10+delta_x | i_1 == taille/2+delta_y & j_1 == taille/2+10+delta_x | i_1 == taille/2+delta_y & j_1 == taille/2-20+delta_x | i_1 == taille/2+delta_y & j_1 == taille/2-30+delta_x | i_1 == taille/2+delta_y & j_1 == taille/2+20+delta_x | i_1 == taille/2+delta_y & j_1 == taille/2+30+delta_x )
		figure(f1);
		hold on;
		plot(j_1,i_1,'r+');
		hold off;

		%figure;
		%imshow(I_2);
		%hold on;
		%plot(j_2_tilde,i_2_tilde,'r+');
		%title("Point sélectionné reprojeté sur l'image 2");
		%hold off;

		figure(f2);
		hold on;
		plot(j_2_tilde_tr,i_2_tilde_tr,'r+');
		hold off;

		[i_1,j_1]
		%P_1_tilde
		P_1_tilde_tr
		%[i_2_tilde,j_2_tilde]
		[i_2_tilde_tr,j_2_tilde_tr]
		%P_2_tilde
		P_2_tilde_tr

	end

	% Évalution de la position du point trouvé dans les dimensions des images
	%condition_image_2 = i_2_tilde > 0.5 & i_2_tilde <= size(masque_2,1) & j_2_tilde > 0.5 & j_2_tilde <= size(masque_2,2);
	condition_image_2_tr = i_2_tilde_tr > 0.5 & i_2_tilde_tr <= size(masque_2_tr,1) & j_2_tilde_tr > 0.5 & j_2_tilde_tr <= size(masque_2_tr,2);
	ens_condition_image_2_tr(i_1,j_1) = condition_image_2_tr;

	% Si le point reprojeté tombe sur le masque de la deuxième image
	%if condition_image_2 & masque_2(round(i_2_tilde),round(j_2_tilde))
	if condition_image_2_tr & masque_2_tr(round(i_2_tilde_tr),round(j_2_tilde_tr))

		%grad_I_2_tilde = [interp2(dx_I_2,j_2_tilde,i_2_tilde); interp2(dy_I_2,j_2_tilde,i_2_tilde)];
		grad_I_2_tilde_tr = [interp2(dx_I_2_tr,j_2_tilde_tr,i_2_tilde_tr); interp2(dy_I_2_tr,j_2_tilde_tr,i_2_tilde_tr)];

		%second_membre = grad_I_1 - R_2(1:2,1:2) * grad_I_2_tilde;
		second_membre_tr = grad_I_1_tr - R_2_tr(1:2,1:2)' * grad_I_2_tilde_tr;
		%coefficient_p_q = R_2(1,3) * grad_I_2_tilde(1) + R_2(2,3) * grad_I_2_tilde(2);
		coefficient_p_q_tr = R_2_tr(1,3) * grad_I_2_tilde_tr(1) + R_2_tr(2,3) * grad_I_2_tilde_tr(2);

		% Si pas de division par 0, on calcule le score associé
		%if abs(coefficient_p_q) > 0
		if abs(coefficient_p_q_tr) > 0
			%p_estime = second_membre(1) / coefficient_p_q;
			%q_estime = second_membre(2) / coefficient_p_q;
			p_estime_tr = second_membre_tr(1) / coefficient_p_q_tr;
			q_estime_tr = second_membre_tr(2) / coefficient_p_q_tr;


			%normale_calculee = [p_estime ; q_estime ; 1];
			normale_calculee_tr = [p_estime_tr ; q_estime_tr ; -1];
			%normale_calculee = normale_calculee / vecnorm(normale_calculee);
			normale_calculee_tr = normale_calculee_tr / vecnorm(normale_calculee_tr);
			%normale_GT = reshape(N_1(i_1,j_1,:),3,1);
			normale_GT_tr = reshape(N_1_tr(i_1,j_1,:),3,1);
			%angles_trouves(i_1,j_1) = angle_normale(normale_GT,normale_calculee);
			angles_trouves(i_1,j_1) = angle_normale(normale_GT_tr,normale_calculee_tr);
			angles_trouves_fronto(i_1,j_1) = angle_normale([0 ; 0 ; -1],normale_calculee_tr);
			%angles_GT_fronto(i_1,j_1) = angle_normale([0 ; 0 ; -1],normale_GT_tr);

			%liste_p_estimes(indice_pixel) = p_estime;
			%liste_q_estimes(indice_pixel) = q_estime;

			%residu_1 = I_1(i_1,j_1) - 1 / sqrt(p_estime^2 + q_estime^2 + 1);
			%residu_2 = I_1(i_1,j_1) - interp2(I_2,i_2_tilde,j_2_tilde);

			%scores(indice_pixel) = residu_2^2;
		end
	end
end


%% Résultats
% Affichage des normales
figure;
imagesc(angles_trouves);
colorbar;
title("Différences d'angles entre les normales VT et les normales estimées")
axis equal;

figure;
imagesc(angles_trouves_fronto);
colorbar;
title("Différences d'angles entre les normales fronto et les normales estimées")
axis equal;


for i = 1:taille
	for j = 1:taille
		normale_GT_tr = reshape(N_1_tr(i,j,:),3,1);
		angles_GT_fronto(i,j) = angle_normale([0 ; 0 ; -1],normale_GT_tr);
	end
end

figure;
imagesc(angles_GT_fronto);
colorbar;
title("Différences d'angles entre les normales fronto et les normales VT")
axis equal;


figure;
s = surf(X_tr,Y_tr,Z_1_tr);
s.EdgeColor = 'none';
colormap gray;
hold on
axis equal;
quiver3(X_tr,Y_tr,Z_1_tr,N_1_tr(:,:,1),N_1_tr(:,:,2),N_1_tr(:,:,3));

format long;
R_2_tr
t_2_tr
format short;
