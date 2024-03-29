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
scores = 10 * ones(nb_pixels_utilises, nb_profondeurs_testees);
liste_p_estimes = zeros(nb_pixels_utilises, nb_profondeurs_testees);
liste_q_estimes = zeros(nb_pixels_utilises, nb_profondeurs_testees);
angles_trouves = zeros(taille,taille);

%% Algorithme
% Sélection d'un pixel
for indice_pixel = 1:size(i_1_liste,1)
	i_1 = i_1_liste(indice_pixel);
	j_1 = j_1_liste(indice_pixel);
	grad_I_1 = [dx_I_1(i_1,j_1); dy_I_1(i_1,j_1)];

	% Sélection d'une profondeur
	for indice_z_tilde = 1:nb_profondeurs_testees
		z_tilde = valeurs_z(indice_z_tilde);
		% Changements de repère
		P_1_tilde = [X(i_1,j_1); Y(i_1,j_1); z_tilde];
		P_2_tilde = R_2' * P_1_tilde;
		i_2_tilde = P_2_tilde(1) + C_x;
		j_2_tilde = P_2_tilde(2) + C_y;

		% Évalution de la position du point trouvé dans les dimensions des images
		condition_image_2 = i_2_tilde > 0 & i_2_tilde <= size(masque_2,1) & j_2_tilde > 0 & j_2_tilde <= size(masque_2,2);

		% Si le point reprojeté tombe sur le masque de la deuxième image
		if condition_image_2 & masque_2(round(i_2_tilde),round(j_2_tilde))

			grad_I_2_tilde = [interp2(dx_I_2,i_2_tilde,j_2_tilde); interp2(dy_I_2,i_2_tilde,j_2_tilde)];

			coefficient_p_q = R_2(3,1) * grad_I_2_tilde(1) + R_2(3,2) * grad_I_2_tilde(2);
			second_membre = grad_I_1 - R_2(1:2,1:2) * grad_I_2_tilde;

			% Si pas de division par 0, on calcule le score associé
			if abs(coefficient_p_q) > 0
				p_estime = second_membre(1) / coefficient_p_q;
				q_estime = second_membre(2) / coefficient_p_q;


				normale_calculee = [p_estime ; q_estime ; 1];
				normale_calculee = normale_calculee / vecnorm(normale_calculee);
				normale_GT = reshape(N_1(i_1,j_1,:),3,1);
				angles_trouves(i_1,j_1) = angle_normale(normale_GT,normale_calculee);

				liste_p_estimes(indice_pixel, indice_z_tilde) = p_estime;
				liste_q_estimes(indice_pixel, indice_z_tilde) = q_estime;

				residu_1 = I_1(i_1,j_1) - 1 / sqrt(p_estime^2 + q_estime^2 + 1);
				residu_2 = I_1(i_1,j_1) - interp2(I_2(i_2_tilde,j_2_tilde));

				scores(indice_pixel, indice_z_tilde) = residu_1^2 + lambda * residu_2^2;
				scores(indice_pixel, indice_z_tilde) = residu_2^2;
			end
		end
	end
end


%% Résultats
% Sélection des profondeurs avec le score minimal
[A, indices_min] = min(scores, [], 2);
z_in = transpose(valeurs_z(indices_min));
z = zeros(256, 256);
z(find(masque_1)) = z_in;

% Affichage
figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
plot3(X,Y,z,'k.');
%plot3(X,Y,Z_1,'k.');
xlabel('$x$','Interpreter','Latex','FontSize',30);
ylabel('$y$','Interpreter','Latex','FontSize',30);
zlabel('$z$','Interpreter','Latex','FontSize',30);
axis equal;
rotate3d;

% Affichage des normales
figure;
imagesc(angles_trouves);
colorbar;


% Estimation de z avec p et q
p_estimes = liste_p_estimes(indices_min);
q_estimes = liste_p_estimes(indices_min);
p = zeros(taille,taille);
q = zeros(taille,taille);
p(find(masque_1)) = p_estimes;
q(find(masque_1)) = q_estimes;
z_estim = integration_SCS(p,q);

