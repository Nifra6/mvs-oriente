%% Clear
clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

%% Données
load ../data/donnees_calotte;
I_1 = I(:,:,1);
I_2 = I(:,:,2);
masque_1 = masque(:,:,1);
masque_2 = masque(:,:,2);
R_2 = R(:,:,1);
t_2 = t(:,1);
[dy_I_2, dx_I_2] = gradient(I_2);

homographie_infini = K * R_2' * K_inv;

%% Paramètres
valeurs_z = 60:1:120;	% Les valeurs de profondeurs utilisées
lambda = 0;				% Paramètre pour l'évaluation du score de la reprojection
deuxieme_image = 0;		% Si oui, utilisation d'une deuxième rotation

%% Variables utiles
[i_1_liste,j_1_liste] = find(masque_1);
nb_pixels_utilises = size(i_1_liste,1);
nb_profondeurs_testees = size(valeurs_z,2);
scores = 10 * ones(nb_pixels_utilises, nb_profondeurs_testees);
liste_p_estimes = zeros(nb_pixels_utilises, nb_profondeurs_testees);
liste_q_estimes = zeros(nb_pixels_utilises, nb_profondeurs_testees);
Omega = zeros(taille,taille);

%% Algorithme
% Sélection d'un pixel
for indice_pixel = 1:size(i_1_liste,1)
	%disp(int2str(indice_pixel) + " sur " + int2str(size(i_1_liste,1)));
	i_1 = i_1_liste(indice_pixel);
	j_1 = j_1_liste(indice_pixel);
	grad_I_1 = [dx_I_1(i_1,j_1); dy_I_1(i_1,j_1)];

	% Sélection d'une profondeur
	for indice_z_tilde = 1:nb_profondeurs_testees
		z_tilde = valeurs_z(indice_z_tilde);
		% Changements de repère
		P_1_tilde = [X(i_1,j_1); Y(i_1,j_1); z_tilde];
		P_2_tilde = R_2' * P_1_tilde;
		i_2_tilde = round(P_2_tilde(1) + C_x);
		j_2_tilde = round(P_2_tilde(2) + C_y);

		% Évalution de la position du point trouvé dans les dimensions des images
		condition_image_2 = i_2_tilde > 0 & i_2_tilde <= size(masque_2,1) & j_2_tilde > 0 & j_2_tilde <= size(masque_2,2);

		% Si le point reprojeté tombe sur le masque de la deuxième image
		if condition_image_2 & masque_2(i_2_tilde,j_2_tilde)

			grad_I_2_tilde = [dx_I_2(i_2_tilde,j_2_tilde); dy_I_2(i_2_tilde,j_2_tilde)];

			coefficient_p_q = R_2(3,1) * grad_I_2_tilde(1) + R_2(3,2) * grad_I_2_tilde(2);
			second_membre = grad_I_1 - R_2(1:2,1:2) * grad_I_2_tilde;

			% Si pas de division par 0, on calcule le score associé
			if abs(coefficient_p_q) > 0
				Omega(i_1,j_1) = 1;
				p_estime = second_membre(1) / coefficient_p_q;
				q_estime = second_membre(2) / coefficient_p_q;

				liste_p_estimes(indice_pixel, indice_z_tilde) = p_estime;
				liste_q_estimes(indice_pixel, indice_z_tilde) = q_estime;


				% Calcul homographie
				normale = (1 / (p_estime^2 + q_estime^2 + 1)) * [p_estime ; q_estime ; 1];
				homographie = homographie_infini + (1/z_tilde) * t_2 * normale' * K_inv;


				if (indice_pixel == 12506 && z_tilde == round(Z_1(i_1,j_1)))
					range = 4;
					i_br = i_1-range:i_1+range;
					j_br = j_1-range:j_1+range;
					[j_liste, i_liste] = meshgrid(i_br,j_br);
					H_ij = homographie * [i_liste(:)' ; j_liste(:)' ; ones(1,size(i_liste,1)*size(i_liste,2))]
					G_ij = H_ij(1:2,:) ./ H_ij(3,:)
					subplot(2,2,1);
					imshow(I_1(i_1-range:i_1+range,j_1-range:j_1+range));
					subplot(2,2,2);
					i_bruh = reshape(interp2(I_2, G_ij(2,:), G_ij(1,:),'cubic'),size(i_liste,1),size(j_liste,1));
					i_bruh
					imshow(i_bruh);
					subplot(2,2,3);
					%i_truc = I_1(i_1-range:i_1+range,j_1-range:j_1+range) - i_bruh
					imshow(I_1);
					axis on
					hold on;
					plot(j_1,i_1, 'r+', 'MarkerSize', 30, 'LineWidth', 2)
					hold off;
					subplot(2,2,4);
					imshow(I_2)
					axis on
					hold on;
					plot(j_2_tilde,i_2_tilde, 'r+', 'MarkerSize', 30, 'LineWidth', 2)
					i_2_pouet = G_ij(1,floor(range^2 / 2)+1);
					j_2_pouet = G_ij(1,floor(range^2 / 2)+1);
					plot(j_2_pouet,i_2_pouet, 'b+', 'MarkerSize', 30, 'LineWidth', 2);
					hold off;
				end
				H_ij = homographie * [i_1 ; j_1 ; 1];
				G_ij = H_ij(1:2) / H_ij(3);


				residu_1 = I_1(i_1,j_1) - 1 / sqrt(p_estime^2 + q_estime^2 + 1);
				%residu_1 = I_1(i_1,j_1) - interp2(I_2,G_ij(2),G_ij(1),'nearest');
				residu_2 = I_1(i_1,j_1) - I_2(i_2_tilde,j_2_tilde);

				scores(indice_pixel, indice_z_tilde) = residu_1^2 + lambda * residu_2^2;
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
%figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
%plot3(X,Y,z,'k.');
%%plot3(X,Y,Z_1,'k.');
%xlabel('$x$','Interpreter','Latex','FontSize',30);
%ylabel('$y$','Interpreter','Latex','FontSize',30);
%zlabel('$z$','Interpreter','Latex','FontSize',30);
%axis equal;
%rotate3d;


% Estimation de z avec p et q
p_estimes = liste_p_estimes(indices_min);
q_estimes = liste_p_estimes(indices_min);
p = zeros(taille,taille);
q = zeros(taille,taille);
p(find(masque_1)) = p_estimes;
q(find(masque_1)) = q_estimes;
z_estim = integration_SCS(p,q);
