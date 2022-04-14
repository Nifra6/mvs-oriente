%% Clear
clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

%% Données
load donnees_calotte;
R2 = inv(R_2);
R3 = inv(R_3);


%% Paramètres
valeurs_z = 0:1:120;	% Les valeurs de profondeurs utilisées
lambda = 0;				% Paramètre pour l'évaluation du score de la reprojection
deuxieme_image = 1;		% Si oui, utilisation d'une deuxième rotation

%% Variables utiles
[i_1_liste,j_1_liste] = find(masque_1);
nb_pixels_utilises = size(i_1_liste,1);
nb_profondeurs_testees = size(valeurs_z,2);
scores = 10 * ones(nb_pixels_utilises, nb_profondeurs_testees);
liste_p_estimes = zeros(nb_pixels_utilises, nb_profondeurs_testees);
liste_q_estimes = zeros(nb_pixels_utilises, nb_profondeurs_testees);
liste_normales = zeros(taille,taille,3);
Omega = zeros(taille,taille);

%% Pixel de vérification
figure;
imshow(I_1);
P = drawpoint;
pos = P.Position;
i_1_test = round(pos(2));
j_1_test = round(pos(1));
close;


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
		P_1 = [X(i_1,j_1); Y(i_1,j_1); z_tilde];
		P_2 = R2 * P_1;
		P_3 = R3 * P_1;
		i_2 = round(P_2(1) + C_x);
		j_2 = round(P_2(2) + C_y);
		i_3 = round(P_3(1) + C_x);
		j_3 = round(P_3(2) + C_y);

		% Évalution de la position du point trouvé dans les dimensions des images
		condition_image_2 = i_2 > 0 & i_2 <= size(masque_2,1) & j_2 > 0 & j_2 <= size(masque_2,2);
		condition_image_3 = i_3 > 0 & i_3 <= size(masque_3,1) & j_3 > 0 & j_3 <= size(masque_3,2);

		if deuxieme_image
			% Si le point reprojeté tombe sur les deux masques des deux images
			if condition_image_2 & masque_2(i_2,j_2) & condition_image_3 & masque_3(i_3,j_3)

				grad_I_2_tilde = [dx_I_2(i_2,j_2); dy_I_2(i_2,j_2)];
				grad_I_3_tilde = [dx_I_3(i_3,j_3); dy_I_3(i_3,j_3)];

				denominateur_pq_2 = R2(1:2,3)' * grad_I_2_tilde;
				numerateur_pq_2 = grad_I_1 - R2(1:2,1:2)' * grad_I_2_tilde;
				denominateur_pq_3 = R3(1:2,3)' * grad_I_3_tilde;
				numerateur_pq_3 = grad_I_1 - R3(1:2,1:2)' * grad_I_3_tilde;


				% Si pas de division par 0, on calcule le score associé
				if abs(denominateur_pq_2) > 0 & abs(denominateur_pq_3) > 0
					Omega(i_1,j_1) = 1;
					p_estime_2 = -numerateur_pq_2(1) / denominateur_pq_2;
					q_estime_2 = -numerateur_pq_2(2) / denominateur_pq_2;
					p_estime_3 = -numerateur_pq_3(1) / denominateur_pq_3;
					q_estime_3 = -numerateur_pq_3(2) / denominateur_pq_3;

					normale_2 = [p_estime_2 ; q_estime_2 ; 1];
					normale_2 = normale_2 / norm(normale_2);
					normale_3 = [p_estime_3 ; q_estime_3 ; 1];
					normale_3 = normale_3 / norm(normale_3);
					if (i_1 == i_1_test & j_1 == j_1_test & round(z_tilde) == round(Z_1(i_1,j_1)))
						normale = reshape(N_1(i_1,j_1,:),3,1)
						disp("Normale avec image 2");
						normale(1) + normale_2(1)
						normale(2) + normale_2(2)
						normale(3) - normale_2(3)
						disp("Normale avec image 3");
						normale - normale_3
					end
					if (z_tilde == round(Z_1(i_1,j_1)))
						liste_normales(i_1,j_1,:) = normale_2;
					end

					%liste_p_estimes(indice_pixel, indice_z_tilde) = p_estime;
					%liste_q_estimes(indice_pixel, indice_z_tilde) = q_estime;

					residu_2 = I_1(i_1,j_1) - 1 / sqrt(p_estime_2^2 + q_estime_2^2 + 1);
					residu_2_2 = I_1(i_1,j_1) - I_2(i_2,j_2);
					residu_3 = I_1(i_1,j_1) - 1 / sqrt(p_estime_3^2 + q_estime_3^2 + 1);
					residu_3_2 = I_1(i_1,j_1) - I_3(i_3,j_3);

					scores(indice_pixel, indice_z_tilde) = residu_2^2 + residu_3^2 + lambda * (residu_2_2^2 + residu_3_2^2);
				end
			end
		else
			% Si le point reprojeté tombe sur le masque de la deuxième image
			if condition_image_2 & masque_2(i_2,j_2)

				grad_I_2_tilde = [dx_I_2(i_2,j_2); dy_I_2(i_2,j_2)];

				denominateur_pq = R_2(3,1) * grad_I_2_tilde(1) + R_2(3,2) * grad_I_2_tilde(2);
				numerateur_pq = grad_I_1 - R_2(1:2,1:2) * grad_I_2_tilde;

				% Si pas de division par 0, on calcule le score associé
				if abs(denominateur_pq) > 0
					Omega(i_1,j_1) = 1;
					p_estime = numerateur_pq(1) / denominateur_pq;
					q_estime = numerateur_pq(2) / denominateur_pq;

					liste_p_estimes(indice_pixel, indice_z_tilde) = p_estime;
					liste_q_estimes(indice_pixel, indice_z_tilde) = q_estime;

					residu_1 = I_1(i_1,j_1) - 1 / sqrt(p_estime^2 + q_estime^2 + 1);
					residu_2 = I_1(i_1,j_1) - I_2(i_2,j_2);

					scores(indice_pixel, indice_z_tilde) = residu_1^2 + lambda * residu_2^2;
				end
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


% Estimation de z avec p et q
p_estimes = liste_p_estimes(indices_min);
q_estimes = liste_p_estimes(indices_min);
p = zeros(taille,taille);
q = zeros(taille,taille);
p(find(masque_1)) = p_estimes;
q(find(masque_1)) = q_estimes;
z_estim = integration_SCS(p,q);


% Affichage du relief reconstruit sur le domaine Omega :
%figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
%z_estim(Omega==0) = nan;
%plot3(X,Y,z_estim,'k.');
%xlabel('$x$','Interpreter','Latex','FontSize',30);
%ylabel('$y$','Interpreter','Latex','FontSize',30);
%zlabel('$z$','Interpreter','Latex','FontSize',30);
%axis equal;
%rotate3d;

figure;
imshow(N_1(:,:,1))

figure;
imshow(liste_normales(:,:,1))
