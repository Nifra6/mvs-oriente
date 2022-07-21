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
% Les gradients des images
%[dx_I_1, dy_I_1] = gradient(I_1);
[G,imask] = make_gradient(masque_1);
Dx_main = G(1:2:end-1,:);
Dy_main = G(2:2:end,:);
clear G;
dx_I_1_calc = Dx_main * I_1(imask);
dy_I_1_calc = Dy_main * I_1(imask);
dx_I_1 = zeros(size(I_1));
dy_I_1 = zeros(size(I_1));
dx_I_1(imask) = dx_I_1_calc;
dy_I_1(imask) = dy_I_1_calc;

%[dx_I_2, dy_I_2] = gradient(I_2);
[G,imask] = make_gradient(masque_2);
Dx_main = G(1:2:end-1,:);
Dy_main = G(2:2:end,:);
clear G;
dx_I_2_calc = Dx_main * I_2(imask);
dy_I_2_calc = Dy_main * I_2(imask);
dx_I_2 = zeros(size(I_2));
dy_I_2 = zeros(size(I_2));
dx_I_2(imask) = dx_I_2_calc;
dy_I_2(imask) = dy_I_2_calc;


% Caractéristiques de la caméra
u_0 = size(I_1,1)/2;
v_0 = size(I_1,2)/2;
pixelSize = 1.5/size(I_1,2);
% Coordonnées
u_1_coord = pixelSize*(1-u_0):pixelSize:pixelSize*(size(I_1,2)-u_0);
v_1_coord = pixelSize*(1-v_0):pixelSize:pixelSize*(size(I_1,1)-v_0);
[Y,X] = meshgrid(v_1_coord,u_1_coord);

%% Paramètres
valeurs_z 		= 1:0.04:2.6;	% Les valeurs de profondeurs utilisées
rayon_voisinage	= 3;			% Voisinage à prendre en compte
affichage_log	= 0;			% Affichage d'informations diverses
interpolation 	= 'nearest';	% Type d'interpolation

%% Variables utiles
[i_1_liste, j_1_liste] = find(masque_1);
nb_pixels_utilises = size(i_1_liste,1);
nb_profondeurs_testees = size(valeurs_z,2);
scores = 10 * ones(nb_pixels_utilises, nb_profondeurs_testees);
liste_p_estimes = zeros(nb_pixels_utilises, nb_profondeurs_testees);
liste_q_estimes = zeros(nb_pixels_utilises, nb_profondeurs_testees);
liste_normales = zeros(size(I_1,1), size(I_1,2), 3);

%% Algorithme

fprintf("\n");
tic;
% Sélection d'un pixel
for indice_pixel = 1:nb_pixels_utilises
	i_1 		= i_1_liste(indice_pixel);
	j_1 		= j_1_liste(indice_pixel);
	grad_I_1	= [dx_I_1(i_1,j_1); dy_I_1(i_1,j_1)];

	% Affichage de la progression des calculs
	if mod(indice_pixel, 100) == 0
		fprintf('\r');
		fprintf("Progression : %d / %d",indice_pixel,nb_pixels_utilises);
	end

	% Sélection de la profondeur
	for indice_z = 1:nb_profondeurs_testees
		z = valeurs_z(indice_z);

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
			if abs(denominateur_pq) > 0.001

				% Estimation de la pente
				p_estime = numerateur_pq(1) / denominateur_pq;
				liste_p_estimes(indice_pixel,indice_z) = p_estime;
				q_estime = numerateur_pq(2) / denominateur_pq;
				liste_q_estimes(indice_pixel,indice_z) = q_estime;

				% Calcul du plan au pixel considéré
				normale = (1 / sqrt(p_estime^2 + q_estime^2 + 1)) * [p_estime ; q_estime ; -1];
				if (affichage_log)
					disp("===== Comparaison des normales")
					normale_theorique = reshape(N_1(i_1,j_1,:),3,1);
					normale_theorique - normale
					(180/pi) * atan2(norm(cross(normale_theorique, normale)), dot(normale_theorique, normale))
					(180/pi) * acos(dot(normale_theorique, normale) / (norm(normale_theorique)*norm(normale)))
				end
				d_equation_plan = -P_1' * normale;

				% Calcul de la transformation géométrique
				u_1_decales = pixelSize*(j_1-u_0-rayon_voisinage):pixelSize:pixelSize*(j_1-u_0+rayon_voisinage);
				v_1_decales = pixelSize*(i_1-v_0-rayon_voisinage):pixelSize:pixelSize*(i_1-v_0+rayon_voisinage);
				[v_1_decales, u_1_decales] = meshgrid(v_1_decales,u_1_decales);
				z_1_decales = -(d_equation_plan + normale(1) * u_1_decales(:) + normale(2) * v_1_decales(:)) / normale(3);
				if (affichage_log)
					disp("===== Profondeurs z")
					z
					reshape(z_1_decales, 2*rayon_voisinage+1, 2*rayon_voisinage+1)
				end

				% Reprojection du voisinage
				P_1_voisinage = [u_1_decales(:)' ; v_1_decales(:)' ; z_1_decales'];
				P_2_voisinage = R_2 * P_1_voisinage + t_2;
				i_2_voisinage = round(P_2_voisinage(2,:)/pixelSize + v_0);
				j_2_voisinage = round(P_2_voisinage(1,:)/pixelSize + u_0);
				if (affichage_log)
					disp("===== Les i du voisinage")
					i_2
					i_2_voisinage_re = reshape(i_2_voisinage, 2*rayon_voisinage+1, 2*rayon_voisinage+1)
					disp("===== Les j du voisinage")
					j_2
					j_2_voisinage_re = reshape(j_2_voisinage, 2*rayon_voisinage+1, 2*rayon_voisinage+1)
					disp("===== Le contour du voisinage")
					i_2_voisinage_re
					[i_2_limites, j_2_limites] = limites_voisinage(i_2_voisinage_re+u_0,j_2_voisinage_re+v_0);
					i_2_limites
				end

				% Récupération des niveaux de gris dans l'image 2 du voisinage	
				I_2_voisinage = reshape(interp2(I_2, j_2_voisinage(:), i_2_voisinage(:),interpolation),2*rayon_voisinage+1,2*rayon_voisinage+1)';
				scores(indice_pixel, indice_z) = sum((I_1(i_1-rayon_voisinage:i_1+rayon_voisinage,j_1-rayon_voisinage:j_1+rayon_voisinage) - I_2_voisinage).^2,'all');

				if (affichage_log)
					disp("===== Différences entre les images")
					diff = I_1(i_1-rayon_voisinage:i_1+rayon_voisinage,j_1-rayon_voisinage:j_1+rayon_voisinage) - I_2_voisinage
					sum(diff,"all")
				end
			end
		end
	end
end
toc;

%% Résultats
% Sélection des profondeurs avec le score minimal
scores_corriges = (scores ~= 0) .* scores + (scores == 0) .* ones(size(scores));
[~, indices_min] = min(scores_corriges, [], 2);
z_in = transpose(valeurs_z(indices_min));
z = zeros(size(I_1,1), size(I_1,2));
z(find(masque_1)) = z_in;

% Affichage
%figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
%plot3(X,Y,-z,'k.');
%xlabel('$x$','Interpreter','Latex','FontSize',30);
%ylabel('$y$','Interpreter','Latex','FontSize',30);
%zlabel('$z$','Interpreter','Latex','FontSize',30);
%axis equal;
%rotate3d;

% Affichage du résultat :
XYZ_test = shapeFromDmOrtho(z/pixelSize, masque_1);

figure('Name', 'z from estimated (p,q)')
surfl(XYZ_test(:,:,1),XYZ_test(:,:,2),-XYZ_test(:,:,3),[0 90])
shading flat
colormap gray
axis ij
axis tight
axis off

% Intégration théorique
z_true = integration_SCS(q_true,p_true);
XYZ_true = shapeFromDmOrtho(z_true/pixelSize, masque_1);

figure('Name', 'z integrated from (p_true,q_true)')
surfl(XYZ_true(:,:,1),XYZ_true(:,:,2),-XYZ_true(:,:,3),[0 90])
shading flat
colormap gray
axis ij
axis tight
axis off

% Intégration théorique
[G,imask] = make_gradient(masque_1);
Dx_main = G(1:2:end-1,:);
Dy_main = G(2:2:end,:);
clear G;
p_calc = Dx_main * Z_1(imask) / pixelSize;
q_calc = Dy_main * Z_1(imask) / pixelSize;
p_true_2 = zeros(size(I_1));
q_true_2 = zeros(size(I_1));
p_true_2(imask) = p_calc;
q_true_2(imask) = q_calc;

[dx_Z_1, dy_Z_1] = gradient(Z_1/pixelSize);
%dx_Z_1(find(dx_Z_1 >= 100)) = 0;
%dy_Z_1(find(dy_Z_1 >= 100)) = 0;
p_true_2 = zeros(size(I_1,1),size(I_1,2));
q_true_2 = zeros(size(p_true_2));
p_true_2(find(masque_1)) = dx_Z_1(find(masque_1));
q_true_2(find(masque_1)) = dy_Z_1(find(masque_1));

z_true_2 = integration_SCS(q_true_2,p_true_2);
XYZ_true_2 = shapeFromDmOrtho(z_true_2/pixelSize, masque_1);

figure('Name', 'z integrated from (p_true_2,q_true_2)')
surfl(XYZ_true_2(:,:,1),XYZ_true_2(:,:,2),-XYZ_true_2(:,:,3),[0 90])
shading flat
colormap gray
axis ij
axis tight
axis off

% Intégration depuis les p et q
p = zeros(size(I_1,1), size(I_1,2));
q = zeros(size(p));
indices_lin = sub2ind(size(liste_p_estimes), (1:nb_pixels_utilises)', indices_min);
p(find(masque_1)) = liste_p_estimes(indices_lin);
q(find(masque_1)) = liste_q_estimes(indices_lin);
z_int = integration_SCS(q,p);

XYZ_int = shapeFromDmOrtho(z_int/pixelSize, masque_1);

figure('Name', 'z integrated from estimated (p,q)')
surfl(XYZ_int(:,:,1),XYZ_int(:,:,2),-XYZ_int(:,:,3),[0 90])
shading flat
colormap gray
axis ij
axis tight
axis off

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
