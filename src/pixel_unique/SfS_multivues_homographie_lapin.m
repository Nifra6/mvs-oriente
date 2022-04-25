%% Clear
clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);
format long;

%% Imports
addpath(genpath('../../../Développement/Ortho/Toolbox/'))

%% Données
load ../../data/data_bunny_ortho;
load ../../data/normales_veritables;
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
R_1_2 = R(:,:,indice_deuxieme_image) * R(:,:,indice_premiere_image)';
t_1_2 = t(:,indice_deuxieme_image) - R_1_2 * t(:,indice_premiere_image);
% Le gradient de l'image 2
[dx_I_1, dy_I_1] = gradient(I_1);
[dx_I_2, dy_I_2] = gradient(I_2);

% Les normales
N_1 = n_true;

%figure()
%displayNormalMap(N_1,'GT normMap');

% Caractéristique de la caméra
u_0 		= size(I_1,2)/2;	% En pix
v_0 		= size(I_1,1)/2;	% En pix
pixelSize	= 1.5/size(I_1,2);	% En m.pix^{-1}

%% Paramètres
valeurs_z 		= 1:0.04:2.6;	% Les valeurs de profondeurs utilisées
rayon_voisinage	= 3;			% Voisinage à prendre en compte
affichage_log	= 1;			% Affichage d'informations diverses

%% Algorithme
while (1)
	% Sélection d'un pixel
	figure;
	title("Cliquez sur le pixel souhaité.")
	imshow(I_1);
	P			= drawpoint;
	pos 		= P.Position;
	i_1 		= round(pos(2))
	j_1 		= round(pos(1))
	%i_1 = 220;
	%j_1 = 175;
	grad_I_1	= [dx_I_1(i_1,j_1); dy_I_1(i_1,j_1)];

	% Récupération de la profondeur
	z = Z_1(i_1,j_1);

	% Changements de repère
	u_1 = j_1 - u_0;
	v_1 = i_1 - v_0;
	P_1 = [pixelSize * u_1 ; pixelSize * v_1 ; z];
	P_2 = R_1_2 * P_1 + t_1_2;
	u_2 = round(P_2(1)/pixelSize);
	v_2 = round(P_2(2)/pixelSize);
	i_2 = v_2 + v_0;
	j_2 = u_2 + u_0;

	% Vérification si pixel hors image
	condition_image = i_2 > 0 & i_2 <= size(masque_2,1) & j_2 > 0 & j_2 <= size(masque_2,2);

	% Si le point reprojeté tombe sur le masque de la deuxième image
	if condition_image & masque_2(i_2,j_2)

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
			if (affichage_log)
				disp("===== Comparaison des normales")
				normale_theorique = reshape(N_1(i_1,j_1,:),3,1);
				normale_theorique = normale_theorique / norm(normale_theorique);
				normale
				normale_theorique
				normale_theorique - normale
				(180/pi) * atan2(norm(cross(normale_theorique,normale)),dot(normale_theorique,normale))
			end
			%normale = (1/norm(normale_theorique)) * normale_theorique;
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
				[Z_1(i_1-rayon_voisinage,j_1-rayon_voisinage),Z_1(i_1-rayon_voisinage,j_1+rayon_voisinage),Z_1(i_1+rayon_voisinage,j_1+rayon_voisinage),Z_1(i_1+rayon_voisinage,j_1-rayon_voisinage)]
			end

			% Reprojection du voisinage
			P_1_voisinage = [u_1_decales(:)' ; v_1_decales(:)' ; z_1_decales'];
			P_2_voisinage = R_1_2 * P_1_voisinage + t_1_2;
			i_2_voisinage = round(P_2_voisinage(2,:)/pixelSize + v_0);
			j_2_voisinage = round(P_2_voisinage(1,:)/pixelSize + u_0);
			i_2_voisinage_re = reshape(i_2_voisinage, 2*rayon_voisinage+1, 2*rayon_voisinage+1);
			j_2_voisinage_re = reshape(j_2_voisinage, 2*rayon_voisinage+1, 2*rayon_voisinage+1);
			if (affichage_log)
				disp("===== Les i du voisinage")
				i_2
				i_2_voisinage_re
				disp("===== Les j du voisinage")
				j_2
				j_2_voisinage_re
				disp("===== Le contour du voisinage")
				i_2_voisinage_re
				[i_2_limites, j_2_limites] = limites_voisinage(i_2_voisinage_re+u_0,j_2_voisinage_re+v_0);
				i_2_limites
			end

			% Récupération des niveaux de gris dans l'image 2 du voisinage	
			I_2_voisinage = reshape(interp2(I_2, j_2_voisinage(:), i_2_voisinage(:),'cubic'),2*rayon_voisinage+1,2*rayon_voisinage+1)';
			if (affichage_log)
				disp("===== Différences entre les images")
				diff = I_1(i_1-rayon_voisinage:i_1+rayon_voisinage,j_1-rayon_voisinage:j_1+rayon_voisinage) - I_2_voisinage
				sum(diff,"all")
			end


			% Affichage des résultats
			figure;

			% Image 1
			subplot(2,2,1);
			imshow(I_1(i_1-rayon_voisinage:i_1+rayon_voisinage,j_1-rayon_voisinage:j_1+rayon_voisinage));
			title("Image 1 au voisinage")

			% Image 2
			subplot(2,2,2);
			imshow(I_2_voisinage);
			title("Image 2 au voisinage")

			% Voisinage 1
			subplot(2,2,3);
			imshow(I_1);
			axis on
			hold on;
			i_1_voisinage = i_1-rayon_voisinage:i_1+rayon_voisinage;
			j_1_voisinage = j_1-rayon_voisinage:j_1+rayon_voisinage;
			[i_1_voisinage, j_1_voisinage] = meshgrid(i_1_voisinage,j_1_voisinage);
			[i_1_limites, j_1_limites] = limites_voisinage(i_1_voisinage,j_1_voisinage);
			%fill(j_1_limites,i_1_limites,'g');
			plot(j_1,i_1, 'r+', 'MarkerSize', 30, 'LineWidth', 2);
			title("Localisation du voisinage sur l'image 1");
			hold off;

			% Voisinage 2
			subplot(2,2,4);
			imshow(I_2)
			axis on
			hold on;
			[i_2_limites, j_2_limites] = limites_voisinage(i_2_voisinage_re,j_2_voisinage_re);
			%fill(j_2_limites,i_2_limites,'g');
			plot(j_2,i_2, 'r+', 'MarkerSize', 30, 'LineWidth', 2);
			title("Localisation du voisinage sur l'image 2");
			hold off;
		end

		figure()
		imshow(I_1);
		axis on
		hold on;
		i_1_voisinage = i_1-rayon_voisinage:i_1+rayon_voisinage;
		j_1_voisinage = j_1-rayon_voisinage:j_1+rayon_voisinage;
		[i_1_voisinage, j_1_voisinage] = meshgrid(i_1_voisinage,j_1_voisinage);
		[i_1_limites, j_1_limites] = limites_voisinage(i_1_voisinage,j_1_voisinage);
		%fill(j_1_limites,i_1_limites,'g');

		P_arrow = P_1 + normale_theorique/100;
		u_arrow = round(P_2(1)/pixelSize);
		v_arrow = round(P_2(2)/pixelSize);
		i_arrow = v_arrow + v_0;
		j_arrow = u_arrow + u_0;
		i_arrow = i_1 + normale_theorique(2)
		j_arrow = j_1 + normale_theorique(1)

		for i_base = 
		X = [j_1 j_arrow];
		Y = [i_1 i_arrow];
		%annotation('arrow',X,Y,'color','r')
		plot(X,Y, 'b-', 'MarkerSize', 30, 'LineWidth', 2);

		plot(j_1,i_1, 'r+', 'MarkerSize', 30, 'LineWidth', 2);
		title("Localisation du voisinage sur l'image 1");
		hold off;
		%axis equal

	end
	disp("Appuyez sur une touche pour recommencer.")
	pause;
	close all;
end

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
