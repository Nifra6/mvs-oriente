%% Clear
clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

%% Données
load ../../data/simulateur_formate.mat;
% Indices des images
indice_premiere_image = 1;
indice_deuxieme_image = 2;
% Les profondeurs
Z_1 = z(:,:,indice_premiere_image);
% Les images
I_1 = I(:,:,indice_premiere_image);
I_2 = I(:,:,indice_deuxieme_image);
% TODO Normales
N_1 = N(:,:,:,indice_premiere_image);
% Les masques
masque_1 = masque(:,:,indice_premiere_image);
masque_2 = masque(:,:,indice_deuxieme_image);
% La pose
R_1_2 = R(:,:,indice_deuxieme_image) * R(:,:,indice_premiere_image)';
t_1_2 = t(:,indice_deuxieme_image) - R_1_2 * t(:,indice_premiere_image);
% Le gradient de l'image 2
[dx_I_1, dy_I_1] = gradient(I_1);
[dx_I_2, dy_I_2] = gradient(I_2);

%% Paramètres
rayon_voisinage		= 1;			% Voisinage carré à prendre en compte
affichage_log		= 1;			% Affichage d'informations diverses
interpolation		= 'nearest';	% Type d'interpolation
seuil_denominateur	= 0;			% Seuil pour accepter la division
%facteur = 460 * (4/3) * (1/3) ;
facteur = 451 * (4/3) * (1/3) ;

%% Algorithme
while (1)
	% Sélection d'un pixel
	figure;
	title("Cliquez sur le pixel souhaité.")
	imshow(I_1);
	P			= drawpoint;
	pos 		= P.Position;
	i_1 		= round(pos(2));
	j_1 		= round(pos(1));
	grad_I_1	= [dx_I_1(i_1,j_1); dy_I_1(i_1,j_1)];

	% Récupération de la profondeur
	z = Z_1(i_1,j_1);

	% Changements de repère
	u_1 = j_1 - u_0;
	v_1 = i_1 - v_0;
	P_1	= [u_1 / facteur ; v_1 / facteur ; z];
	P_2 = R_1_2 * P_1 + t_1_2;
	u_2 = P_2(1) * facteur;
	v_2 = P_2(2) * facteur;
	i_2 = v_2 + v_0;
	j_2 = u_2 + u_0;
	if (affichage_log)
		disp("===== Changement de repère")
		[i_1 j_1]
		[i_2 j_2]
	end

	% Vérification si pixel hors image
	condition_image = i_2 > 0 & i_2 <= nombre_lignes & j_2 > 0 & j_2 <= nombre_colonnes;

	% Si le point reprojeté tombe sur le masque de la deuxième image
	if condition_image & masque_2(round(i_2),round(j_2))

		grad_I_2 		= [interp2(dx_I_2,j_2,i_2,interpolation); interp2(dy_I_2,j_2,i_2,interpolation)];
		numerateur_pq 	= grad_I_1 - R_1_2(1:2,1:2)' * grad_I_2;
		denominateur_pq = R_1_2(1:2,3)' * grad_I_2

		% Si pas de division par 0, on continue
		if (abs(denominateur_pq) > seuil_denominateur)

			% Estimation de la pente
			p_estime = numerateur_pq(1) / denominateur_pq;
			q_estime = numerateur_pq(2) / denominateur_pq;

			% Calcul du plan au pixel considéré
			normale = (1 / sqrt(p_estime^2 + q_estime^2 + 1)) * [p_estime ; q_estime ; -1];
			if (affichage_log)
				disp("===== Pente estimée")
				[p_estime q_estime]
				disp("===== Comparaison des normales")
				normale_theorique = reshape(N_1(i_1,j_1,:),3,1);
				normale_theorique - normale;
				angle_normales = (180/pi) * atan2(norm(cross(normale_theorique,normale)),dot(normale_theorique,normale))
			end
			d_equation_plan = -P_1' * normale;

			% Calcul de la transformation géométrique
			u_1_decales = (j_1-u_0-rayon_voisinage)/facteur:(1/facteur):(j_1-u_0+rayon_voisinage)/facteur;
			v_1_decales = (i_1-v_0-rayon_voisinage)/facteur:(1/facteur):(i_1-v_0+rayon_voisinage)/facteur;
			[u_1_decales, v_1_decales] = meshgrid(u_1_decales,v_1_decales);
			z_1_decales = -(d_equation_plan + normale(1) * u_1_decales(:) + normale(2) * v_1_decales(:)) / normale(3);
			if (affichage_log)
				disp("===== Profondeurs z")
				z
				reshape(z_1_decales, 2*rayon_voisinage+1, 2*rayon_voisinage+1)
			end

			% Reprojection du voisinage
			P_1_voisinage = [u_1_decales(:)' ; v_1_decales(:)' ; z_1_decales'];
			P_2_voisinage = R_1_2 * P_1_voisinage + t_1_2;
			i_2_voisinage = P_2_voisinage(2,:) * facteur + v_0;
			j_2_voisinage = P_2_voisinage(1,:) * facteur + u_0;
			if (affichage_log)
				disp("===== Les i du voisinage")
				i_2_voisinage_re = reshape(i_2_voisinage, 2*rayon_voisinage+1, 2*rayon_voisinage+1)
				disp("===== Les j du voisinage")
				j_2_voisinage_re = reshape(j_2_voisinage, 2*rayon_voisinage+1, 2*rayon_voisinage+1)
				disp("===== Le contour du voisinage")
				i_2_voisinage_re
				[i_2_limites, j_2_limites] = limites_voisinage(i_2_voisinage_re+u_0,j_2_voisinage_re+v_0);
				i_2_limites
			end

			% Récupération des niveaux de gris dans l'image 2 du voisinage	
			I_2_voisinage = reshape(interp2(I_2, j_2_voisinage(:), i_2_voisinage(:),interpolation),2*rayon_voisinage+1,2*rayon_voisinage+1);
			if (affichage_log)
				disp("===== Les images récupérées")
				I_1(i_1-rayon_voisinage:i_1+rayon_voisinage,j_1-rayon_voisinage:j_1+rayon_voisinage)
				I_2_voisinage
				indices_2_voisi = sub2ind([nombre_lignes nombre_colonnes], round(i_2_voisinage), round(j_2_voisinage));
				I_2_voisinage_nearest = reshape(I_2(indices_2_voisi),2*rayon_voisinage+1,2*rayon_voisinage+1);
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
			fill(j_1_limites,i_1_limites,'g');
			plot(j_1,i_1, 'r+', 'MarkerSize', 30, 'LineWidth', 2);
			title("Localisation du voisinage sur l'image 1");
			hold off;

			% Voisinage 2
			subplot(2,2,4);
			imshow(I_2)
			axis on
			hold on;
			[i_2_limites, j_2_limites] = limites_voisinage(i_2_voisinage_re,j_2_voisinage_re);
			fill(j_2_limites,i_2_limites,'g');
			plot(j_2,i_2, 'r+', 'MarkerSize', 30, 'LineWidth', 2);
			title("Localisation du voisinage sur l'image 2");
			hold off;
		end
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
