%% Clear
clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

%% Données
load ../../data/donnees_calotte;
% Les images
I_1 = I(:,:,1);
I_2 = I(:,:,2);
% Les masques des images
masque_1 = masque(:,:,1);
masque_2 = masque(:,:,2);
% La pose
R_2 = R(:,:,1);
t_2 = t(:,1);
% Le gradient de l'image 2
%dx_I_2 = dx_I(:,:,2);
%dy_I_2 = dy_I(:,:,2);
[dy_I_1, dx_I_1] = gradient(I_1);
[dy_I_2, dx_I_2] = gradient(I_2);

%% Paramètres
valeurs_z 		= 60:.1:120;	% Les valeurs de profondeurs utilisées
range			= 4;		% Voisinage à prendre en compte
affichage_log	= 0;		% Affichage d'informations diverses

%% Algorithme
nb_profondeurs = size(valeurs_z,2);
erreurs_mvs	= zeros(nb_profondeurs,1);
erreurs_sfs	= zeros(nb_profondeurs,1);
erreurs_pq = zeros(nb_profondeurs,1);

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
	for indice_z = 1:nb_profondeurs
		z = valeurs_z(indice_z);

		% Changements de repère
		P_1	= [i_1 - u_0; j_1 - v_0; z];
		P_2 = R_2 * P_1;
		i_2 = round(P_2(1) + u_0);
		j_2 = round(P_2(2) + v_0);

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
				p_estime = -numerateur_pq(1) / denominateur_pq;
				q_estime = -numerateur_pq(2) / denominateur_pq;

				% Calcul du plan au pixel considéré
				normale = (1 / (p_estime^2 + q_estime^2 + 1)) * [p_estime ; q_estime ; 1];
				if (affichage_log)
					disp("===== Comparaison des normales")
					normale_theorique = reshape(N_1(i_1,j_1,:),3,1);
					normale_theorique - normale
					(180/pi) * atan2(norm(cross(normale_theorique,normale)),dot(normale_theorique,normale))
					(180/pi) * acos(dot(normale_theorique,normale)/(norm(normale_theorique)*norm(normale)))
				end
				d_equation_plan = -P_1' * normale;

				% Calcul de la transformation géométrique
				i_1_decales = i_1-u_0-range:i_1-u_0+range;
				j_1_decales = j_1-v_0-range:j_1-v_0+range;
				[i_1_decales, j_1_decales] = meshgrid(i_1_decales,j_1_decales);
				z_1_decales = -(d_equation_plan + normale(1) * i_1_decales(:) + normale(2) * j_1_decales(:)) / normale(3);
				if (affichage_log)
					disp("===== Profondeurs z")
					z
					reshape(z_1_decales, 2*range+1, 2*range+1)
				end

				% Reprojection du voisinage
				P_1_voisinage = [i_1_decales(:)' ; j_1_decales(:)' ; z_1_decales'];
				P_1_voisinage_mvs = [i_1_decales(:)' ; j_1_decales(:)' ; z * ones(1,size(i_1_decales,1)^2)];
				P_2_voisinage = R_2 * P_1_voisinage;
				P_2_voisinage_mvs = R_2 * P_1_voisinage_mvs;
				i_2_voisinage = round(P_2_voisinage(1,:) + u_0);
				i_2_voisinage_mvs = round(P_2_voisinage_mvs(1,:) + u_0);
				j_2_voisinage = round(P_2_voisinage(2,:) + v_0);
				j_2_voisinage_mvs = round(P_2_voisinage_mvs(2,:) + v_0);
				if (affichage_log)
					disp("===== Les i du voisinage")
					i_2
					i_2_voisinage_re = reshape(i_2_voisinage, 2*range+1, 2*range+1)
					disp("===== Les j du voisinage")
					j_2
					j_2_voisinage_re = reshape(j_2_voisinage, 2*range+1, 2*range+1)
					disp("===== Le contour du voisinage")
					i_2_voisinage_re
					[i_2_limites, j_2_limites] = limites_voisinage(i_2_voisinage_re+u_0,j_2_voisinage_re+v_0);
					i_2_limites
				end

				% Récupération des niveaux de gris dans l'image 2 du voisinage	
				I_2_voisinage = reshape(interp2(I_2, j_2_voisinage(:), i_2_voisinage(:),'nearest'),2*range+1,2*range+1)';
				I_2_voisinage_mvs = reshape(interp2(I_2, j_2_voisinage_mvs(:), i_2_voisinage_mvs(:),'nearest'),2*range+1,2*range+1)';
				if (affichage_log)
					disp("===== Différences entre les images")
					diff = I_1(i_1-range:i_1+range,j_1-range:j_1+range) - I_2_voisinage
					sum(diff,"all")
				end

				% Calcul de l'erreur
				I_1_voisinage = I_1(i_1-range:i_1+range,j_1-range:j_1+range);
				%I_2_voisinage_mvs = I_2(i_2-range:i_2+range,j_2-range:j_2+range);
				erreurs_mvs(indice_z) = (1/(2*range+1)^2) * sum((I_1_voisinage - I_2_voisinage_mvs).^2,'all');
				erreurs_sfs(indice_z) = (1/(2*range+1)^2) * sum((I_1_voisinage - I_2_voisinage).^2,'all');
				erreurs_pq(indice_z) = (I_1(i_1,j_1) - 1 / sqrt(p_estime^2 + q_estime^2 + 1))^2;


				if (round(z) == round(Z_1(i_1,j_1)) & 0)
					% Affichage des résultats
					figure;

					% Image 1
					subplot(2,2,1);
					imshow(I_1(i_1-range:i_1+range,j_1-range:j_1+range));
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
					i_1_voisinage = i_1-range:i_1+range;
					j_1_voisinage = j_1-range:j_1+range;
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
		end
	end

	% Afficher les erreurs
	figure;
	plot(valeurs_z',erreurs_mvs,'g')
	hold on
	plot(valeurs_z',erreurs_sfs,'r')
	%plot(valeurs_z',erreurs_pq,'b')
	plot([Z_1(i_1,j_1) ; Z_1(i_1,j_1)],[0 ; max(erreurs_mvs)],'b')
	legend('MVS Voisinage', 'SfS Voisinage')
	hold off

	% Meilleures profondeur
	[~, indice_mvs] = min(erreurs_mvs);
	[~, indice_sfs] = min(erreurs_sfs);
	disp("==============================");
	disp("Valeur réelle :")
	Z_1(i_1,j_1)
	disp("Valeur MVS :")
	valeurs_z(indice_mvs)
	disp("Valeur SFS :")
	valeurs_z(indice_sfs)

	% Demander le nouveau pixel
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
