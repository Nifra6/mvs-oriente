%% Clear
clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

%% Données
load ../../data/donnees_calotte;
% Indices des images
indice_premiere_image = 1;
indice_deuxieme_image = 4;
% Les images
I_1 = I(:,:,indice_premiere_image);
I_2 = I(:,:,indice_deuxieme_image);
% Les tailles
[nombre_lignes, nombre_colonnes] = size(I_1);
% Les masques des images
masque_1 = masque(:,:,indice_premiere_image);
masque_2 = masque(:,:,indice_deuxieme_image);
% La pose
R_1_2 = R(:,:,indice_deuxieme_image) * R(:,:,indice_premiere_image)';
t_1_2 = t(:,indice_deuxieme_image) - R_1_2 * t(:,indice_premiere_image);
% Le gradient de l'image 2
[dx_I_1, dy_I_1] = gradient(I_1);
[dx_I_2, dy_I_2] = gradient(I_2);

%% Paramètres
valeurs_z 			= 80:.1:140;	% Les valeurs de profondeurs testées
rayon_voisinage		= 4;			% Voisinage carré à prendre en compte
affichage_log		= 1;			% Affichage d'informations diverses
interpolation		= 'nearest';	% Type d'interpolation
seuil_denominateur	= 0;			% Seuil pour accepter la division

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
	%i_1 = 199;
	%j_1 = 142;
	grad_I_1	= [dx_I_1(i_1,j_1); dy_I_1(i_1,j_1)];

	% Récupération de la profondeur
	for indice_z = 1:nb_profondeurs
		z = valeurs_z(indice_z);
		profondeur_reelle = round(z*10) == round(Z_1(i_1,j_1)*10);

		% Changements de repère
		u_1 = j_1 - u_0;
		v_1 = i_1 - u_0;
		P_1	= [u_1; v_1; z];
		P_2 = R_1_2 * P_1 + t_1_2;
		u_2 = P_2(1);
		v_2 = P_2(2);
		i_2 = v_2 + v_0;
		j_2 = u_2 + u_0;
		if (affichage_log && profondeur_reelle)
			disp("===== Changement de repère")
			[i_1 j_1]
			[i_2 j_2]
		end

		% Vérification si pixel hors image
		condition_image = i_2 > 0 & i_2 <= nombre_lignes & j_2 > 0 & j_2 <= nombre_colonnes;

		% Si le point reprojeté tombe sur le masque de la deuxième image
		if (condition_image && masque_2(round(i_2),round(j_2)))

			if (affichage_log && profondeur_reelle)
				disp("===== On est dans l'image")
			end
			grad_I_2 		= [interp2(dx_I_2,j_2,i_2); interp2(dy_I_2,j_2,i_2)];
			numerateur_pq 	= grad_I_1 - R_1_2(1:2,1:2)' * grad_I_2;
			denominateur_pq = R_1_2(1:2,3)' * grad_I_2;

			% Si pas de division par 0, on continue
			if (abs(denominateur_pq) > seuil_denominateur)

				% Estimation de la pente
				p_estime = numerateur_pq(1) / denominateur_pq;
				q_estime = numerateur_pq(2) / denominateur_pq;

				% Calcul du plan au pixel considéré
				normale = (1 / sqrt(p_estime^2 + q_estime^2 + 1)) * [p_estime ; q_estime ; -1];
				if (affichage_log && profondeur_reelle)
					disp("===== Comparaison des normales")
					normale_theorique = reshape(N_1(i_1,j_1,:),3,1);
					normale_theorique - normale;
					angle_normales = (180/pi) * atan2(norm(cross(normale_theorique,normale)),dot(normale_theorique,normale))
				end
				d_equation_plan = -P_1' * normale;

				% Calcul de la transformation géométrique
				u_1_decales = j_1-u_0-rayon_voisinage:j_1-u_0+rayon_voisinage;
				v_1_decales = i_1-v_0-rayon_voisinage:i_1-v_0+rayon_voisinage;
				[u_1_decales, v_1_decales] = meshgrid(u_1_decales,v_1_decales);
				z_1_decales = -(d_equation_plan + normale(1) * u_1_decales(:) + normale(2) * v_1_decales(:)) / normale(3);
				if (affichage_log && profondeur_reelle)
					disp("===== Profondeurs z")
					z
					reshape(z_1_decales, 2*rayon_voisinage+1, 2*rayon_voisinage+1)
				end

				% Reprojection du voisinage
				P_1_voisinage = [u_1_decales(:)' ; v_1_decales(:)' ; z_1_decales'];
				P_1_voisinage_mvs = [u_1_decales(:)' ; v_1_decales(:)' ; z * ones(1,size(u_1_decales,1)^2)];
				P_2_voisinage = R_1_2 * P_1_voisinage + t_1_2;
				P_2_voisinage_mvs = R_1_2 * P_1_voisinage_mvs + t_1_2;
				i_2_voisinage = P_2_voisinage(2,:) + v_0;
				i_2_voisinage_mvs = P_2_voisinage_mvs(2,:) + v_0;
				j_2_voisinage = P_2_voisinage(1,:) + u_0;
				j_2_voisinage_mvs = P_2_voisinage_mvs(1,:) + u_0;
				if (affichage_log && profondeur_reelle)
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
				I_2_voisinage_mvs = reshape(interp2(I_2, j_2_voisinage_mvs(:), i_2_voisinage_mvs(:),interpolation),2*rayon_voisinage+1,2*rayon_voisinage+1);
				if (affichage_log && profondeur_reelle)
					disp("===== Les images récupérées")
					I_1(i_1-rayon_voisinage:i_1+rayon_voisinage,j_1-rayon_voisinage:j_1+rayon_voisinage)
					I_2_voisinage
					indices_2_voisi = sub2ind([nombre_lignes nombre_colonnes], round(i_2_voisinage), round(j_2_voisinage));
					I_2_voisinage_nearest = reshape(I_2(indices_2_voisi),2*rayon_voisinage+1,2*rayon_voisinage+1);
					disp("===== Différences entre les images")
					diff = I_1(i_1-rayon_voisinage:i_1+rayon_voisinage,j_1-rayon_voisinage:j_1+rayon_voisinage) - I_2_voisinage
					sum(diff,"all")
				end

				% Calcul de l'erreur
				I_1_voisinage = I_1(i_1-rayon_voisinage:i_1+rayon_voisinage,j_1-rayon_voisinage:j_1+rayon_voisinage);
				%I_2_voisinage_mvs = I_2(i_2-rayon_voisinage:i_2+rayon_voisinage,j_2-rayon_voisinage:j_2+rayon_voisinage);
				erreurs_mvs(indice_z) = (1/(2*rayon_voisinage+1)^2) * sum((I_1_voisinage - I_2_voisinage_mvs).^2,'all');
				erreurs_sfs(indice_z) = (1/(2*rayon_voisinage+1)^2) * sum((I_1_voisinage - I_2_voisinage).^2,'all');
				erreurs_pq(indice_z) = (I_1(i_1,j_1) - 1 / sqrt(p_estime^2 + q_estime^2 + 1))^2;


				if (round(z) == round(Z_1(i_1,j_1)) & 0)
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
		end
	end

	% Afficher les erreurs
	figure;
	plot(valeurs_z',erreurs_mvs,'g')
	hold on
	plot(valeurs_z',erreurs_sfs,'r')
	plot([Z_1(i_1,j_1) ; Z_1(i_1,j_1)],[0 ; max(erreurs_mvs)],'b')
	legend('MVS Voisinage', 'SfS Voisinage')

	% Meilleures profondeur
	erreurs_mvs_corrige = (erreurs_mvs ~= 0) .* erreurs_mvs + (erreurs_mvs == 0) .* ones(size(erreurs_mvs));
	erreurs_sfs_corrige = (erreurs_sfs ~= 0) .* erreurs_sfs + (erreurs_sfs == 0) .* ones(size(erreurs_sfs));
	[~, indice_mvs] = min(erreurs_mvs_corrige);
	[~, indice_sfs] = min(erreurs_sfs_corrige);
	disp("==============================");
	disp("Valeur réelle :")
	Z_1(i_1,j_1)
	disp("Valeur MVS :")
	valeurs_z(indice_mvs)
	disp("Valeur SFS :")
	valeurs_z(indice_sfs)

	plot([valeurs_z(indice_mvs) ; valeurs_z(indice_mvs)],[0 ; max(erreurs_mvs)],'g')
	plot([valeurs_z(indice_sfs) ; valeurs_z(indice_sfs)],[0 ; max(erreurs_mvs)],'r')
	hold off

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
