% Convertir les données du simulateur lambertien (https://github.com/bbrument/lambertianRendering_v1)
% vers des données directement exploitables par nos codes.
% Il est possible d'enlever les points qui ne sont pas visibles dans toutes les images pour
% analyser les résultats avec la vérité terrain.

%% Paramètres
filtrage_points	= 1;	% Retire tous les points 3D qui ne sont pas reprojeté dans chacune des images
masque_calotte	= 0;

%% Chargement des données
load ../../data/perspectif/simulateur.mat;

%% Mise en forme des données
% Les nombres
[nb_lignes, nb_colonnes, nb_canaux, nb_images] = size(renderedImages);

% Les caractéristiques de la caméra
K = params.K;
u_0 = K(1,3);
v_0 = K(2,3);
f = K(1,1);

% Image de référence
indice_image_reference = 5;

% L'éclairage
s = params.lightIntensity * params.lightSource;

% Les images
I_to_sort = squeeze(renderedImages) / max(renderedImages,[],'all');
I = zeros(size(I_to_sort));
I(:,:,1) = I_to_sort(:,:,indice_image_reference);
I(:,:,2:indice_image_reference) = I_to_sort(:,:,1:indice_image_reference-1);
I(:,:,indice_image_reference+1:end) = I_to_sort(:,:,indice_image_reference+1:end);

% Les masques
masque = ones(size(I));
if (masque_calotte)
	masque = I > 0;
end

% Les poses
R_to_sort = params.w2cPoses(:,1:3,:);
t_to_sort = squeeze(params.w2cPoses(:,4,:));
R = zeros(size(R_to_sort));
R(:,:,1) = R_to_sort(:,:,indice_image_reference);
R(:,:,2:indice_image_reference) = R_to_sort(:,:,1:indice_image_reference-1);
R(:,:,indice_image_reference+1:end) = R_to_sort(:,:,indice_image_reference+1:end);
t = zeros(size(t_to_sort));
t(:,1) = t_to_sort(:,indice_image_reference);
t(:,2:indice_image_reference) = t_to_sort(:,1:indice_image_reference-1);
t(:,indice_image_reference+1:end) = t_to_sort(:,indice_image_reference+1:end);


% Les profondeurs
%z_to_sort = distMaps;
z_to_sort = depthMaps;
z = zeros(size(z_to_sort));
z(:,:,1) = z_to_sort(:,:,indice_image_reference);
z(:,:,2:indice_image_reference) = z_to_sort(:,:,1:indice_image_reference-1);
z(:,:,indice_image_reference+1:end) = z_to_sort(:,:,indice_image_reference+1:end);


% Les normales
N_to_sort = normalMaps;
N = zeros(size(N_to_sort));
N(:,:,:,1) = N_to_sort(:,:,:,indice_image_reference);
N(:,:,:,2:indice_image_reference) = N_to_sort(:,:,:,1:indice_image_reference-1);
N(:,:,:,indice_image_reference+1:end) = N_to_sort(:,:,:,indice_image_reference+1:end);


%% Optionel - Retirer les points masqués dans certaines images
if (filtrage_points)
	grille_pixels = 2;
	rayon_voisinage = 15;
	offset = 0.5;
	taille_patch = (2 * rayon_voisinage + 1)^2;
	liste_z_a_regarder = [min(z(:,:,1),[],'all'), max(z(:,:,1),[],'all')];
	K_inv = inv(K);

	for i_z = 1:size(liste_z_a_regarder,2)
		z_a_regarder = liste_z_a_regarder(i_z); 

		% Calcul des positions 3D des points dans l'image 1
		masque_1 = masque(:,:,1);
		[i_k, j_k] = find(masque_1);
		ind_1 = sub2ind([nb_lignes nb_colonnes], i_k, j_k);
		if (grille_pixels > 0)
			indices_grilles = (mod(i_k,grille_pixels) == 1) & (mod(j_k,grille_pixels) == 1);
			ind_1 = ind_1(find(indices_grilles));
			i_k = i_k(find(indices_grilles));
			j_k = j_k(find(indices_grilles));
		end
		nb_pixels_etudies = size(ind_1,1);
		%Z_1 = z(:,:,1);
		p_1 = [j_k - offset , i_k - offset , ones(size(i_k))]';
		P_k = zeros(3, nb_pixels_etudies, nb_images);
		P_k(:,:,1) = z_a_regarder * (K_inv * p_1);

		% Les poses relatives
		R_1_k = zeros(3,3,nb_images-1);
		t_1_k = zeros(3,nb_images-1);
		for k = 1:nb_images-1
			R_1_k(:,:,k) = R(:,:,k+1) * R(:,:,1)';
			t_1_k(:,k) = t(:,k+1) - R_1_k(:,:,k) * t(:,1);
		end

		% Création du patch
		normale = [N(ind_1)' ; N(ind_1+nb_lignes*nb_colonnes)' ; N(ind_1+2*nb_lignes*nb_colonnes)'];
		voisinage_ligne = -rayon_voisinage*nb_lignes:nb_lignes:rayon_voisinage*nb_lignes;
		voisinage_colonne = -rayon_voisinage:rayon_voisinage;
		grille_voisinage = voisinage_ligne + voisinage_colonne';
		grille_voisinage = grille_voisinage';

		% Calcul du plan considéré
		d_equation_plan = sum(-P_k(:,:,1) .* normale,1);

		% Calcul de la transformation géométrique
		ind_decales = ind_1 + grille_voisinage(:)'; % Création de matrice avec 2 vecteurs
		[i_1_decales, j_1_decales] = ind2sub([nb_lignes, nb_colonnes], ind_decales);
		u_1_decales = j_1_decales - offset;
		v_1_decales = i_1_decales - offset;

		% Reprojection du voisinage
		i_2_voisinage = zeros(nb_pixels_etudies, taille_patch, nb_images-1);
		j_2_voisinage = zeros(nb_pixels_etudies, taille_patch, nb_images-1);
		u_1_decales_vec = reshape(u_1_decales',1,nb_pixels_etudies*taille_patch);
		v_1_decales_vec = reshape(v_1_decales',1,nb_pixels_etudies*taille_patch);
		for k = 1:nb_images-1
			for pixel = 1:nb_pixels_etudies
				homographie = K * (R_1_k(:,:,k) - t_1_k(:,k) * normale(:,pixel)' / d_equation_plan(pixel)) * K_inv;
				p_2_voisinage = homographie * z_a_regarder * [u_1_decales(pixel,:) ; v_1_decales(pixel,:) ; ones(1,taille_patch)];
				u_2_voisinage = p_2_voisinage(1,:) ./ p_2_voisinage(3,:);
				v_2_voisinage = p_2_voisinage(2,:) ./ p_2_voisinage(3,:);
				i_2_voisinage(pixel,:,k) = v_2_voisinage + offset;
				j_2_voisinage(pixel,:,k) = u_2_voisinage + offset;
			end
		end

		% Vérification des pixels hors images
		condition_image = prod(prod( i_2_voisinage > 0.5 & i_2_voisinage <= nb_lignes & j_2_voisinage > 0.5 & j_2_voisinage <= nb_colonnes ,3),2);
		condition_image = ones(nb_pixels_etudies,1);
		for k = 1:nb_images-1
			for l = 1:taille_patch
				condition_image = condition_image & i_2_voisinage(:,l,k) > 0.5 & i_2_voisinage(:,l,k) <= nb_lignes & j_2_voisinage(:,l,k) > 0.5 & j_2_voisinage(:,l,k) <= nb_colonnes;
			end
		end

		masque(:,:,1) = 0;
		masque(ind_1) = condition_image;
	end
end


%% Sauvegardes des données exploitables
save('../../data/perspectif/simulateur_formate.mat','nb_images', 'nb_lignes', 'nb_colonnes', 'indice_image_reference', 's', 'K', 'u_0', 'v_0', 'f', 'I', 'masque', 'R', 't', 'z', 'N');
