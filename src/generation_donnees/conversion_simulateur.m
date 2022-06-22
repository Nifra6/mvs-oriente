% Convertir les données du simulateur lambertien (https://github.com/bbrument/lambertianRendering_v1)
% vers des données directement exploitables par nos codes.

%% Chargement des données
load ../../data/simulateur.mat;

%% Mise en forme des données
% Les nombres
[nombre_lignes, nombre_colonnes, nombre_canaux, nombre_images] = size(renderedImages);

% Les caractéristiques de la caméra
K = params.K;
u_0 = K(1,3);
v_0 = K(2,3);
f = K(1,1);
facteur_k = 451*(4/3^2);

% Image de référence
indice_image_reference = 2;

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


%% Retirer les points masqués dans certaines images
% Calcul des positions 3D des points dans l'image 1
masque_1 = masque(:,:,1);
[i_k, j_k] = find(masque_1);
ind_1 = sub2ind([nombre_lignes nombre_colonnes], i_k, j_k);
Z_1 = z(:,:,1);
P_k = zeros(3, size(ind_1,1), nombre_images);
P_k(:,:,1) = [(j_k - u_0) / facteur_k, (i_k - v_0) / facteur_k, Z_1(:)].';

% Les poses relatives
R_1_k = zeros(3,3,nombre_images-1);
t_1_k = zeros(3,nombre_images-1);
for k = 1:nombre_images-1
	R_1_k(:,:,k) = R(:,:,k+1) * R(:,:,1)';
	t_1_k(:,k) = t(:,k+1) - R_1_k(:,:,k) * t(:,1);
end

% Déprojection et reprojection dans les autres images
for k = 1:nombre_images-1
	P_k(:,:,k+1) = R_1_k(:,:,k) * P_k(:,:,1) + t_1_k(:,k);
	i_k(:,k+1) = (P_k(2,:,k+1) * facteur_k + v_0).';
	j_k(:,k+1) = (P_k(1,:,k+1) * facteur_k + u_0).';
end

% Vérification des pixels hors images
rayon_voisinage = 3;
condition_image = ones(size(ind_1,1),1);
for k = 1:nombre_images-1
	condition_image = condition_image & i_k(:,k+1) > 0.5 + rayon_voisinage & i_k(:,k+1) <= nombre_lignes + 0.5 - rayon_voisinage & j_k(:,k+1) > 0.5 + rayon_voisinage & j_k(:,k+1) <= nombre_colonnes + 0.5 - rayon_voisinage;
end

masque(:,:,1) = reshape(condition_image, nombre_lignes, nombre_colonnes);


%% Sauvegardes des données exploitables
save('../../data/simulateur_formate.mat','nombre_images', 'nombre_lignes', 'nombre_colonnes', 'indice_image_reference', 's', 'K', 'u_0', 'v_0', 'f', 'facteur_k', 'I', 'masque', 'R', 't', 'z', 'N');
