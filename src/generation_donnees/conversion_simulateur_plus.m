load ../../data/simulateur.mat;

% Les nombres
[nombre_lignes, nombre_colonnes, nombre_canaux, nombre_images] = size(renderedImages);

% Les caractéristiques de la caméra
K = params.K;
u_0 = K(1,3);
v_0 = K(2,3);
f = K(1,1);
facteur_k = 451*(4/3^2);

% Image de référence
indice_image_reference = 5;

% Les images
I_to_sort = squeeze(renderedImages) / max(renderedImages,[],'all');
size(I_to_sort)
I = zeros(size(I_to_sort));
I(:,:,1) = I_to_sort(:,:,indice_image_reference);
I(:,:,2:indice_image_reference) = I_to_sort(:,:,1:indice_image_reference-1);
I(:,:,indice_image_reference+1:end) = I_to_sort(:,:,indice_image_reference+1:end);

% Les masques
masque = ones(size(I));

% Les poses
R_to_sort = params.w2cPoses(:,1:3,:);
size(R_to_sort)
t_to_sort = squeeze(params.w2cPoses(:,4,:));
size(t_to_sort)
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
size(z_to_sort)
z = zeros(size(z_to_sort));
z(:,:,1) = z_to_sort(:,:,indice_image_reference);
z(:,:,2:indice_image_reference) = z_to_sort(:,:,1:indice_image_reference-1);
z(:,:,indice_image_reference+1:end) = z_to_sort(:,:,indice_image_reference+1:end);


% Les normales
N = normalMaps;


save('../../data/simulateur_formate.mat','nombre_images', 'nombre_lignes', 'nombre_colonnes', 'indice_image_reference', 'K', 'u_0', 'v_0', 'f', 'facteur_k', 'I', 'masque', 'R', 't', 'z', 'N');
