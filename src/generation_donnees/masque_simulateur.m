load '../../data/simulateur_formate.mat'

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

save('../../data/simulateur_formate.mat','nombre_images', 'nombre_lignes', 'nombre_colonnes', 'indice_image_reference', 'K', 'u_0', 'v_0', 'f', 'facteur_k', 'I', 'masque', 'R', 't', 'z', 'N');
