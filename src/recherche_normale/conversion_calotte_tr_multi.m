%% Chargement de la calotte mal conventionnée
clear;
load donnees_calotte_multi.mat;

%% Modification de la convention
% Les caractéristiques de la caméra
C_x_tr = C_y;
C_y_tr = C_x;

% Les images
I_tr = zeros(size(I));
for k = 1:nombre_images
	I_tr(:,:,k) = I(:,:,k)';
end

% Les coordonnées
X_tr = Y;
Y_tr = X;
dist_cam2cal = 200;
z_tr = zeros(size(z));
for k = 1:nombre_images
	z_tr(:,:,k) = dist_cam2cal - z(:,:,k)';
end

% Les masques
masque_tr = zeros(size(masque));
for k = 1:nombre_images
	masque_tr(:,:,k) = masque(:,:,k)';
end

% Les poses qui marchent pour le cas sans décalage
R_tr = R;
t_new = repmat([0 ; 0 ; dist_cam2cal],1,nombre_images);
t_tr = zeros(3,nombre_images-1);
for k = 1:nombre_images-1
	t_tr(:,k) = t_new(:,k+1) - R_tr(:,:,k) * t_new(:,1); 
end

% Les normales
N_1_tr = permute(N_1,[2 1 3]);
N_1_tr(:,:,3) = -N_1_tr(:,:,3);

save('donnees_calotte_tr_multi.mat');
