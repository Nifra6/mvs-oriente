%% Chargement de la calotte mal conventionnée
clear;
load donnees_calotte.mat;

%% Modification de la convention
% Les caractéristiques de la caméra
C_x_tr = t_sur_2;
C_y_tr = t_sur_2;
taille = t;
delta_x = 0;
delta_y = 0;

% Les images
I_1_tr = I_k(:,:,1);
I_2_tr = I_k(:,:,2);

% Les coordonnées
X_tr = X_1;
Y_tr = Y_1;
Z_1_tr = Z_k(:,:,1);
dist_cam2cal = 2*R;

% Les masques
masque_1_tr = masque_k(:,:,1);
masque_2_tr = masque_k(:,:,2);

% Les poses qui marchent pour le cas sans décalage
R_1 = R_k(:,:,1)';
R_2 = R_k(:,:,2)';
R_2_tr = R_2 * R_1';
t_1_new = - R_1 * t_k(:,1);
t_2_new = - R_2 * t_k(:,2);
t_2_tr = t_2_new - R_2_tr * t_1_new;

% Les normales
N_1_tr = N_1;

save('donnees_calotte_tr.mat');
