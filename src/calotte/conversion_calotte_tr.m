%% Chargement de la calotte mal conventionnée
clear;
load donnees_calotte.mat;

%% Modification de la convention
% Les caractéristiques de la caméra
C_x_tr = C_y;
C_y_tr = C_x;

% Les images
I_1_tr = I_1';
I_2_tr = I_2';
I_3_tr = I_3';

% Les coordonnées
X_tr = Y;
Y_tr = X;
dist_cam2cal = 200;
Z_1_tr = dist_cam2cal - Z_1';

% Les masques
masque_1_tr = masque_1';
masque_2_tr = masque_2';
masque_3_tr = masque_3';

% Les poses qui marchent pour le cas sans décalage
R_2_tr = R_2;
R_3_tr = R_3;
t_1_new = [0 ; 0 ; dist_cam2cal];
t_2_new = [0 ; 0 ; dist_cam2cal];
t_3_new = [0 ; 0 ; dist_cam2cal];
t_2_tr = t_2_new - R_2_tr * t_1_new;
t_3_tr = t_3_new - R_3_tr * t_1_new;

% Les normales
N_1_tr = permute(N_1,[2 1 3]);
N_1_tr(:,:,3) = -N_1_tr(:,:,3);

save('donnees_calotte_tr.mat');
