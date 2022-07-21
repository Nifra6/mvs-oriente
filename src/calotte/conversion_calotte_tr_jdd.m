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
I_1_tr = I_1;
I_2_tr = I_2;

% Les coordonnées
X_tr = X_1;
Y_tr = Y_1;
Z_1_tr = Z_1;
dist_cam2cal = t;

% Les masques
masque_1_tr = masque_1;
masque_2_tr = masque_2;

% Les poses qui marchent pour le cas sans décalage
R_2_tr = Rotation';
t_1_new = O_1;
t_2_new = O_2;
t_2_tr = t_2_new - R_2_tr * t_1_new;

% Les normales
N_1_tr = N_1;

save('donnees_calotte_tr.mat');
