%% Chargement des données
clear;
load donnees_calotte_tr;

%% Mise en forme des données
% Les nombres
[nombre_lignes, nombre_colonnes] = size(I_1_tr);
nombre_images = 3;

% Les caractéristiques de la caméra
u_0 = C_x_tr;
v_0 = C_y_tr;
facteur_k = 1;

% Image de référence
indice_image_reference = 1;

% L'éclairage
s = [0 ; 0 ; 1];

% Les images
I = zeros(nombre_lignes, nombre_colonnes, nombre_images);
I(:,:,1) = I_1_tr;
I(:,:,2) = I_2_tr;
I(:,:,3) = I_3_tr;

% Les masques
masque = ones(size(I));
masque(:,:,1) = masque_1_tr;
masque(:,:,2) = masque_2_tr;
masque(:,:,3) = masque_3_tr;

% Les poses
R = zeros(3,3,nombre_images);
R(:,:,1) = eye(3);
R(:,:,2) = R_2_tr;
R(:,:,3) = R_3_tr;
t = zeros(3,nombre_images);
t(:,1) = t_1_new;
t(:,2) = t_2_new;
t(:,3) = t_3_new;

% Les profondeurs
z = zeros(size(I));
z(:,:,1) = Z_1_tr;
%z(:,:,2) = ?
%z(:,:,3) = ?

% Les normales
N = zeros(nombre_lignes,nombre_colonnes,3,nombre_images);
N(:,:,:,1) = N_1_tr;
%N(:,:,:,2) = ?
%N(:,:,:,3) = ?


%% Sauvegardes des données exploitables
save('../../data/simulateur_formate.mat','nombre_images', 'nombre_lignes', 'nombre_colonnes', 'indice_image_reference', 's', 'u_0', 'v_0', 'facteur_k', 'I', 'masque', 'R', 't', 'z', 'N');
