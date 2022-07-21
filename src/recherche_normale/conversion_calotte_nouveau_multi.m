%% Chargement des données
clear;
load donnees_calotte_tr_multi;

%% Mise en forme des données
% Les nombres
[nombre_lignes, nombre_colonnes, nombre_images] = size(I_tr);

% Les caractéristiques de la caméra
u_0 = C_x_tr;
v_0 = C_y_tr;
facteur_k = 1;

% Image de référence
indice_image_reference = 1;

% L'éclairage
s = [0 ; 0 ; 1];

% Les images
I = I_tr;

% Les masques
masque = masque_tr;

% Les poses
R = zeros(3,3,nombre_images);
R(:,:,1) = eye(3);
R(:,:,2:end) = R_tr;
t = t_new;

% Les profondeurs
z = z_tr;

% Les normales
N = zeros(nombre_lignes,nombre_colonnes,3,nombre_images);
N(:,:,:,1) = N_1_tr;
%N(:,:,:,2) = ?
%N(:,:,:,3) = ?


%% Sauvegardes des données exploitables
save('../../data/simulateur_formate.mat','nombre_images', 'nombre_lignes', 'nombre_colonnes', 'indice_image_reference', 's', 'u_0', 'v_0', 'facteur_k', 'I', 'masque', 'R', 't', 'z', 'N');
