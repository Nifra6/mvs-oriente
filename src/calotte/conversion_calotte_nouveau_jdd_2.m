%% Chargement des données
clear;
load donnees_calotte_tr;

%% Mise en forme des données
% Les nombres
[nb_lignes, nb_colonnes] = size(I_1_tr);
nb_images = 2;

% Les caractéristiques de la caméra
u_0 = C_x_tr;
v_0 = C_y_tr;
facteur_k = 1;

% Image de référence
indice_image_reference = 1;

% L'éclairage
s = [0 ; 0 ; 1];

% Les images
I = zeros(nb_lignes, nb_colonnes, nb_images);
I(:,:,1) = I_1_tr;
I(:,:,2) = I_2_tr;
%I(:,:,3) = I_3_tr;

% Les masques
masque = ones(size(I));
masque(:,:,1) = masque_1_tr;
masque(:,:,2) = masque_2_tr;
%masque(:,:,3) = masque_3_tr;

% Les poses
R = zeros(3,3,nb_images);
R(:,:,1) = R_1;
R(:,:,2) = R_2;
t = zeros(3,nb_images);
t(:,1) = t_1_new;
t(:,2) = t_2_new;

% Les profondeurs
z = zeros(size(I));
z(:,:,1) = Z_1_tr;
%z(:,:,2) = ?
%z(:,:,3) = ?

% Les normales
N = zeros(nb_lignes,nb_colonnes,3,nb_images);
N(:,:,:,1) = N_1_tr;
%N(:,:,:,2) = ?
%N(:,:,:,3) = ?


%% Sauvegardes des données exploitables
save('../../data/simulateur_formate.mat','nb_images', 'nb_lignes', 'nb_colonnes', 'indice_image_reference', 's', 'u_0', 'v_0', 'facteur_k', 'I', 'masque', 'R', 't', 'z', 'N');
