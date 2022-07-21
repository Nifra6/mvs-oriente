% Bruiter les images issues du simulateur avec un bruit blanc uniforme,
% ayant une probabilité réglable de survenir par pixel.

%% Chargement des données
load ../../data/simulateur_formate.mat;

%% Paramètres
snr = 6;		% Force du bruit en 8 bit (à régler entre 0 et 255)

%% Bruitage des données
bruit_g = (snr/255) * (rand(size(I))-0.5);
I = I + bruit_g;
I = (I <= 1) .* I + (I > 1) .* 1;
I = (I >= 0) .* I + (I < 0) .* 0;


%% Sauvegardes des données bruitées
save('../../data/simulateur_formate.mat','nombre_images', 'nombre_lignes', 'nombre_colonnes', 'indice_image_reference', 's', 'K', 'u_0', 'v_0', 'f', 'facteur_k', 'I', 'masque', 'R', 't', 'z', 'N');
