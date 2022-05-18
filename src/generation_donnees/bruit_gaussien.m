% Bruiter les images issues du simulateur avec un bruit de type poivre et sel,
% ayant une probabilité réglable de survenir par pixel.

%% Chargement des données
load ../../data/simulateur_formate.mat;

%% Bruitage des données
snr = 10;
bruit_g = (snr/255) * rand(size(I));
I = I + bruit_g;


%% Sauvegardes des données bruitées
save('../../data/simulateur_formate.mat','nombre_images', 'nombre_lignes', 'nombre_colonnes', 'indice_image_reference', 's', 'K', 'u_0', 'v_0', 'f', 'facteur_k', 'I', 'masque', 'R', 't', 'z', 'N');
