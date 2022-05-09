% Bruiter les images issues du simulateur avec un bruit de type poivre et sel,
% ayant une probabilité réglable de survenir par pixel.

%% Chargement des données
load ../../data/simulateur.mat;

%% Bruitage des données
seuil = 0.02;
tirage = rand(size(I));
bruit_poivre = tirage < seuil;
bruit_sel = tirage > (1 - seuil);
pas_bruit = ones(size(I)) - bruit_sel - bruit_poivre;
I = pas_bruit .* I + bruit_poivre .* zeros(size(I)) + bruit_sel .* ones(size(I));


%% Sauvegardes des données bruitées
save('../../data/simulateur_formate.mat','nombre_images', 'nombre_lignes', 'nombre_colonnes', 'indice_image_reference', 'K', 'u_0', 'v_0', 'f', 'facteur_k', 'I', 'masque', 'R', 't', 'z', 'N');
