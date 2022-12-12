% Bruiter les images issues du simulateur avec un bruit blanc uniforme,
% ayant une probabilité réglable de survenir par pixel.


%% Paramètres
snr = 6;		% Force du bruit en 8 bit (à régler entre 0 et 255)

%% Chargement des données
load ../../data/perspectif/simulateur_formate.mat;

%% Bruitage des données
bruit_g = (snr/255) * (rand(size(I))-0.5);
I = I + bruit_g;
I = (I <= 1) .* I + (I > 1) .* 1;
I = (I >= 0) .* I + (I < 0) .* 0;


%% Sauvegardes des données bruitées
clear snr bruit_g
save('../../data/perspectif/simulateur_formate.mat');
