%% Données
pixelSize = 2 / 1920;
s = [0 ; 0 ; -1];

% Première image
[image, colorMap] = imread('../../data/lapin_stanford/im1.png');

Im = ind2gray(image, colorMap);
mask = Im > 0;

R = eyes(3);
t = [0 ; 0 ; 2];

% Deuxième image

[image, colorMap] = imread('../../data/lapin_stanford/im2.png');

Im(:,:,2) = ind2gray(image, colorMap);
mask(:,:,2) = Im(:,:,2) > 0;

angle = (pi/180) * 7;
R(:,:,2) = [cos(angle) 0 sin(angle); 0 1 0 ; -sin(angle) 0 cos(angle)];
t(:,2) = [0.3 ; 0 ; 2];

save ../../data/lapin_ortho
