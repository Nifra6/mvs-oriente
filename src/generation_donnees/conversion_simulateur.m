load ../../data/simulateur.mat;

% Les nombres
[nombre_lignes, nombre_colonnes, nombre_canaux, nombre_images] = size(renderedImages);

% Les caractéristiques de la caméra
K = params.K;
u_0 = K(1,3);
v_0 = K(2,3);
f = K(1,1);

% Les images
I = squeeze(renderedImages) / max(renderedImages,[],'all');

% Les masques
masque = ones(size(I));

% Les poses
R = params.w2cPoses(:,1:3,:);
t = squeeze(params.w2cPoses(:,4,:));

% Les profondeurs
z = depthMaps;


save('../../data/simulateur_formate.mat','nombre_images', 'nombre_lignes', 'nombre_colonnes', 'K', 'u_0', 'v_0', 'f', 'I', 'masque', 'R', 't', 'z');
