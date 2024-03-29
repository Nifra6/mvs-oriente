%% Clear
clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

%% Imports
addpath(genpath('../toolbox/'));

%% Données
load ../../data/data_bunny_ortho.mat;
% Les profondeurs
Z_1 = z(:,:,1);
% Les images
I_1 = Im(:,:,1);
I_2 = Im(:,:,2);
% Les masques des images
masque_1 = mask(:,:,1);
masque_2 = mask(:,:,2);
% La pose
R_2 = R(:,:,2) * R(:,:,1)';
t_2 = t(:,2) - R_2 * t(:,1);
% Le gradient de l'image 2
[dy_I_1, dx_I_1] = gradient(I_1);
[dy_I_2, dx_I_2] = gradient(I_2);
% Ratio pixel distance
u_0 = size(I_1,1)/2;
v_0 = size(I_1,2)/2;
pixelSize = 1.5/540;
% Coordonnées
u_1_coord = pixelSize*(1-u_0):pixelSize:pixelSize*(size(I_1,2)-u_0);
v_1_coord = pixelSize*(1-v_0):pixelSize:pixelSize*(size(I_1,1)-v_0);
[Y,X] = meshgrid(v_1_coord,u_1_coord);

%% Paramètres
valeurs_z 		= 1:0.04:2.6;	% Les valeurs de profondeurs utilisées
range			= 3;		% Voisinage à prendre en compte
affichage_log	= 0;		% Affichage d'informations diverses
interpolation 	= 'nearest';

%% Variables utiles
[i_1_liste, j_1_liste] = find(masque_1);
nb_pixels_utilises = size(i_1_liste,1);
nb_profondeurs_testees = size(valeurs_z,2);
scores = 10 * ones(nb_pixels_utilises, nb_profondeurs_testees);

%% Algorithme

fprintf("\n");
tic;
% Sélection d'un pixel
for indice_pixel = 1:nb_pixels_utilises
	i_1 		= i_1_liste(indice_pixel);
	j_1 		= j_1_liste(indice_pixel);

	% Affichage de la progression des calculs
	if mod(indice_pixel, 100) == 0
		fprintf('\r');
		fprintf("Progression : %d / %d",indice_pixel,nb_pixels_utilises);
	end

	% Sélection de la profondeur
	for indice_z = 1:nb_profondeurs_testees
		z = valeurs_z(indice_z);

		% Changements de repère
		P_1	= [pixelSize*(j_1 - u_0); pixelSize*(i_1 - v_0); z];
		P_2 = R_2 * P_1 + t_2;
		i_2 = round(P_2(2)/pixelSize + u_0);
		j_2 = round(P_2(1)/pixelSize + v_0);

		% Vérification si pixel hors image
		condition_image = i_2 > 0 & i_2 <= size(masque_2,1) & j_2 > 0 & j_2 <= size(masque_2,2);

		% Si le point reprojeté tombe sur le masque de la deuxième image
		if condition_image & masque_2(i_2,j_2)

			% Calcul de la transformation géométrique
			u_1_decales = pixelSize*(j_1-u_0-range):pixelSize:pixelSize*(j_1-u_0+range);
			v_1_decales = pixelSize*(i_1-v_0-range):pixelSize:pixelSize*(i_1-v_0+range);
			[v_1_decales, u_1_decales] = meshgrid(v_1_decales,u_1_decales);

			% Reprojection du voisinage
			P_1_voisinage = [u_1_decales(:)' ; v_1_decales(:)' ; z * ones(1,size(u_1_decales,1)^2)];
			P_2_voisinage = R_2 * P_1_voisinage + t_2;
			i_2_voisinage = round(P_2_voisinage(2,:)/pixelSize + v_0);
			j_2_voisinage = round(P_2_voisinage(1,:)/pixelSize + u_0);

			% Récupération des niveaux de gris dans l'image 2 du voisinage	
			I_2_voisinage = reshape(interp2(I_2, j_2_voisinage(:), i_2_voisinage(:),interpolation),2*range+1,2*range+1)';
			scores(indice_pixel, indice_z) = sum((I_1(i_1-range:i_1+range,j_1-range:j_1+range) - I_2_voisinage).^2,'all');

		end
	end
end
toc;

%% Résultats
% Sélection des profondeurs avec le score minimal
[A, indices_min] = min(scores, [], 2);
z_in = transpose(valeurs_z(indices_min));
z = zeros(size(I_1,1), size(I_1,2));
z(find(masque_1)) = z_in;

% Affichage
%figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
%plot3(X,Y,-z,'k.');
%xlabel('$x$','Interpreter','Latex','FontSize',30);
%ylabel('$y$','Interpreter','Latex','FontSize',30);
%zlabel('$z$','Interpreter','Latex','FontSize',30);
%axis equal;
%rotate3d;

% Affichage du resultat :
XYZ_test = shapeFromDmOrtho(z/pixelSize, masque_1);

figure('Name', 'z from estimated (p,q)')
surfl(XYZ_test(:,:,1),XYZ_test(:,:,2),-XYZ_test(:,:,3),[0 90])
shading flat
colormap gray
axis ij
axis tight
axis off


%% Fonctions annexes

function [i_limite, j_limite] = limites_voisinage(i_voisinage, j_voisinage)
	i_limite = i_voisinage(1,:)';
	i_limite = [i_limite ; i_voisinage(2:end-1,end)];
	i_limite = [i_limite ; wrev(i_voisinage(end,:))'];
	i_limite = [i_limite ; wrev(i_voisinage(2:end-1,1))];
	j_limite = j_voisinage(1,:)';
	j_limite = [j_limite ; j_voisinage(2:end-1,end)];
	j_limite = [j_limite ; wrev(j_voisinage(end,:))'];
	j_limite = [j_limite ; wrev(j_voisinage(2:end-1,1))];
end
