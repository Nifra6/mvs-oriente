%% Clear
clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

%% Données
load ../../data/perspectif/simulateur_formate.mat;
% Indices des images
indice_premiere_image = 2;
indice_deuxieme_image = 4;
% Les profondeurs
Z_1 = z(:,:,indice_premiere_image);
% Les images
I_1 = I(:,:,indice_premiere_image);
I_2 = I(:,:,indice_deuxieme_image);
% TODO Normales
%N_1 = N(:,:,:,indice_premiere_image);
% Les masques
masque_1 = masque(:,:,indice_premiere_image);
masque_2 = masque(:,:,indice_deuxieme_image);
% La pose
R_1_2 = R(:,:,indice_deuxieme_image) * R(:,:,indice_premiere_image)';
t_1_2 = t(:,indice_deuxieme_image) - R_1_2 * t(:,indice_premiere_image);
K_inv = inv(K);
% Le gradient de l'image 2
[dx_I_1, dy_I_1] = gradient(I_1);
[dx_I_2, dy_I_2] = gradient(I_2);


figure; hold on;
rotate3d;
axis equal;
plotCamera('Orientation',R(:,:,indice_premiere_image)','Location',-R(:,:,indice_premiere_image)'*t(:,indice_premiere_image),'Size',0.2);
plotCamera('Orientation',R(:,:,indice_deuxieme_image)','Location',-R(:,:,indice_deuxieme_image)'*t(:,indice_deuxieme_image),'Size',0.2);
box on;


%% Paramètres
rayon_voisinage		= 1;			% Voisinage carré à prendre en compte
affichage_log		= 1;			% Affichage d'informations diverses
interpolation		= 'nearest';	% Type d'interpolation
seuil_denominateur	= 0;			% Seuil pour accepter la division
offset 				= 0.5;

%% Algorithme
% Sélection d'un pixel
figure;
title("Cliquez sur le pixel souhaité.")
imshow(I_1);
P			= drawpoint;
pos 		= P.Position;
i_1 		= round(pos(2));
j_1 		= round(pos(1));
%i_1 = 292;
%j_1 = 191;
grad_I_1	= [dx_I_1(i_1,j_1); dy_I_1(i_1,j_1)];


f2 = figure;
imshow(I_2);
hold on;

% Récupération de la profondeur
nb_z = 100;
liste_z = linspace(0.01,5,nb_z);
for ind_z = 1:nb_z
	z = liste_z(ind_z);

	% Changements de repère
	p_1 = [j_1 - offset ; i_1 - offset ; 1];
	P_1 = z * K_inv * p_1
	P_2 = R_1_2 * P_1 + t_1_2
	p_2 = (K * P_2) / P_2(3);
	u_2 = p_2(1) + offset;
	v_2 = p_2(2) + offset;
	i_2 = v_2;
	j_2 = u_2;
	if (affichage_log)
		disp("===== Changement de repère")
		[i_1 j_1]
		[i_2 j_2]
	end

	plot(j_2,i_2,'+b');

end

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
