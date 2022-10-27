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
R_1_m = R(:,:,indice_premiere_image)'; 
R_2_m = R(:,:,indice_deuxieme_image)';
t_1_m = -R_1_m * t(:,indice_premiere_image);
t_2_m = -R_2_m * t(:,indice_deuxieme_image);
%pose1 = rigid3d(R_1_m,t_1_m');
%pose2 = rigid3d(R_2_m,t_2_m');
cam1 = plotCamera('Orientation', R_1_m', 'Location', t_1_m, 'Opacity', 0, 'Size', 0.2, 'Color', [0 0 1], 'Label', 'Référence', 'AxesVisible', 1);
cam2 = plotCamera('Orientation', R_2_m', 'Location', t_2_m, 'Opacity', 0, 'Size', 0.2, 'Color', [1 0 0], 'Label', 'Témoin', 'AxesVisible', 1);
%plotCamera('Orientation',R(:,:,indice_premiere_image)','Location',R(:,:,indice_premiere_image)'*t(:,indice_premiere_image),'Size',0.2);
%plotCamera('Orientation',R(:,:,indice_deuxieme_image)','Location',R(:,:,indice_deuxieme_image)'*t(:,indice_deuxieme_image),'Size',0.2);
box on;
xlabel('x');
ylabel('y');
zlabel('z');


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
title('Image de référence');
P			= drawpoint;
pos 		= P.Position;
i_1 		= round(pos(2));
j_1 		= round(pos(1));
%i_1 = 128;
%j_1 = 250;
grad_I_1	= [dx_I_1(i_1,j_1); dy_I_1(i_1,j_1)];


f2 = figure;
imshow(I_2);
title('Image témoin');
hold on;

% Récupération de la profondeur
nb_z = 100;
liste_z = linspace(1,2.5,nb_z);
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
