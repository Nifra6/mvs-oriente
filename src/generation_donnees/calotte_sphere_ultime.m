%% Clear
clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

%% Paramètres

% Paramètres des images
taille 	= 256;							% Taille des images
I 		= zeros(taille, taille, 3);		% Ensemble des images
masque 	= zeros(taille, taille, 3);		% Ensemble des masques

% Paramètres de la calotte de sphère :
rayon_sphere	= 0.9*taille/2;			% Rayon de la sphère
C_x				= taille/2;				% Abscisse du centre de la sphère
C_y				= taille/2;				% Ordonnée du centre de la sphère
alpha			= 0.7;					% Proportion entre les rayons des silhouettes
rayon_calotte 	= alpha*rayon_sphere;	% Rayon de la calotte

% Paramètres de la caméra
eloignement_camera = rayon_sphere + 100;
f 		= 1;							% Distance focale
u_0 	= C_x;							% Coordonnées du point principal
v_0		= C_y;
K 		= [f 0 u_0 ; 0 f v_0 ; 0 0 1];	% Matrice de calibrage
K_inv 	= [1/f 0 -u_0/f ; 0 1/f -v_0/f ; 0 0 1];


%% Préparation de la scène

% Matrices de coordonnées (abscisses orientées vers le bas, ordonnées vers la droite) :
[Y,X]		= meshgrid(1:taille,1:taille);
X			= X-C_x;
Y			= Y-C_y;
valeurs_X	= (1:taille)-C_x;
valeurs_Y	= (1:taille)-C_y;

% Éclairage (peut être modifié) :
S = [ 0 ; 0 ; 1 ];

% Angles de rotations
%angles = [-pi/6 0 ; 0 -pi/6 ; pi/12 -pi/17; -pi/19 pi/7];
%angles = [-pi/6 0 ; 0 -pi/6];

% Tirage aléatoire des angles
nombre_tirage = 10;
angles = zeros(nombre_tirage,2);
for k = 1: nombre_tirage
	angles(k,:) = 2 * (rand(1,2) - 0.5) * (pi / 6);
end

% Matrices de stockage
I 		= zeros(taille, taille, size(angles,1)+1);		% Ensemble des images
masque 	= zeros(taille, taille, size(angles,1)+1);		% Ensemble des masques
R       = zeros(3,3,size(angles,1)+1);
t		= zeros(3,size(angles,1)+1);

% Translation et rotation de la première caméra
R(:,:,1) = eye(3);
t(:,1) = [0 ; 0 ; eloignement_camera];

% Rotation entre les deux poses de la caméra, supposée orthogonale :
for k = 1:size(angles,1)
	angle_k = angles(k,:);
	R(:,:,k+1) = [cos(angle_k(2))                , 0               ,  sin(angle_k(2));...
			sin(angle_k(1))  * sin(angle_k(2)) , cos(angle_k(1)) , -sin(angle_k(1)) * cos(angle_k(2));...
			-cos(angle_k(1)) * sin(angle_k(2)) , sin(angle_k(1)) ,  cos(angle_k(1)) * cos(angle_k(2))];
	t(:,k+1) = [0 ; 0 ; eloignement_camera];
end

%% Calcul de l'image alignée

% Calcul de la première image et de son gradient :
Z_1 = 1e5*ones(taille,taille);		% Hauteur de la calotte
dx_I_1 = zeros(taille,taille);		% Dérivée en x de I_1
dy_I_1 = zeros(taille,taille);		% Dérivée en y de I_1
N_1 = zeros(taille,taille,3);		% Normale
for i_1 = 1:taille
	for j_1 = 1:taille
		x_1 = X(i_1,j_1);
		y_1 = Y(i_1,j_1);
		rho_carre_1 = x_1^2+y_1^2;
		if rho_carre_1 <= rayon_calotte^2
			masque(i_1,j_1,1) = 1;
			z_1 = sqrt(rayon_sphere^2-rho_carre_1);
			Z_1(i_1,j_1) = eloignement_camera - z_1;
			P_1 = [ x_1 ; y_1 ; z_1 ];
			n_1 = P_1/rayon_sphere;
			ombrage = n_1'*S;
			N_1(i_1,j_1,:) = n_1 / norm(n_1);
			%temp = N_1(i_1,j_1,2);
			%N_1(i_1,j_1,2) = N_1(i_1,j_1,1);
			%N_1(i_1,j_1,1) = N_1(i_1,j_1,2);
			N_1(i_1,j_1,3) = -N_1(i_1,j_1,3);
			if ombrage < 0
				disp('Attention : ombres propres !');
				return;
			end
			I(i_1,j_1,1) = ombrage;
			dx_I_1(i_1,j_1) = -x_1/(rayon_sphere^2*ombrage);
			dy_I_1(i_1,j_1) = -y_1/(rayon_sphere^2*ombrage);
		end
	end
end

% Affichage du relief de la calotte :
figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
plot3(X,Y,Z_1,'k.');
xlabel('$x$','Interpreter','Latex','FontSize',30);
ylabel('$y$','Interpreter','Latex','FontSize',30);
zlabel('$z$','Interpreter','Latex','FontSize',30);
axis equal;
rotate3d;

% Affichage de la première image :
figure('Name','Image 1','Position',[0.33*L,0,0.33*L,0.5*H]);
imagesc(valeurs_Y,valeurs_X,I(:,:,1));
colormap gray;
set(gca,'FontSize',20);
xlabel('$y$','Interpreter','Latex','FontSize',30);
ylabel('$x$','Interpreter','Latex','FontSize',30);
axis equal;

%% Calcul des autres images

dx_I = dx_I_1;
dy_I = dy_I_1;

% Calcul de la k^ème image et de son gradient :
for k = 1:size(angles,1)
	dx_I_k 	= zeros(taille,taille);		% Dérivée en x de I_k
	dy_I_k 	= zeros(taille,taille);		% Dérivée en y de I_k
	angle = angles(k,:);
	for i_k = 1:taille
		for j_k = 1:taille
			x_k 			= X(i_k,j_k);
			y_k 			= Y(i_k,j_k);
			rho_carre_k 	= x_k^2 + y_k^2;
			if rho_carre_k <= rayon_sphere^2
				z_k 		= sqrt(rayon_sphere^2 - rho_carre_k);
				P_k 		= [ x_k ; y_k ; z_k ];
				P_1 		= R(:,:,k+1) * P_k;
				rho_carre_1 = P_1(1)^2 + P_1(2)^2;
				if rho_carre_1 <= rayon_calotte^2
					if P_1(3) < 0
						disp('Attention : parties cachees !');
						return;
					end
					masque(i_k,j_k,2) 	= 1;
					n_1 				= P_1 / rayon_sphere;
					ombrage 			= n_1' * S;
					if ombrage < 0
						disp('Attention : ombres propres !');
						return;
					end
					racine 			= sqrt(1 - (x_k^2+y_k^2) / rayon_sphere^2);
					I(i_k,j_k,k+1) 	= - cos(angle(1)) * sin(angle(2)) * x_k / rayon_sphere ...
						+ sin(angle(1)) * y_k / rayon_sphere ...
						+ cos(angle(1)) * cos(angle(2)) * racine;
					dx_I_k(i_k,j_k) = - cos(angle(1)) * sin(angle(2)) / rayon_sphere ...
						- cos(angle(1)) * cos(angle(2)) * x_k / (rayon_sphere^2 * racine);
					dy_I_k(i_k,j_k) = sin(angle(1)) / rayon_sphere ...
						- cos(angle(1)) * cos(angle(2)) * y_k / (rayon_sphere^2 * racine);
				end
			end
		end
	end
	dx_I(:,:,k+1) = dx_I_k;
	dy_I(:,:,k+1) = dy_I_k;

	% Affichage de la k^ème image :
	if (k < 3)
		figure('Name','Image 2','Position',[0.66*L,0,0.33*L,0.5*H]);
		imagesc(valeurs_Y,valeurs_X,I(:,:,k+1));
		colormap gray;
		set(gca,'FontSize',20);
		xlabel('$y$','Interpreter','Latex','FontSize',30);
		ylabel('$x$','Interpreter','Latex','FontSize',30);
		axis equal;
	end

end


%% Préparation des données pour usage
for k = 1:size(angles,1)+1
	I(:,:,k) = transpose(I(:,:,k));
	masque(:,:,k) = transpose(masque(:,:,k));
end
for k = 1:3
	N_1(:,:,k) = transpose(N_1(:,:,k));
end
Z_1 = Z_1';

save ../../data/donnees_calotte;
