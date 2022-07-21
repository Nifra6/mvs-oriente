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

% Paramètres du repère monde
delta_x = 20;							% Déplacement du repère monde selon x
delta_y = 30;							% Déplacement du repère monde selon y
deplacement_monde = [delta_x ; delta_y ; 0];

% Paramètres de la calotte de sphère :
rayon_sphere	= 0.9*taille/2;			% Rayon de la sphère
C_x				= taille/2;				% Abscisse du centre de la sphère
C_y				= taille/2;				% Ordonnée du centre de la sphère
alpha			= 0.7;					% Proportion entre les rayons des silhouettes
rayon_calotte 	= alpha*rayon_sphere;	% Rayon de la calotte

%% Préparation de la scène

% Matrices de coordonnées (abscisses orientées vers le bas, ordonnées vers la droite) :
[Y,X]		= meshgrid(1:taille,1:taille);
X			= X - C_x - delta_x;
Y			= Y - C_y - delta_y;
valeurs_X	= (1:taille) - C_x - delta_x;
valeurs_Y	= (1:taille) - C_y - delta_y;

% Éclairage (peut être modifié) :
S = [ 0 ; 0 ; 1 ];

% Angles de rotations
%angles = [-pi/6 0 ; 0 -pi/6 ; pi/12 -pi/17; -pi/19 pi/7];
angles = [-pi/6 0 ; 0 -pi/6];
%angles = [];

% Tirage aléatoire des angles
nombre_predefinis = size(angles,1);
nombre_tirage = 10;
angles = [angles ; zeros(nombre_tirage-nombre_predefinis, 2)];
for k = nombre_predefinis+1:nombre_tirage
	angles(k,:) = 2 * (rand(1,2) - 0.5) * (pi / 6);
end

% Matrices de stockage
nombre_images = size(angles,1) + 1;
I 		 = zeros(taille, taille, nombre_images);		% Ensemble des images
masque 	 = zeros(taille, taille, nombre_images);		% Ensemble des masques
z 		 = zeros(taille, taille, nombre_images);		% Ensemble des profondeurs
N 		 = zeros(taille, taille, 3, nombre_images);	% Ensemble des normales
R        = zeros(3,3,size(angles,1));
R_part_1 = zeros(3,3,size(angles,1));
R_part_2 = zeros(3,3,size(angles,1));
t		 = zeros(3,size(angles,1));

% Rotation entre les deux poses de la caméra, supposée orthogonale :
for k = 1:size(angles,1)
	angle_k = angles(k,:);
	R_part_1(:,:,k) = [1 , 0               , 0               ;...
	  				   0 , cos(angle_k(1)) , -sin(angle_k(1));...
	  				   0 , sin(angle_k(1)) , cos(angle_k(1))];
	R_part_2(:,:,k) = [cos(angle_k(2)) , -sin(angle_k(2)) , 0;...
				       sin(angle_k(2)) , cos(angle_k(2))  , 0;...
				       0               , 0                , 1];
	R(:,:,k) = R_part_1(:,:,k) * R_part_2(:,:,k);
end

%% Calcul de l'image alignée

% Calcul de la première image et de son gradient :
Z_1 = zeros(taille,taille);			% Hauteur de la calotte
dx_I_1 = zeros(taille,taille);		% Dérivée en x de I_1
dy_I_1 = zeros(taille,taille);		% Dérivée en y de I_1
N_1 = zeros(taille,taille,3);		% Normale
for i_1 = 1:taille
	for j_1 = 1:taille
		x_1 = X(i_1,j_1);
		y_1 = Y(i_1,j_1);
		rho_carre_1 = x_1^2 + y_1^2;
		if rho_carre_1 <= rayon_calotte^2
			masque(i_1,j_1,1) = 1;
			z_1 = sqrt(rayon_sphere^2 - rho_carre_1);
			Z_1(i_1,j_1) = z_1;
			z(i_1,j_1,1) = z_1;
			P_1 = [ x_1 ; y_1 ; z_1 ];
			n_1 = P_1 / rayon_sphere;
			ombrage = n_1' * S;
			N_1(i_1,j_1,:) = n_1 / vecnorm(n_1);
			N(i_1,j_1,:,1) = n_1 / vecnorm(n_1);
			if ombrage < 0
				disp('Attention : ombres propres !');
				return;
			end
			I(i_1,j_1,1) = ombrage;
			dx_I_1(i_1,j_1) = -x_1 / (rayon_sphere^2 * ombrage);
			dy_I_1(i_1,j_1) = -y_1 / (rayon_sphere^2 * ombrage);
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
			deplacement 	= deplacement_monde - R(:,:,k) * deplacement_monde; 
			x_k 			= X(i_k,j_k) + deplacement(1);
			y_k 			= Y(i_k,j_k) + deplacement(2);
			rho_carre_k 	= x_k^2 + y_k^2;
			if rho_carre_k <= rayon_sphere^2
				z_k 		= sqrt(rayon_sphere^2 - rho_carre_k);
				z(i_k,j_k,k+1) = z_k;
				P_k 		= [ x_k ; y_k ; z_k ];
				P_1 		= R(:,:,k) * P_k;
				rho_carre_1 = P_1(1)^2 + P_1(2)^2;
				if rho_carre_1 <= rayon_calotte^2
					if P_1(3) < 0
						disp('Attention : parties cachees !');
						return;
					end
					masque(i_k,j_k,k+1) = 1;
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

for k = 1:size(angles,1)
	R(:,:,k) = R_part_2(:,:,k)' * R_part_1(:,:,k);
end



save donnees_calotte_multi;
