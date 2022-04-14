clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

taille = 256;				% Taille des images

% Paramètres de la calotte de sphère :
rayon_sphere = 0.9*taille/2;		% Rayon de la sphère
C_x = taille/2;				% Abscisse du centre de la sphère
C_y = taille/2;				% Ordonnée du centre de la sphère
alpha = 0.7;				% Proportion entre les rayons des silhouettes
rayon_calotte = alpha*rayon_sphere;	% Rayon de la calotte

% Matrices de coordonnées (abscisses orientées vers le bas, ordonnées vers la droite) :
[Y,X] = meshgrid(1:taille,1:taille);
X = X-C_x;
Y = Y-C_y;
valeurs_X = (1:taille)-C_x;
valeurs_Y = (1:taille)-C_y;

% Éclairage (peut être modifié) :
S = [ 0 ; 0 ; 1 ];

% Rotation entre les deux poses de la caméra, supposée orthogonale :
theta_2 = -pi/6;
cos_theta_2 = cos(theta_2);
sin_theta_2 = sin(theta_2);
R_2 = [ 1 0 0 ; 0 cos_theta_2 -sin_theta_2 ; 0 sin_theta_2 cos_theta_2 ];

% Deuxième rotation entre les deux poses de la caméra, supposée orthogonale :
theta_3 = -pi/6;
cos_theta_3 = cos(theta_3);
sin_theta_3 = sin(theta_3);
R_3 = [ cos_theta_3 0 sin_theta_3 ; 0 1 0 ; -sin_theta_3 0 cos_theta_3 ];

% Calcul de la première image et de son gradient :
masque_1 = zeros(taille,taille);	% Masque de la calotte
Z_1 = zeros(taille,taille);		% Hauteur de la calotte
I_1 = zeros(taille,taille);		% Image de la calotte
dx_I_1 = zeros(taille,taille);		% Dérivée en x de I_1
dy_I_1 = zeros(taille,taille);		% Dérivée en y de I_1
for i_1 = 1:taille
	for j_1 = 1:taille
		x_1 = X(i_1,j_1);
		y_1 = Y(i_1,j_1);
		rho_carre_1 = x_1^2+y_1^2;
		if rho_carre_1 <= rayon_calotte^2
			masque_1(i_1,j_1) = 1;
			z_1 = sqrt(rayon_sphere^2-rho_carre_1);
			Z_1(i_1,j_1) = z_1;
			P_1 = [ x_1 ; y_1 ; z_1 ];
			n_1 = P_1/rayon_sphere;
			ombrage = n_1'*S;
			if ombrage < 0
				disp('Attention : ombres propres !');
				return;
			end
			I_1(i_1,j_1) = ombrage;
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
imagesc(valeurs_Y,valeurs_X,I_1);
colormap gray;
set(gca,'FontSize',20);
xlabel('$y$','Interpreter','Latex','FontSize',30);
ylabel('$x$','Interpreter','Latex','FontSize',30);
axis equal;

% Calcul de la deuxième image et de son gradient :
masque_2 = zeros(taille,taille);	% Masque de la calotte
I_2 = zeros(taille,taille);			% Image de la calotte
dx_I_2 = zeros(taille,taille);		% Dérivée en x de I_2
dy_I_2 = zeros(taille,taille);		% Dérivée en y de I_2
for i_2 = 1:taille
	for j_2 = 1:taille
		x_2 = X(i_2,j_2);
		y_2 = Y(i_2,j_2);
		rho_carre_2 = x_2^2+y_2^2;
		if rho_carre_2 <= rayon_sphere^2
			z_2 = sqrt(rayon_sphere^2-rho_carre_2);
			P_2 = [ x_2 ; y_2 ; z_2 ];
			P_1 = R_2*P_2;
			rho_carre_1 = P_1(1)^2+P_1(2)^2;
			if rho_carre_1 <= rayon_calotte^2
				if P_1(3) < 0
					disp('Attention : parties cachees !');
					return;
				end
				masque_2(i_2,j_2) = 1;
				n_1 = P_1/rayon_sphere;
				ombrage = n_1'*S;
				if ombrage < 0
					disp('Attention : ombres propres !');
					return;
				end
				I_2(i_2,j_2) = ombrage;
				racine = sqrt(1-(x_2^2+y_2^2)/rayon_sphere^2);
				I_2(i_2,j_2) = sin_theta_2*y_2/rayon_sphere+cos_theta_2*racine;
				dx_I_2(i_2,j_2) = -cos_theta_2*x_2/(rayon_sphere^2*racine);
				dy_I_2(i_2,j_2) = sin_theta_2/rayon_sphere-...
							cos_theta_2*y_2/(rayon_sphere^2*racine);
			end
		end
	end
end


% Affichage de la deuxième image :
figure('Name','Image 2','Position',[0.66*L,0,0.33*L,0.5*H]);
imagesc(valeurs_Y,valeurs_X,I_2);
colormap gray;
set(gca,'FontSize',20);
xlabel('$y$','Interpreter','Latex','FontSize',30);
ylabel('$x$','Interpreter','Latex','FontSize',30);
axis equal;


% Calcul de la troisième image et de son gradient :
masque_3 = zeros(taille,taille);	% Masque de la calotte
I_3 = zeros(taille,taille);			% Image de la calotte
dx_I_3 = zeros(taille,taille);		% Dérivée en x de I_3
dy_I_3 = zeros(taille,taille);		% Dérivée en y de I_3
for i_3 = 1:taille
	for j_3 = 1:taille
		x_3 = X(i_3,j_3);
		y_3 = Y(i_3,j_3);
		rho_carre_3 = x_3^2+y_3^2;
		if rho_carre_3 <= rayon_sphere^2
			z_3 = sqrt(rayon_sphere^2-rho_carre_3);
			P_3 = [ x_3 ; y_3 ; z_3 ];
			P_1 = R_3*P_3;
			rho_carre_1 = P_1(1)^2+P_1(2)^2;
			if rho_carre_1 <= rayon_calotte^2
				if P_1(3) < 0
					disp('Attention : parties cachees !');
					return;
				end
				masque_3(i_3,j_3) = 1;
				n_1 = P_1/rayon_sphere;
				ombrage = n_1'*S;
				if ombrage < 0
					disp('Attention : ombres propres !');
					return;
				end
				I_3(i_3,j_3) = ombrage;
				racine = sqrt(1-(x_3^2+y_3^2)/rayon_sphere^2);
				I_3(i_3,j_3) = -sin_theta_3*x_3/rayon_sphere+cos_theta_3*racine;
				dx_I_3(i_3,j_3) = -sin_theta_3/rayon_sphere-...
							cos_theta_3*x_3/(rayon_sphere^2*racine);
				dy_I_3(i_3,j_3) = -cos_theta_3*y_3/(rayon_sphere^2*racine);
			end
		end
	end
end


% Affichage de la troisième image :
figure('Name','Image 3','Position',[0.66*L,0,0.33*L,0.5*H]);
imagesc(valeurs_Y,valeurs_X,I_3);
colormap gray;
set(gca,'FontSize',20);
xlabel('$y$','Interpreter','Latex','FontSize',30);
ylabel('$x$','Interpreter','Latex','FontSize',30);
axis equal;



save donnees_calotte;
