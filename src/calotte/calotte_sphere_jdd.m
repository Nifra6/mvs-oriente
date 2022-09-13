clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

load pose.mat;

t = 200;				% Taille des images
t_sur_2 = t/2;

% Matrices de coordonnées (abscisses orientées vers la droite, ordonnées vers le bas) :
[X_1,Y_1] = meshgrid(1:t,1:t);
X_1 = X_1-t_sur_2;			% On choisit l'origine du repère au centre de l'image
Y_1 = Y_1-t_sur_2;
valeurs_X_1 = (1:t)-t_sur_2;		% Listes de valeurs utiles pour les affichages
valeurs_Y_1 = (1:t)-t_sur_2;

% Paramètres de la calotte de sphère :
R = 0.9*t_sur_2;			% Rayon de la sphère
O_1 = [ 0 ; 0 ; t ];		% Coordonnées du centre de la sphère
alpha = 0.7;				% Proportion entre les rayons des silhouettes
r = alpha*R;				% Rayon de la calotte

% Éclairage (peut être modifié) :
S = [ 0 ; 0 ; -1 ];			% S_z < 0 car ce vecteur vise vers l'appareil photo

% Calcul de la première image et de son gradient :
masque_1 = zeros(t,t);			% Masque de la calotte
Z_1 = zeros(t,t);			% Fonction de profondeur de la calotte
I_1 = zeros(t,t);			% Image de la calotte
N_1 = zeros(t,t,3);
dx_I_1 = zeros(t,t);			% Dérivée en x de I_1
dy_I_1 = zeros(t,t);			% Dérivée en y de I_1
for i_1 = 1:t
	for j_1 = 1:t
		x_1 = X_1(i_1,j_1);
		y_1 = Y_1(i_1,j_1);
		rho_1_carre = (x_1-O_1(1))^2+(y_1-O_1(2))^2;	% Distance au carré à la projection du centre
		if rho_1_carre <= r^2
			masque_1(i_1,j_1) = 1;
			z_1 = O_1(3)-sqrt(R^2-rho_1_carre);
			Z_1(i_1,j_1) = z_1;
			P_1 = [ x_1 ; y_1 ; z_1 ];
			n_1 = (P_1-O_1)/R;
			N_1(i_1,j_1,:) = n_1;
			ombrage = n_1'*S;
			if ombrage < 0
				disp('Attention : ombres propres !');
				return;
			end
			I_1(i_1,j_1) = ombrage;
			dx_I_1(i_1,j_1) = (P_1(1)-O_1(1))/(P_1(3)-O_1(3))/R;
			dy_I_1(i_1,j_1) = (P_1(2)-O_1(2))/(P_1(3)-O_1(3))/R;
		end
	end
end

% Affichage du relief de la calotte :
figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
plot3(X_1,Y_1,Z_1,'k.');
xlabel('$x$','Interpreter','Latex','FontSize',30);
ylabel('$y$','Interpreter','Latex','FontSize',30);
zlabel('$z$','Interpreter','Latex','FontSize',30);
axis equal;
rotate3d;

% Affichage de la deuxième image :
figure('Name','Image 1','Position',[0.33*L,0,0.33*L,0.5*H]);
imagesc(valeurs_X_1,valeurs_Y_1,I_1);
colormap gray;
set(gca,'FontSize',20);
xlabel('$x$','Interpreter','Latex','FontSize',30);
ylabel('$y$','Interpreter','Latex','FontSize',30);
axis equal;




% Matrices de coordonnées (abscisses orientées vers la droite, ordonnées vers le bas) :
[X_2,Y_2] = meshgrid(1:t,1:t);
X_2 = X_2-t_sur_2;			% On choisit l'origine du repère au centre de l'image
Y_2 = Y_2-t_sur_2;
valeurs_X_2 = (1:t)-t_sur_2;		% Listes de valeurs utiles pour les affichages
valeurs_Y_2 = (1:t)-t_sur_2;

% Rotation et translation entre les deux poses de la caméra, supposée orthogonale :
theta = pi/6;
c_1 = cos(theta);
s_1 = sin(theta);
Rotation_1 = [ c_1 0 -s_1 ; 0 1 0 ; s_1 0 c_1 ];
psi = pi/8;
c_2 = cos(psi);
s_2 = sin(psi);
Rotation_2 = [ 1 0 0 ; 0 c_2 -s_2 ; 0 s_2 c_2 ];
Rotation = Rotation_2*Rotation_1;
translation = [ s_1*O_1(3) ; 0 ; (1-c_1)*O_1(3)+50 ];
%translation = [0 ; 0 ; 0];



% Coordonnées du centre de la sphère dans le deuxième repère :
O_2 = Rotation'*(O_1-translation);

% Calcul de la deuxième image et de son gradient :
masque_2 = zeros(t,t);		% Masque de la calotte
I_2 = zeros(t,t);		% Image de la calotte
%dx_I_2 = zeros(t,t);		% Dérivée en x de I_2
%dy_I_2 = zeros(t,t);		% Dérivée en y de I_2
for i_2 = 1:t
	for j_2 = 1:t
		x_2 = X_2(i_2,j_2);
		y_2 = Y_2(i_2,j_2);
		rho_2_carre = (x_2-O_2(1))^2+(y_2-O_2(2))^2;	% Distance au carré à la projection du centre
		if rho_2_carre <= R^2
			z_2 = O_2(3)-sqrt(R^2-rho_2_carre);
			P_2 = [ x_2 ; y_2 ; z_2 ];
			P_1 = Rotation*P_2+translation;
			rho_1_carre = (P_1(1)-O_1(1))^2+(P_1(2)-O_1(2))^2;
			if rho_1_carre <= r^2 & P_1(3) <= O_1(3)
				masque_2(i_2,j_2) = 1;
				n_1 = (P_1-O_1)/R;
				ombrage = n_1'*S;
				I_2(i_2,j_2) = ombrage;
%				dx_I_2(i_2,j_2) = (-s+c*(P_2(1)-O_2(1))/(P_2(3)-O_2(3)))/R;
%				dy_I_2(i_2,j_2) = c/R*(P_2(2)-O_2(2))/(P_2(3)-O_2(3));
			end
		end
	end
end

% Affichage de la deuxième image :
figure('Name','Image 1','Position',[0.66*L,0,0.33*L,0.5*H]);
imagesc(valeurs_X_2,valeurs_Y_2,I_2);
colormap gray;
set(gca,'FontSize',20);
xlabel('$x$','Interpreter','Latex','FontSize',30);
ylabel('$y$','Interpreter','Latex','FontSize',30);
axis equal;

save donnees_calotte;
