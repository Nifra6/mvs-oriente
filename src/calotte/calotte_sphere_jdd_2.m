clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

t = 200;				% Taille des images
t_sur_2 = t/2;

% Paramètres de la calotte de sphère :
R = 0.9*t_sur_2;			% Rayon de la sphère
alpha = 0.7;				% Proportion entre les rayons des silhouettes
r = alpha*R;				% Rayon de la calotte

% Coordonnées du centre de la sphère dans le repère monde :
C_monde = [ 0 ; 0 ; 0 ];

% Éclairage (peut être modifié) :
S_monde = [ 0 ; 0 ; 1 ];		% S_z > 0 car la source est située à l'infini vers les z > 0

% Calcul des N images, avec une caméra supposée orthogonale :
N = 2;
liste_theta_y = [ pi ; pi ; pi ];
liste_theta_x = [ 0 ; pi/6 ; pi/12];
liste_translation_x = [ 0 ; 0 ; 0.5*R ; 0.2*R ];
liste_translation_y = [ 0 ; -0.5*R ; 0.5*R ; 0.5*R ];
liste_translation_z = [ 2*R ; 2*R ; 2*R ; R ];



load pose.mat;


I_k = zeros(t,t,N);
masque_k = zeros(t,t,N);
Z_k = zeros(t,t,N);
R_k = zeros(3,3,N);
t_k = zeros(3,N);
N_1 = zeros(t,t,3);
for k = 1:N

	% Position du centre de la caméra, relativement au repère monde :
	translation_k = [ liste_translation_x(k) ; liste_translation_y(k) ; liste_translation_z(k) ];
	switch k
		case 1
			translation_k = poseMat1(:,4);
		case 2
			translation_k = poseMat2(:,4);
	end
	t_k(:,k) = translation_k;

	% Rotation de la caméra, relativement au repère monde :
	theta_y_k = liste_theta_y(k);
	c_y_k = cos(theta_y_k);
	s_y_k = sin(theta_y_k);
	Rotation_y_k = [ c_y_k 0 -s_y_k ; 0 1 0 ; s_y_k 0 c_y_k ];
	theta_x_k = liste_theta_x(k);
	c_x_k = cos(theta_x_k);
	s_x_k = sin(theta_x_k);
	Rotation_x_k = [ 1 0 0 ; 0 c_x_k -s_x_k ; 0 s_x_k c_x_k ];
	Rotation_k = Rotation_x_k*Rotation_y_k;
	switch k
		case 1
			Rotation_k = poseMat1(:,1:3);
		case 2
			Rotation_k = poseMat2(:,1:3);
	end
	R_k(:,:,k) = Rotation_k;

	% Coordonnées du centre de la sphère dans le repère caméra :
	C_k = Rotation_k'*(C_monde-translation_k);

	% Matrices de coordonnées (abscisses orientées vers la droite, ordonnées vers le bas) :
	[X_k,Y_k] = meshgrid(1:t,1:t);
	X_k = X_k-t_sur_2;			% Origine du repère au centre de l'image
	Y_k = Y_k-t_sur_2;
	valeurs_X_k = (1:t)-t_sur_2;		% Listes de valeurs utiles pour les affichages
	valeurs_Y_k = (1:t)-t_sur_2;
	if (k == 1)
		X_1 = X_k;
		Y_1 = Y_k;
	end

	% Calcul de l'image :
	for i = 1:t
		for j = 1:t
			x_k = X_k(i,j);
			y_k = Y_k(i,j);
			rho_k_carre = (x_k-C_k(1))^2+(y_k-C_k(2))^2;	% Distance au carré à la projection du centre
			if rho_k_carre <= R^2
				z_k = C_k(3)-sqrt(R^2-rho_k_carre);
				Z_k(i,j,k) = z_k;
				P_k = [ x_k ; y_k ; z_k ];
				P_monde = Rotation_k*P_k+translation_k;
				rho_monde_carre = (P_monde(1)-C_monde(1))^2+(P_monde(2)-C_monde(2))^2;
				if rho_monde_carre <= r^2 & P_k(3) <= C_k(3)
					masque_k(i,j,k) = 1;
					n_monde = (P_monde-C_monde)/R;
					if (k == 1)
						N_1(i,j,:) = n_monde;
					end
					ombrage = n_monde'*S_monde;
					I_k(i,j,k) = ombrage;
				end
			end
		end
	end

	% Affichage de l'image :
	figure('Name','Image k','Position',[(k-1)*0.33*L,0,0.33*L,0.5*H]);
	imagesc(valeurs_X_k,valeurs_Y_k,I_k(:,:,k));
	colormap gray;
	set(gca,'FontSize',20);
	xlabel('$x$','Interpreter','Latex','FontSize',30);
	ylabel('$y$','Interpreter','Latex','FontSize',30);
	axis equal;

end

save donnees_calotte;
