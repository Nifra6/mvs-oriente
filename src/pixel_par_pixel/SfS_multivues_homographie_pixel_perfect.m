%% Clear
clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

%% Données
load ../../data/donnees_calotte;
I_1 = I(:,:,1);
I_2 = I(:,:,2);
masque_1 = masque(:,:,1);
masque_2 = masque(:,:,2);
R_2 = R(:,:,1)';
t_2 = t(:,1);
dx_I_2 = dx_I(:,:,2);
dy_I_2 = dy_I(:,:,2);
%[dy_I_2, dx_I_2] = gradient(I_2);

%% Paramètres
valeurs_z = 60:1:120;	% Les valeurs de profondeurs utilisées
range = 1;				% Voisinage à prendre en compte

%% Algorithme
while (1)
	% Sélection d'un pixel
	figure;
	imshow(I_1);
	P = drawpoint;
	pos = P.Position;
	i_1 = round(pos(2));
	j_1 = round(pos(1));
	grad_I_1 = [dy_I_1(i_1,j_1); dx_I_1(i_1,j_1)];

	% Récupération de la profondeur
	z = Z_1(i_1,j_1);

	% Changements de repère
	P_1 = [i_1 - u_0; j_1 - v_0; z];
	P_2 = R_2 * P_1;
	i_2 = round(P_2(1) + u_0);
	j_2 = round(P_2(2) + v_0);

	% Vérification si pixel hors image
	condition_image = i_2 > 0 & i_2 <= size(masque_2,1) & j_2 > 0 & j_2 <= size(masque_2,2);

	% Si le point reprojeté tombe sur le masque de la deuxième image
	if condition_image & masque_2(i_2,j_2)

		grad_I_2 = [dx_I_2(i_2,j_2); dy_I_2(i_2,j_2)];

		denominateur_pq = R_2(1:2,3)' * grad_I_2;
		numerateur_pq = grad_I_1 - R_2(1:2,1:2)' * grad_I_2;

		% Si pas de division par 0, on continue
		if abs(denominateur_pq) > 0

			% Estimation de la pente
			p_estime = numerateur_pq(1) / denominateur_pq;
			q_estime = numerateur_pq(2) / denominateur_pq;

			% Calcul du plan au pixel considéré
			normale = (1 / (p_estime^2 + q_estime^2 + 1)) * [p_estime ; q_estime ; 1];
			disp("Comparaison des normales")
			N_1(i_1,j_1) - normale 
			d_equation_plan = -P_1' * normale;

			% Calcul de la transformation géométrique
			v_1_decales = i_1-v_0-range:i_1-v_0+range;
			u_1_decales = j_1-u_0-range:j_1-u_0+range;
			[u_1_decales, v_1_decales] = meshgrid(u_1_decales,v_1_decales); %ok
			z_1_decales = -(d_equation_plan + normale(1) * u_1_decales(:) + normale(2) * v_1_decales(:)) / normale(3);
			z;

			% Reprojection du voisinage
			P_1_voisinage = [u_1_decales(:)' ; v_1_decales(:)' ; z_1_decales'];
			P_2_voisinage = R_2' * P_1_voisinage;
			i_2_voisinage = round(P_2_voisinage(2,:) + v_0);
			j_2_voisinage = round(P_2_voisinage(1,:) + u_0);
			i_2;
			j_2;



			%if (indice_pixel == 12475 && z == round(Z_1(i_1,j_1)))
			if (z == round(Z_1(i_1,j_1)))
				z_1_decales
				z
				reshape(i_2_voisinage, 2*range + 1, 2*range + 1) - i_2
				reshape(j_2_voisinage, 2*range + 1, 2*range + 1) - j_2
				figure;
				subplot(2,2,1);
				imshow(I_1(i_1-range:i_1+range,j_1-range:j_1+range));
				subplot(2,2,2);
				i_bruh = reshape(interp2(I_2, i_2_voisinage(:), j_2_voisinage(:),'nearest'),2*range+1,2*range+1);
				i_bruh
				imshow(i_bruh);
				subplot(2,2,3);
				imshow(I_1);
				axis on
				hold on;
				i_1_limite = [i_1 - range ; i_1 - range ; i_1 + range ; i_1 + range];
				j_1_limite = [j_1 - range ; j_1 + range ; j_1 + range ; j_1 - range];
				fill(i_1_limite,j_1_limite,'g')
				plot(i_1,j_1, 'r+', 'MarkerSize', 30, 'LineWidth', 2)
				hold off;
				subplot(2,2,4);
				imshow(I_2)
				axis on
				hold on;
				i_2_voisinage = reshape(i_2_voisinage, 2*range+1, 2*range+1);
				j_2_voisinage = reshape(j_2_voisinage, 2*range+1, 2*range+1);
				i_2_limite = [i_2_voisinage(1,1) ; i_2_voisinage(1,end) ; i_2_voisinage(end,end) ; i_2_voisinage(end,1)];
				j_2_limite = [j_2_voisinage(1,1) ; j_2_voisinage(1,end) ; j_2_voisinage(end,end) ; j_2_voisinage(end,1)];
				fill(i_2_limite,j_2_limite,'g')
				plot(i_2,j_2, 'r+', 'MarkerSize', 30, 'LineWidth', 2)
				hold off;
			end
			%H_ij = homographie * [i_1 ; j_1 ; 1];
			%G_ij = H_ij(1:2) / H_ij(3);


			residu_1 = I_1(i_1,j_1) - 1 / sqrt(p_estime^2 + q_estime^2 + 1);
			%residu_1 = I_1(i_1,j_1) - interp2(I_2,G_ij(2),G_ij(1),'nearest');
			%residu_2 = I_1(i_1,j_1) - I_2(i_2,j_2);
			residu_2 = 1;

			%scores(indice_pixel, indice_z) = residu_1^2 + lambda * residu_2^2;
		end
	end

	disp("Appuyez sur une touche pour recommencer.")
	pause;
	close all;
end


%% Résultats
% Sélection des profondeurs avec le score minimal
[A, indices_min] = min(scores, [], 2);
z_in = transpose(valeurs_z(indices_min));
z = zeros(256, 256);
z(find(masque_1)) = z_in;

% Affichage
%figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
%plot3(X,Y,z,'k.');
%%plot3(X,Y,Z_1,'k.');
%xlabel('$x$','Interpreter','Latex','FontSize',30);
%ylabel('$y$','Interpreter','Latex','FontSize',30);
%zlabel('$z$','Interpreter','Latex','FontSize',30);
%axis equal;
%rotate3d;


% Estimation de z avec p et q
p_estimes = liste_p_estimes(indices_min);
q_estimes = liste_p_estimes(indices_min);
p = zeros(taille,taille);
q = zeros(taille,taille);
p(find(masque_1)) = p_estimes;
q(find(masque_1)) = q_estimes;
z_estim = integration_SCS(p,q);



%% Fonction annexe

function images_decalees = decaler(image)
	[nombre_lignes, nombre_colonnes, nombres_images] = size(image);
	images_decalees = zeros(nombre_lignes, nombre_colonnes, nombres_images, 9);
	for k = 1:nombres_images
		images_decalees(:,:,k,:) = decalage(image(:,:,k));
	end
end
