%% Trucs de Matlab
% Clear
clear;
close all;
% Paramètres d'affichage
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);
% Imports de fonctions utiles
addpath(genpath('../toolbox/'));

%% Paramètres
nombre_z		= 11;				% Valeurs de profondeurs testées
interpolation 	= 'linear';			% Type d'interpolation
estimateur		= 'MSE';			% Estimateur utilisé pour l'évaluation des erreurs
affichage 		= 'Iteration';		% Type d'affichage de la progression
affichage_debug = 0;				% Affichage d'informations diverses
rayon_voisinage = 1;				% Rayon du voisinage carré à prendre en compte
filtrage 		= 1;				% Utilisation d'un filtrage gaussien
sigma_filtre 	= 5;			    % Écart type du filtre gaussien

%% Données
% Fichier des données
load ../../data/simulateur_formate.mat;
% Taille des images
nombre_pixels = nombre_lignes * nombre_colonnes;
taille_patch  = (2*rayon_voisinage + 1)^2;	% Nombre de pixels dans un patch
% Les profondeurs
Z_1 = z(:,:,1);
valeurs_z = linspace(min(Z_1,[],'all'),max(Z_1,[],'all'),nombre_z);	% Valeurs de profondeurs testées
% Les normales (pour debug)
N_1 = N(:,:,:,1);
% Les poses relatives
R_1_k = zeros(3,3,nombre_images-1);
t_1_k = zeros(3,nombre_images-1);
for k = 1:nombre_images-1
	R_1_k(:,:,k) = R(:,:,k+1) * R(:,:,1)';
	t_1_k(:,k) = t(:,k+1) - R_1_k(:,:,k) * t(:,1);
end
% TODO un truc plus propre pour ce qui suit
% Modifications du masque (pour correspondre aux patchs utilisés)
for k = 1:1
	masque(1:rayon_voisinage,:,k) = 0;
	masque(end-rayon_voisinage:end,:,k) = 0;
	masque(:,1:rayon_voisinage,k) = 0;
	masque(:,end-rayon_voisinage:end,k) = 0;
end
% Filtrage des pixels considérés par le masque
[i_k, j_k]  = find(masque(:,:,1));
ind_1		= sub2ind([nombre_lignes nombre_colonnes], i_k, j_k);
nombre_pixels_etudies = size(ind_1,1);
P_k 		= zeros(3,nombre_pixels_etudies,nombre_images);
P_k(:,:,1) 	= [(j_k - u_0) / facteur_k, (i_k - v_0) / facteur_k, zeros(length(i_k), 1)].';


%% Calcul du filtre
if (filtrage)
	% Analytique
	rayon_masque = sigma_filtre * 4;
	taille_masque = (2*rayon_masque+1)^2;
	[x,y] = meshgrid(-rayon_masque:rayon_masque,-rayon_masque:rayon_masque);
	u_x = 0; u_y = 0;
	filtre = 1./(2*pi*sigma_filtre^2) .* exp(((x-u_x).^2+(y-u_y).^2)./(2*sigma_filtre^2));
	filtre = filtre / sum(filtre(:));
	dx_filtre = -(x-u_x)./(2*pi*sigma_filtre^4) .* exp(((x-u_x).^2+(y-u_y).^2)./(2*sigma_filtre^2));
	dy_filtre = -(y-u_y)./(2*pi*sigma_filtre^4) .* exp(((x-u_x).^2+(y-u_y).^2)./(2*sigma_filtre^2));
	dx_filtre = dx_filtre;
	dy_filtre = dy_filtre;

	% Fonctions Matlab
	%u_x = 0; u_y = 0;
	%cote_masque = ceil(sqrt(8*sigma_filtre));
	%filtre = fspecial('gauss',cote_masque,sigma_filtre);
	%filtre = filtre / sum(filtre(:));
	%dx_conv = [0 0 0 ; 1 0 -1 ; 0 0 0];
	%dy_conv = [0 1 0 ; 0 0 0 ; 0 -1 0];
	%dx_filtre = conv2(filtre,dx_conv);
	%dy_filtre = conv2(filtre,dy_conv);
	%dx_filtre = dx_filtre / sum(dx_filtre(:));
	%dy_filtre = dy_filtre / sum(dy_filtre(:));

	% Filtrage de l'image
	I_filtre = zeros(size(I));
	for k = 1:nombre_images
		I_filtre(:,:,k) = conv2(I(:,:,k),filtre,'same');
	end
else
	I_filtre = I;
end

%% Calcul des gradients
dx_I_k = zeros(size(I));
dy_I_k = zeros(size(I));
for k = 1:nombre_images
	if (filtrage)
		% Gradient filtré
		dx_I = conv2(I(:,:,k),dx_filtre,'same');
		dy_I = conv2(I(:,:,k),dy_filtre,'same');
	else
		% Gradient non filtré
		[dx_I, dy_I] = gradient(I(:,:,k));
		%[dx_I, dy_I] = gradient_correct(I(:,:,k),masque(:,:,k),1);
	end
	% Sauvegarde des gradients
	dx_I_k(:,:,k) = dx_I;
	dy_I_k(:,:,k) = dy_I;
end
dx_I_1 = dx_I_k(:,:,1);
dy_I_1 = dy_I_k(:,:,1);
grad_I_1 	= [ dx_I_1(ind_1) , dy_I_1(ind_1) ].';
grad_I_x	= [dx_I_1(ind_1)'];
grad_I_y    = [dy_I_1(ind_1)'];

%% Construction du voisinage
voisinage_ligne = -rayon_voisinage*nombre_lignes:nombre_lignes:rayon_voisinage*nombre_lignes;
voisinage_colonne = -rayon_voisinage:rayon_voisinage;
grille_voisinage = voisinage_ligne + voisinage_colonne';
grille_voisinage = grille_voisinage';

%% Mise en forme des normales
normale_theorique = [N_1(ind_1)' ; N_1(ind_1 + nombre_pixels)' ; N_1(ind_1 + 2*nombre_pixels)'];

%% Boucle de reconstruction
erreurs	= 10*ones(nombre_pixels_etudies, nombre_z);
n_estimes = zeros(3, nombre_pixels_etudies, nombre_z);

tic
fprintf("\n")
for indice_z = 1:nombre_z

	% Affichage de la progression des calculs
	switch (affichage)
		case 'Iteration'
			fprintf('\r');
			fprintf("Progression : %d / %d",indice_z,nombre_z);
		case 'Pourcentage'
			if mod(indice_z,round(nombre_z/25)) == 0
				disp("Progression à " + int2str(indice_z/nombre_z*100) + "%");
			end
	end

	% Sélection d'une profondeur
	valeur_z 	= valeurs_z(indice_z);
	P_k(3,:,1) 	= valeur_z;

	% Changements de repère
	for k = 1:nombre_images-1
		P_k(:,:,k+1) = R_1_k(:,:,k) * P_k(:,:,1) + t_1_k(:,k);
		i_k(:,k+1) = (P_k(2,:,k+1) * facteur_k + v_0).';
		j_k(:,k+1) = (P_k(1,:,k+1) * facteur_k + u_0).';
	end

	% Vérification des pixels hors images
	condition_image = ones(nombre_pixels_etudies,nombre_images-1);
	for k = 1:nombre_images-1
		condition_image(:,k) = i_k(:,k+1) > 0.5 & i_k(:,k+1) <= nombre_lignes & j_k(:,k+1) > 0.5 & j_k(:,k+1) <= nombre_colonnes;
	end

	% Calcul des gradients
	for k = 1:nombre_images-1
		i_k(:,k+1) = (ones(nombre_pixels_etudies,1) - condition_image(:,k)) + condition_image(:,k) .* i_k(:,k+1);
		j_k(:,k+1) = (ones(nombre_pixels_etudies,1) - condition_image(:,k)) + condition_image(:,k) .* j_k(:,k+1);
		% P'têt une carabistouille ci-dessus
		grad_I_x(k+1,:) = interp2(dx_I_k(:,:,k+1),j_k(:,k+1),i_k(:,k+1),interpolation)';
		grad_I_y(k+1,:) = interp2(dy_I_k(:,:,k+1),j_k(:,k+1),i_k(:,k+1),interpolation)';
		i_k(:,k+1) = round(i_k(:,k+1));
		j_k(:,k+1) = round(j_k(:,k+1));
		ind(:,k+1) = sub2ind([nombre_lignes nombre_colonnes], i_k(:,k+1), j_k(:,k+1));
	end

	% Calcul des numérateurs et dénominateurs
	denominateur = [];
	numerateur_x = [];
	numerateur_y = [];
	for k = 1:nombre_images-1
		numerateur = [grad_I_x(1,:); grad_I_y(1,:)] - R_1_k(1:2,1:2,k)' * [grad_I_x(k+1,:); grad_I_y(k+1,:)];
		numerateur_x(k,:) = numerateur(1,:);
		numerateur_y(k,:) = numerateur(2,:);
		denominateur(k,:) = R_1_k(1:2,3,k)' * [grad_I_x(k+1,:); grad_I_y(k+1,:)];
	end
	clear numerateur

	% Calcul des coefficients p et q
	p_q = 0;	
	for k = 1:nombre_images-1
		p_q = p_q + denominateur(k,:) .* [numerateur_x(k,:); numerateur_y(k,:)];
	end
	p_q 	= p_q ./ sum(denominateur.^2, 1);
	p_estim = p_q(1, :);
	q_estim = p_q(2, :);
	clear denominateur  numerateur_x numerateur_y;

	% Calcul de la normale
	normale = [p_estim ; q_estim ; -ones(1,nombre_pixels_etudies)] ./ sqrt(p_estim.^2 + q_estim.^2 + ones(1,nombre_pixels_etudies));
	n_estimes(:,:,indice_z) = normale;

	% Calcul du plan considéré
	d_equation_plan = sum(-P_k(:,:,1) .* normale,1);

	% Calcul de la transformation géométrique
	ind_decales = ind_1 + grille_voisinage(:)'; % Création de matrice avec 2 vecteurs
	[i_1_decales, j_1_decales] = ind2sub([nombre_lignes, nombre_colonnes], ind_decales);
	u_1_decales = (j_1_decales-u_0) / facteur_k ;
	v_1_decales = (i_1_decales-v_0) / facteur_k;

	normale_1 = repmat(normale(1,:)',1,taille_patch);
	normale_2 = repmat(normale(2,:)',1,taille_patch);
	normale_3 = repmat(normale(3,:)',1,taille_patch);
	z_1_decales = -(d_equation_plan' + normale_1.*u_1_decales + normale_2.*v_1_decales)./normale_3;

	% Reprojection du voisinage
	i_2_voisinage = zeros(nombre_pixels_etudies, taille_patch, nombre_images-1);
	j_2_voisinage = zeros(nombre_pixels_etudies, taille_patch, nombre_images-1);
	u_1_decales_vec = reshape(u_1_decales',1,nombre_pixels_etudies*taille_patch);
	v_1_decales_vec = reshape(v_1_decales',1,nombre_pixels_etudies*taille_patch);
	z_1_decales_vec = reshape(z_1_decales',1,nombre_pixels_etudies*taille_patch);
	P_1_voisinage = [u_1_decales_vec ; v_1_decales_vec ; z_1_decales_vec];
	for k = 1:nombre_images-1
		P_2_voisinage = R_1_k(:,:,k) * P_1_voisinage + t_1_k(:,k);
		P_2_voisinage_ok = cell2mat(mat2cell(P_2_voisinage,3,repmat(taille_patch,1,nombre_pixels_etudies))');
		i_2_voisinage(:,:,k) = round(P_2_voisinage_ok(2:3:end,:) * facteur_k + v_0);
		j_2_voisinage(:,:,k) = round(P_2_voisinage_ok(1:3:end,:) * facteur_k + u_0);
	end

	% Calcul de l'erreur
	I_1_voisinage = interp2(I_filtre(:,:,1),j_1_decales,i_1_decales,interpolation);
	erreur_k = zeros(nombre_pixels_etudies, nombre_images-1);
	for k = 1:nombre_images-1
		I_k_voisinage = interp2(I_filtre(:,:,k+1),j_2_voisinage(:,:,k),i_2_voisinage(:,:,k),interpolation);
		erreur_k(:,k) = prod(condition_image,2).*sum((I_1_voisinage-I_k_voisinage),2);
	end
	%erreur_k = erreur_k + 10 * (1 - condition_image);
	switch (estimateur)
		case 'MSE'
			erreurs(:,indice_z) = (1 / nombre_images) * sum(erreur_k.^2,2);
		case 'Robuste'
			erreurs(:,indice_z) = (1 / nombre_images) * (1 - exp(-sum(erreur_k.^2,2)/0.2^2));
	end


	% Affichage debug
	ind_debug = 191500;
	if (affichage_debug && round(valeur_z) == round(Z_1(i_k(ind_debug),j_k(ind_debug))))
		[i_k(ind_debug) j_k(ind_debug)]
		Z_1(i_k(ind_debug),j_k(ind_debug))
		i_1_decales(ind_debug,:)
		j_1_decales(ind_debug,:)
		disp("===== Pente et normale")
		[p_estim(ind_debug) q_estim(ind_debug)]
		normale(:,ind_debug)
		N_1(i_k(ind_debug), j_k(ind_debug), :)
		erreur_angulaire(ind_debug)
		-P_k(:,ind_debug,1)
		d_equation_plan(ind_debug)
		z_1_decales(ind_debug,:)
		disp("===== Voisinages sur les images témoins")
		i_2_voisinage(ind_debug,:)
		j_2_voisinage(ind_debug,:)
		I_1_voisinage(ind_debug,:)
		I_k_voisinage(ind_debug,:)
		size(j_1_decales)
		size(j_2_voisinage(:,:,1))
	end


end

fprintf('\n');
toc

%% Résultats
format long;

% Sélections des profondeurs avec l'erreur minimale
erreurs_corrigees = (erreurs ~= 0) .* erreurs + 1 * (erreurs == 0) .* ones(size(erreurs));
%erreurs_corrigees = erreurs;
[~,indices_min] = min(erreurs_corrigees,[],2);
z_in = transpose(valeurs_z(indices_min));
z_estime = nan(nombre_lignes, nombre_colonnes);
z_estime(ind_1) = z_in;

% Sélection des normales
n_totales = zeros(3, nombre_pixels);
n_totales_ind = zeros(3, nombre_pixels_etudies);
for k = 1:3
	n_totales(k,ind_1) = n_estimes(sub2ind(size(n_estimes), k * ones(nombre_pixels_etudies,1), transpose(1:nombre_pixels_etudies), indices_min));
	n_totales_ind(k,:) = n_totales(k,ind_1);
end

% Sélections des erreurs angulaires
angles = zeros(nombre_lignes, nombre_colonnes);
angles(ind_1) = abs((180/pi) * atan2(vecnorm(cross(normale_theorique,n_totales_ind)),dot(normale_theorique,n_totales_ind)));

% Mesures
disp("==============")
disp("Mesure relative de profondeur")
sum(abs(Z_1(ind_1) - z_estime(ind_1)),'all') / nombre_pixels_etudies
ecart_moyen = sum(Z_1(ind_1) - z_estime(ind_1)) / size(z_in,1);
disp("Mesure relative de forme")
sum(abs(Z_1(ind_1) - (z_estime(ind_1) + ecart_moyen)),'all') / size(z_in,1)

% Préparation
X = 1:nombre_colonnes;
X = (X - u_0) / facteur_k;
Y = 1:nombre_lignes;
Y = (Y - v_0) / facteur_k;
[X,Y] = meshgrid(X,Y);

% Affichage
figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
sl = surfl(X,Y,-z_estime,s);
sl.EdgeColor = 'none';
grid off;
colormap gray;
axis equal;

% Affichage des erreurs angulaires
max(angles,[],'all')
erreur_angle_moy = mean(angles(:))
erreur_angle_med = median(angles(:))
%save('../../data/angles_cubic.mat','angles');
angles_norm_max = angles / max(angles,[],'all');
angles_norm_180 = angles / 180;
figure;
imshow(angles_norm_max);

figure;
imshow(angles_norm_180);

save('../../data/normales_simulateur_estimes.mat','n_totales');

