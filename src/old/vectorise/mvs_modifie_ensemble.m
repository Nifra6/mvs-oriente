clear;
close all;

%% Paramètres
interpolation 	= 'linear';			% Type d'interpolation utilisée
estimateur		= 'MSE';			% Estimateur utilisé pour l'évaluation des erreurs photométriques
affichage 		= 'Iteration';		% Type d'affichage de la progression de l'algorithme
offset 			= 0.5;				% Décalage spatial entre les indices des pixels et leur coordonnées
surface 		= 'calotte_sim_flou2';
nb_vues			= 9;
rayon_voisinage = 0;


%% Données
% Chargement des fonctions utiles
addpath(genpath("../toolbox/"));
% Chargement des données
path = "../../data/";
nom_fichier = "simulateur_" + surface + "_formate.mat";
load(path+nom_fichier);
% Nombres d'images et de pixels considérés
nb_pixels = nb_lignes * nb_colonnes;
nb_images = nb_vues;
taille_patch  = (2*rayon_voisinage + 1)^2;		% Nombre de pixels dans un patch
% Les profondeurs
Z_VT = z(:,:,1);
nb_z = 1;
% Les normales (pour analyser des résultats avec vérité terrain)
N_VT = N(:,:,:,1);
% Les poses relatives des caméras
R_1_k = zeros(3,3,nb_images-1);
t_1_k = zeros(3,nb_images-1);
for k = 1:nb_images-1
	R_1_k(:,:,k) = R(:,:,k+1) * R(:,:,1)';
	t_1_k(:,k) = t(:,k+1) - R_1_k(:,:,k) * t(:,1);
end

f2 = figure;
f3 = figure;
f4 = figure;


% Filtrage des pixels considérés par le masque
[i_k, j_k] = find(masque(:,:,1));
ind_1		= sub2ind([nb_lignes nb_colonnes], i_k, j_k);
nb_pixels_etudies = size(ind_1,1);
P_k 		= zeros(3,nb_pixels_etudies,nb_images);
P_k(:,:,1) 	= [(j_k - offset - u_0) / facteur_k, (i_k - offset - v_0) / facteur_k, zeros(length(i_k), 1)].';
z_grossiers_estimes = zeros(nb_pixels_etudies,1);


%% Calcul des gradients
dx_I_k = zeros(size(I));
dy_I_k = zeros(size(I));
for k = 1:nb_images
	% Gradient non filtré
	[dx_I, dy_I] = gradient(I(:,:,k));
	% Sauvegarde des gradients
	dx_I_k(:,:,k) = dx_I;
	dy_I_k(:,:,k) = dy_I;
end
dx_I_1 = dx_I_k(:,:,1);
dy_I_1 = dy_I_k(:,:,1);
grad_I_x	= [dx_I_1(ind_1)'];
grad_I_y    = [dy_I_1(ind_1)'];

%% Construction du voisinage
voisinage_ligne = -rayon_voisinage*nb_lignes:nb_lignes:rayon_voisinage*nb_lignes;
voisinage_colonne = -rayon_voisinage:rayon_voisinage;
grille_voisinage = voisinage_ligne + voisinage_colonne';
grille_voisinage = grille_voisinage';

%% Mise en forme des normales
normale_theorique = [N_VT(ind_1)' ; N_VT(ind_1 + nb_pixels)' ; N_VT(ind_1 + 2*nb_pixels)'];

%% Boucle de reconstruction
erreurs	= 10*ones(nb_pixels_etudies, nb_z);
n_estimes = zeros(3, nb_pixels_etudies, nb_z);

tic
for indice_z = 1:nb_z

	% Affichage de la progression des calculs
	switch (affichage)
		case 'Iteration'
			fprintf('\r');
			fprintf("Progression : %d / %d",indice_z,nb_z);
		case 'Pourcentage'
			if mod(indice_z,round(nb_z/25)) == 0
				disp("Progression à " + int2str(indice_z/nb_z*100) + "%");
			end
	end

	% Sélection d'une profondeur
	valeur_z 	= Z_VT(ind_1);
	P_k(3,:,1) 	= Z_VT(ind_1);

	% Changements de repère
	for k = 1:nb_images-1
		P_k(:,:,k+1) = R_1_k(:,:,k) * P_k(:,:,1) + t_1_k(:,k);
		i_k(:,k+1) = (P_k(2,:,k+1) * facteur_k + offset + v_0).';
		j_k(:,k+1) = (P_k(1,:,k+1) * facteur_k + offset + u_0).';
	end

	% Vérification des pixels hors images
	condition_image = ones(nb_pixels_etudies,nb_images-1);
	for k = 1:nb_images-1
		condition_image(:,k) = i_k(:,k+1) > 0.5 & i_k(:,k+1) <= nb_lignes & j_k(:,k+1) > 0.5 & j_k(:,k+1) <= nb_colonnes;
	end

	% Calcul des gradients
	for k = 1:nb_images-1
		i_k(:,k+1) = ~condition_image(:,k) + condition_image(:,k) .* i_k(:,k+1);
		j_k(:,k+1) = ~condition_image(:,k) + condition_image(:,k) .* j_k(:,k+1);
		grad_I_x(k+1,:) = interp2(dx_I_k(:,:,k+1),j_k(:,k+1),i_k(:,k+1),interpolation)';
		grad_I_y(k+1,:) = interp2(dy_I_k(:,:,k+1),j_k(:,k+1),i_k(:,k+1),interpolation)';
	end

	% Calcul des numérateurs et dénominateurs
	numerateur_x = [];
	numerateur_y = [];
	denominateur = [];
	for k = 1:nb_images-1
		numerateur = [grad_I_x(1,:); grad_I_y(1,:)] ...
			- R_1_k(1:2,1:2,k)' * [grad_I_x(k+1,:); grad_I_y(k+1,:)];
		numerateur_x(k,:) = numerateur(1,:);
		numerateur_y(k,:) = numerateur(2,:);
		denominateur(k,:) = R_1_k(1:2,3,k)' * [grad_I_x(k+1,:); grad_I_y(k+1,:)];
	end
	clear numerateur;

	% Estimation de la normale
	p_truc = numerateur_x ./ denominateur;
	q_truc = numerateur_y ./ denominateur;
	normale = normales_medianes(p_truc,q_truc);
	n_estimes(:,:,indice_z) = normale;

	% Calcul du plan considéré
	d_equation_plan = sum(-P_k(:,:,1) .* normale,1);

	% Calcul de la transformation géométrique
	ind_decales = ind_1 + grille_voisinage(:)'; % Création de matrice avec 2 vecteurs
	[i_1_decales, j_1_decales] = ind2sub([nb_lignes, nb_colonnes], ind_decales);
	u_1_decales = (j_1_decales - offset - u_0) / facteur_k ;
	v_1_decales = (i_1_decales - offset - v_0) / facteur_k;

	normale_1 = repmat(normale(1,:)',1,taille_patch);
	normale_2 = repmat(normale(2,:)',1,taille_patch);
	normale_3 = repmat(normale(3,:)',1,taille_patch);
	z_1_decales = -(d_equation_plan' + normale_1.*u_1_decales + normale_2.*v_1_decales)./normale_3;

	% Reprojection du voisinage
	i_k_voisinage = zeros(nb_pixels_etudies, taille_patch, nb_images-1);
	j_k_voisinage = zeros(nb_pixels_etudies, taille_patch, nb_images-1);
	u_1_decales_vec = reshape(u_1_decales',1,nb_pixels_etudies*taille_patch);
	v_1_decales_vec = reshape(v_1_decales',1,nb_pixels_etudies*taille_patch);
	z_1_decales_vec = reshape(z_1_decales',1,nb_pixels_etudies*taille_patch);
	P_1_voisinage = [u_1_decales_vec ; v_1_decales_vec ; z_1_decales_vec];
	for k = 1:nb_images-1
		P_k_voisinage = R_1_k(:,:,k) * P_1_voisinage + t_1_k(:,k);
		P_k_voisinage_ok = cell2mat(mat2cell(P_k_voisinage,3,repmat(taille_patch,1,nb_pixels_etudies))');
		i_k_voisinage(:,:,k) = P_k_voisinage_ok(2:3:end,:) * facteur_k + offset + v_0;
		j_k_voisinage(:,:,k) = P_k_voisinage_ok(1:3:end,:) * facteur_k + offset + u_0;
	end

	% Calcul de l'erreur
	I_1_voisinage = interp2(I(:,:,1),j_1_decales,i_1_decales,interpolation);
	erreur_k = zeros(nb_pixels_etudies, nb_images-1);
	I_k = zeros(nb_lignes,nb_colonnes);
	for k = 1:nb_images-1
		I_k_voisinage = interp2(I(:,:,k+1),j_k_voisinage(:,:,k),i_k_voisinage(:,:,k),interpolation);
		erreur_k(:,k) = condition_image(:,k).*sum((I_1_voisinage-I_k_voisinage),2);
		if (k == 2)
			I_k(ind_1) = I_k_voisinage;
		end
	end
	switch (estimateur)
		case 'MSE'
			erreurs(:,indice_z) = (1 ./ sum(condition_image,2)) .* sum(erreur_k.^2,2);
		case 'Robuste'
			erreurs(:,indice_z) = (1 ./ sum(condition_image,2)) .* (1 - exp(-sum(erreur_k.^2,2)/0.2^2));
	end

end

fprintf('\n');
toc

%% Résultats

taille_voisi = 2*rayon_voisinage + 1;

figure(f2);
imagesc(I(:,:,1));
axis equal;
title("Image 1");

figure(f3);
imagesc(I_k);
ax = gca;
c.CLim = [1e-7 1e-1];
axis equal;
title("Image k");

figure(f4);
imagesc(I(:,:,1)-I_k);
axis equal;
ax = gca;
c.CLim = [1e-7 1e-1];
title("Différence entre 2 images");
colorbar;




















% Sélections des profondeurs avec l'erreur minimale
erreurs_corrigees = erreurs;
[~,indices_min] = min(erreurs_corrigees,[],2);
z_in = Z_VT(ind_1);
z_estime = nan(nb_lignes, nb_colonnes);
z_estime(ind_1) = z_in;

% Calcul des erreurs de reconstruction
Z_1_in = Z_VT(ind_1);
erreur_z = abs(Z_1_in - z_in);

% Sélection des normales
n_totales = zeros(3, nb_pixels);
n_totales_ind = zeros(3, nb_pixels_etudies);
for k = 1:3
	n_totales(k,ind_1) = n_estimes(sub2ind(size(n_estimes), k * ones(nb_pixels_etudies,1), transpose(1:nb_pixels_etudies), indices_min));
	n_totales_ind(k,:) = n_totales(k,ind_1);
end

% Sélections des erreurs angulaires
angles = zeros(nb_lignes, nb_colonnes);
angles(ind_1) = abs((180/pi) * atan2(vecnorm(cross(normale_theorique,n_totales_ind)),dot(normale_theorique,n_totales_ind)));
angles_ind = angles(ind_1);

% Affichage des erreurs angulaires
erreur_angle_moy = mean(angles_ind(:));
erreur_angle_med = median(angles_ind(:));
espace_z_suivant = 0;
