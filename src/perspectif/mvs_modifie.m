% Reconstruire une surface via l'algorithme de Multi-View Stereo modifié (ou MVSm).
% Utilisé par lancement_test.m

function [z_estime,erreur_z,espace_z_suivant,n_totales_ind,erreur_angle_moy,erreur_angle_med] = mvs_modifie(premiere_iteration,surface,nb_vues,rayon_voisinage,sigma_filtre_I,sigma_filtre_grad,nb_z,z_precedent,espace_z,utilisation_profondeurs_GT,utilisation_normale_GT,grille_pixels)

	%% Paramètres
	interpolation 	= 'linear';			% Type d'interpolation utilisée
	estimateur		= 'MSE';			% Estimateur utilisé pour l'évaluation des erreurs photométriques
	affichage 		= 'Iteration';		% Type d'affichage de la progression de l'algorithme
	offset 			= 0.5;				% Décalage spatial entre les indices des pixels et leur coordonnées


	%% Données
	% Chargement des fonctions utiles
	addpath(genpath("../toolbox/"));
	% Chargement des données
	path = "../../data/perspectif/";
	nom_fichier = "simulateur_" + surface + "_formate.mat";
	load(path+nom_fichier);
	% Nombres d'images et de pixels considérés
	nb_pixels = nb_lignes * nb_colonnes;
	nb_images = nb_vues;
	taille_patch  = (2*rayon_voisinage + 1)^2;		% Nombre de pixels dans un patch
	% Les profondeurs
	Z_VT = z(:,:,1);
	if (premiere_iteration)
		if (surface == "plan_bis")
			valeurs_z = linspace(4,6,nb_z);
		else
			valeurs_z = linspace(min(Z_VT,[],'all'),max(Z_VT,[],'all'),nb_z);
		end
	else
		valeurs_z = linspace(-espace_z,espace_z,nb_z);	% Valeurs de profondeurs testées
	end
	if (utilisation_profondeurs_GT)
		nb_z = 1;
	end
	% Les normales (pour analyser des résultats avec vérité terrain)
	N_VT = N(:,:,:,1);
	% Les poses relatives des caméras
	R_1_k = zeros(3,3,nb_images-1);
	t_1_k = zeros(3,nb_images-1);
	for k = 1:nb_images-1
		R_1_k(:,:,k) = R(:,:,k+1) * R(:,:,1)';
		t_1_k(:,k) = t(:,k+1) - R_1_k(:,:,k) * t(:,1);
	end
	% La matrice inverse de calibrage
	K_inv = inv(K);
	% Modifications du masque (pour correspondre aux patchs utilisés)
	masque(1:rayon_voisinage,:,1) = 0;
	masque(end-rayon_voisinage:end,:,1) = 0;
	masque(:,1:rayon_voisinage,1) = 0;
	masque(:,end-rayon_voisinage:end,1) = 0;
	% Filtrage des pixels considérés par le masque
	[i_k, j_k]  = find(masque(:,:,1));
	ind_1		= sub2ind([nb_lignes nb_colonnes], i_k, j_k);
	% Utilisation d'une grille régulière de pixels
	if (grille_pixels > 0)
		indices_grilles = (mod(i_k,grille_pixels) == 1) & (mod(j_k,grille_pixels) == 1);
		ind_1 = ind_1(find(indices_grilles));
		i_k = i_k(find(indices_grilles));
		j_k = j_k(find(indices_grilles));
	end
	nb_pixels_etudies = size(ind_1,1);
	P_k 		= zeros(3,nb_pixels_etudies,nb_images);
	u_k 		= zeros(nb_pixels_etudies,nb_images);
	v_k 		= zeros(nb_pixels_etudies,nb_images);
	u_k(:,1)	= j_k - offset;
	v_k(:,1)	= i_k - offset;
	p_1			= [u_k(:,1) , v_k(:,1) , ones(nb_pixels_etudies,1)]';
	if (premiere_iteration)
		z_grossiers_estimes = zeros(nb_pixels_etudies,1);
	else
		z_grossiers_estimes = z_precedent(ind_1);
	end


	%% Calcul du filtre
	filtrage = sigma_filtre_grad >= 0 | sigma_filtre_I >= 0;
	if (filtrage)
		u_x = 0; u_y = 0;
		if (sigma_filtre_I > 0)
			cote_masque_I = ceil(4*sigma_filtre_I);
			filtre_I = fspecial('gauss',cote_masque_I,sigma_filtre_I);
			filtre_I = filtre_I / sum(filtre_I(:));
			I_filtre = zeros(size(I));
			for k = 1:nb_images
				I_filtre(:,:,k) = conv2(I(:,:,k),filtre_I,'same');
			end
		else
			I_filtre = I;
		end
		if (sigma_filtre_grad > 0)
			cote_masque_grad = ceil(4*sigma_filtre_grad);
			filtre_grad = fspecial('gauss',cote_masque_grad,sigma_filtre_grad);
			filtre_grad = filtre_grad / sum(filtre_grad(:));
			dx_conv = [0 0 0 ; 1 0 -1 ; 0 0 0];
			dy_conv = [0 1 0 ; 0 0 0 ; 0 -1 0];
			dx_filtre = conv2(filtre_grad,dx_conv);
			dy_filtre = conv2(filtre_grad,dy_conv);
		end
	else
		I_filtre = I;
	end


	%% Calcul des gradients
	dx_I_k = zeros(size(I));
	dy_I_k = zeros(size(I));
	for k = 1:nb_images
		if (sigma_filtre_grad > 0)
			% Gradient filtré
			dx_I = conv2(I(:,:,k),dx_filtre,'same');
			dy_I = conv2(I(:,:,k),dy_filtre,'same');
		else
			% Gradient non filtré
			[dx_I, dy_I] = gradient(I(:,:,k));
		end
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

	%% Mise en forme des normales théoriques
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
		valeur_z 	= z_grossiers_estimes + valeurs_z(indice_z);
		if (utilisation_profondeurs_GT)
			Z = repmat(Z_VT(ind_1)',3,1);
		else
			Z = repmat(valeur_z',3,1);
		end
		P_k(:,:,1) = Z .* (K_inv * p_1);

		% Changements de repère
		for k = 1:nb_images-1
			P_k(:,:,k+1) = R_1_k(:,:,k) * P_k(:,:,1) + t_1_k(:,k);
			p_k = (K * P_k(:,:,k+1)) ./ P_k(3,:,k+1);
			u_k(:,k+1) = p_k(1,:)';
			v_k(:,k+1) = p_k(2,:)';
			i_k(:,k+1) = v_k(:,k+1) + offset;
			j_k(:,k+1) = u_k(:,k+1) + offset;
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
		grad_I_1 = [grad_I_x(1,:); grad_I_y(1,:)];
		%P_monde = R(:,:,1)' * P_k(:,:,1) - R(:,:,1)' * t(:,1);
		P_monde = P_k(:,:,1);
		for k = 1:nb_images-1
			%w_i1k = cross(repmat(R(:,1,k+1),1,nb_pixels_etudies),P_k(:,:,k+1));
			%w_i2k = cross(repmat(R(:,2,k+1),1,nb_pixels_etudies),P_k(:,:,k+1));
			%w_i3k = cross(repmat(R(:,3,k+1),1,nb_pixels_etudies),P_k(:,:,k+1));

			w_i1k = cross(repmat(R_1_k(:,1,k),1,nb_pixels_etudies),P_k(:,:,k+1));
			w_i2k = cross(repmat(R_1_k(:,2,k),1,nb_pixels_etudies),P_k(:,:,k+1));
			w_i3k = cross(repmat(R_1_k(:,3,k),1,nb_pixels_etudies),P_k(:,:,k+1));

			grad_I_k = [grad_I_x(k+1,:); grad_I_y(k+1,:)];

			matrice_numerateur = zeros(2,2,nb_pixels_etudies);
			matrice_numerateur(1,1,:) = w_i1k(2,:);
			matrice_numerateur(1,2,:) = -w_i1k(1,:);
			matrice_numerateur(2,1,:) = w_i2k(2,:);
			matrice_numerateur(2,2,:) = -w_i2k(1,:);


			numerateur_bis = zeros(2,nb_pixels_etudies);
			for pixel = 1:nb_pixels_etudies
				numerateur_bis(:,pixel) = matrice_numerateur(:,:,pixel) * grad_I_k(:,pixel);
			end

			numerateur = (P_k(3,:,k+1).^2) .* P_monde(3,:) .* grad_I_1 + (P_monde(3,:).^2) .* numerateur_bis;
			numerateur_x(k,:) = numerateur(1,:);
			numerateur_y(k,:) = numerateur(2,:);			

			denominateur(k,:) = P_monde(3,:).^2 .* dot([-w_i3k(2,:) ; w_i3k(1,:)],grad_I_k) ...
				+ P_k(3,:,k+1).^2 .* dot(P_monde(1:2,:),grad_I_1);
		end
		clear coeff_z_1 coeff_z_k grad_I_1 grad_I_k numerateur;

		% Estimation de la normale
		p_truc = numerateur_x ./ denominateur;
		q_truc = numerateur_y ./ denominateur;
		normale = normales_medianes(p_truc,q_truc);

		if (utilisation_normale_GT)
			normale = normale_theorique;
		end
		n_estimes(:,:,indice_z) = normale;

		% Calcul du plan considéré
		d_equation_plan = sum(-P_k(:,:,1) .* normale,1);

		% Calcul de la transformation géométrique
		ind_decales = ind_1 + grille_voisinage(:)'; % Création de matrice avec 2 vecteurs
		[i_1_decales, j_1_decales] = ind2sub([nb_lignes, nb_colonnes], ind_decales);
		u_1_decales = j_1_decales - offset;
		v_1_decales = i_1_decales - offset;

		% Reprojection du voisinage
		i_k_voisinage = zeros(nb_pixels_etudies, taille_patch, nb_images-1);
		j_k_voisinage = zeros(nb_pixels_etudies, taille_patch, nb_images-1);
		u_1_decales_vec = reshape(u_1_decales',1,nb_pixels_etudies*taille_patch);
		v_1_decales_vec = reshape(v_1_decales',1,nb_pixels_etudies*taille_patch);
		for k = 1:nb_images-1
			for pixel = 1:nb_pixels_etudies

				homographie = K * (R_1_k(:,:,k) - t_1_k(:,k) * normale(:,pixel)' / d_equation_plan(pixel)) * K_inv;	
				p_k_voisinage = homographie * Z(1,pixel) * [u_1_decales(pixel,:) ; v_1_decales(pixel,:) ; ones(1,taille_patch)];
				u_k_voisinage = p_k_voisinage(1,:) ./ p_k_voisinage(3,:);
				v_k_voisinage = p_k_voisinage(2,:) ./ p_k_voisinage(3,:);

				i_k_voisinage(pixel,:,k) = v_k_voisinage + offset;
				j_k_voisinage(pixel,:,k) = u_k_voisinage + offset;
			end
		end

		% Calcul de l'erreur
		I_1_voisinage = interp2(I_filtre(:,:,1),j_1_decales,i_1_decales,interpolation);
		erreur_k = zeros(nb_pixels_etudies, nb_images-1);
		for k = 1:nb_images-1
			I_k_voisinage = interp2(I_filtre(:,:,k+1),j_k_voisinage(:,:,k),i_k_voisinage(:,:,k),interpolation);
			erreur_k(:,k) = condition_image(:,k).*sum((I_1_voisinage-I_k_voisinage),2);
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

	% Sélections des profondeurs avec l'erreur minimale
	erreurs_corrigees = erreurs;
	[~,indices_min] = min(erreurs_corrigees,[],2);
	if (utilisation_profondeurs_GT)
		z_in = Z_VT(ind_1);
	else
		z_in = z_grossiers_estimes + transpose(valeurs_z(indices_min));
	end
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

	% Sauvegarde de la zone de recherche
	if (utilisation_profondeurs_GT)
		espace_z_suivant = 0;
	else
		espace_z_suivant = abs(valeurs_z(2) - valeurs_z(1));
	end

end
