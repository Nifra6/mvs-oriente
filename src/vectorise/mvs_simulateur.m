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
rayon_voisinage = 3;			% Rayon du voisinage carré à prendre en compte

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
rayon_masque = 1;
taille_masque = (2*rayon_masque+1)^2;
[x,y] = meshgrid(-rayon_masque:rayon_masque,-rayon_masque:rayon_masque);
u_x = 0; u_y = 0; sigma = taille_masque/4;
filtre = 1./(2*pi*sigma^2) .* exp(((x-u_x).^2+(y-u_y).^2)./(2*sigma^2));
filtre = filtre / norm(filtre);
I_filtre = zeros(size(I));
for k = 1:nombre_images
	I_filtre(:,:,k) = conv2(I(:,:,k),filtre,'same');
end
%I_filtre = I;


%% Construction du voisinage
voisinage_ligne = -rayon_voisinage*nombre_lignes:nombre_lignes:rayon_voisinage*nombre_lignes;
voisinage_colonne = -rayon_voisinage:rayon_voisinage;
grille_voisinage = voisinage_ligne + voisinage_colonne';
grille_voisinage = grille_voisinage';

%% Boucle de reconstruction
erreurs	= 10*ones(nombre_pixels_etudies, nombre_z);

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
		i_k(:,k+1) = round(i_k(:,k+1));
		j_k(:,k+1) = round(j_k(:,k+1));
		ind(:,k+1) = sub2ind([nombre_lignes nombre_colonnes], i_k(:,k+1), j_k(:,k+1));
	end

	% Calcul de la transformation géométrique
	ind_decales = ind_1 + grille_voisinage(:)'; % Création de matrice avec 2 vecteurs
	[i_1_decales, j_1_decales] = ind2sub([nombre_lignes, nombre_colonnes], ind_decales);
	u_1_decales = (j_1_decales-u_0) / facteur_k ;
	v_1_decales = (i_1_decales-v_0) / facteur_k;
	
	z_1_decales = valeur_z * ones(size(u_1_decales));

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
z = nan(nombre_lignes, nombre_colonnes);
z(ind_1) = z_in;

% Mesures
disp("==============")
disp("Mesure relative de profondeur")
sum(abs(Z_1(ind_1) - z(ind_1)),'all') / nombre_pixels_etudies
ecart_moyen = sum(Z_1(ind_1) - z(ind_1)) / size(z_in,1);
disp("Mesure relative de forme")
sum(abs(Z_1(ind_1) - (z(ind_1) + ecart_moyen)),'all') / size(z_in,1)

% Préparation
X = 1:nombre_colonnes;
X = (X - u_0) / facteur_k;
Y = 1:nombre_lignes;
Y = (Y - v_0) / facteur_k;
[X,Y] = meshgrid(X,Y);

% Affichage
figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
sl = surfl(X,Y,-z,s);
sl.EdgeColor = 'none';
grid off;
colormap gray;
axis equal;
