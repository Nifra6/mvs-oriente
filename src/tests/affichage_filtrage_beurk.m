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
liste_surface = ["gaussienne_1_bruitee_10", "gaussienne_2_bruitee_10", "sinc_1_bruitee_10", "gaussienne_1_pepper_bruitee_10"];
liste_rayon_voisinage = 4;
nombre_iteration = 1;
liste_ecart_type = 6;
%liste_ecart_type = [3 4]; % temp
liste_nombre_vues = [4];
%liste_nombre_vues = [9]; % temp
liste_nombre_profondeur_iteration = [5000];

%% Variables
nb_surface = size(liste_surface,2);
nb_rayon_voisinage = size(liste_rayon_voisinage,2);
nb_nombre_profondeur = size(liste_nombre_profondeur_iteration,2);
nb_nombre_vues = size(liste_nombre_vues,2);
nb_ecart_type = size(liste_ecart_type,2);

%% Affichage des données fixes
disp("=================");
if (nb_surface == 1)
	disp("Surface utilisée : " + liste_surface(1));
end
if (nb_nombre_vues == 1)
	disp("Nombre de vues utilisées : " + int2str(liste_nombre_vues(1)));
end
if (nb_ecart_type == 1)
	disp("Écart type utilisé : " + int2str(liste_ecart_type(1)));
end

if (nb_rayon_voisinage == 1)
	rayon_voisinage = liste_rayon_voisinage(1);
	taille_patch = 2*rayon_voisinage + 1;
	disp("Taille de patch voisinage : " + int2str(taille_patch) + "x" + int2str(taille_patch));
end
if (nb_nombre_profondeur == 1)
	disp("Nombre de profondeurs testées : " + int2str(liste_nombre_profondeur_iteration(1)));
end
disp("=================");

%% Algorithme
for i_ecart_type = 1:nb_ecart_type
	ecart_type = liste_ecart_type(i_ecart_type);
	if (nb_ecart_type > 1)
		disp("----------------- Écart type utilisé : " + ecart_type);
	end
	if (ecart_type >= 0)
		fichier_bruite = "__bruite_" + int2str(10) + "__filtre_" + int2str(ecart_type);
	end

	for i_surface = 1:nb_surface
		surface = liste_surface(i_surface);
		if (nb_surface > 1)
			disp("-------------- Surface utilisée : " + surface);
		end
		for i_nombre_vues = 1:nb_nombre_vues
			nombre_vues = liste_nombre_vues(i_nombre_vues);
			if (nb_nombre_vues > 1)
				disp("----------- Nombre de vues utilisées : " + int2str(nombre_vues));
			end
			for i_rayon_voisinage = 1:nb_rayon_voisinage
				rayon_voisinage = liste_rayon_voisinage(i_rayon_voisinage);
				taille_patch = 2*rayon_voisinage + 1;
				if (nb_rayon_voisinage > 1)
					disp("-------- Taille de patch voisinage : " + int2str(taille_patch) + "x" + int2str(taille_patch));
				end
				for i_nombre_profondeur = 1:nb_nombre_profondeur
					nombre_profondeur_iteration = liste_nombre_profondeur_iteration(i_nombre_profondeur);
					if (nb_nombre_profondeur > 1)
						disp("----- Nombre de profondeurs testées : " + int2str(nombre_profondeur_iteration));
					end
					for i = 1:nombre_iteration
						if (nombre_iteration > 1)
							disp("-- Itération " + int2str(i) + " / " + int2str(nombre_iteration));
						end
						% Paramétrage
						premiere_iteration = (i == 1);
						if (premiere_iteration)
							z_estime_mvs = 0;
							z_estime_mvsm = 0;
							nombre_z = nombre_profondeur_iteration + 1;
							espace_z = 0;
						else
							nombre_z = 2 * nombre_profondeur_iteration + 1;
						end
						% Exécution
						%% Paramètres
						interpolation 	= 'linear';			% Type d'interpolation
						estimateur		= 'MSE';			% Estimateur utilisé pour l'évaluation des erreurs
						affichage 		= 'Iteration';		% Type d'affichage de la progression
						affichage_debug = 0;				% Affichage d'informations diverses
						%rayon_voisinage = 1;				% Rayon du voisinage carré à prendre en compte
						filtrage 		= ecart_type >= 0;% Utilisation d'un filtrage gaussien
						grille_pixels	= 10;

						%% Données
						% Fichier des données
						path = "../../data/";
						fichier_surface = "simulateur_" + surface + "_formate.mat";
						load(path+fichier_surface);
						% Taille des images
						nombre_pixels = nombre_lignes * nombre_colonnes;
						nombre_images = nombre_vues;
						taille_patch  = (2*rayon_voisinage + 1)^2;	% Nombre de pixels dans un patch
						% Les profondeurs
						Z_1 = z(:,:,1);
						if (premiere_iteration)
							valeurs_z = linspace(min(Z_1,[],'all'),max(Z_1,[],'all'),nombre_z);
						else
							Z_1_estime = z_precedent;
							valeurs_z = linspace(-espace_z,espace_z,nombre_z);	% Valeurs de profondeurs testées
						end
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
						% Utilisation d'une grille régulière de pixels
						if (grille_pixels > 0)
							indices_grilles = (mod(i_k,grille_pixels) == 1) & (mod(j_k,grille_pixels) == 1);
							ind_1 = ind_1(find(indices_grilles));
							i_k = i_k(find(indices_grilles));
							j_k = j_k(find(indices_grilles));
						end
						nombre_pixels_etudies = size(ind_1,1);
						P_k 		= zeros(3,nombre_pixels_etudies,nombre_images);
						P_k(:,:,1) 	= [(j_k - u_0) / facteur_k, (i_k - v_0) / facteur_k, zeros(length(i_k), 1)].';
						if (premiere_iteration)
							z_grossiers_estimes = zeros(nombre_pixels_etudies,1);
						else
							z_grossiers_estimes = Z_1_estime(ind_1);
						end


						%% Calcul du filtre
						if (filtrage)
							% Analytique
							%rayon_masque = floor(ecart_type * 4);
							%taille_masque = (2*rayon_masque+1)^2;
							%[x,y] = meshgrid(-rayon_masque:rayon_masque,-rayon_masque:rayon_masque);
							%u_x = 0; u_y = 0;
							%filtre = 1./(2*pi*ecart_type^2) .* exp(((x-u_x).^2+(y-u_y).^2)./(2*ecart_type^2));
							%filtre = filtre / sum(filtre(:));
							%dx_filtre = -(x-u_x)./(2*pi*ecart_type^4) .* exp(((x-u_x).^2+(y-u_y).^2)./(2*ecart_type^2));
							%dy_filtre = -(y-u_y)./(2*pi*ecart_type^4) .* exp(((x-u_x).^2+(y-u_y).^2)./(2*ecart_type^2));
							%dx_filtre = dx_filtre;
							%dy_filtre = dy_filtre;

							% Fonctions Matlab
							u_x = 0; u_y = 0;
							cote_masque = ceil(4*ecart_type);
							filtre = fspecial('gauss',cote_masque,ecart_type);
							filtre = filtre / sum(filtre(:));
							dx_conv = [0 0 0 ; 1 0 -1 ; 0 0 0];
							dy_conv = [0 1 0 ; 0 0 0 ; 0 -1 0];
							dx_filtre = conv2(filtre,dx_conv);
							dy_filtre = conv2(filtre,dy_conv);
							dx_filtre = dx_filtre;
							dy_filtre = dy_filtre;

							% Filtrage de l'image
							I_filtre = zeros(size(I));
							for k = 1:nombre_images
								I_filtre(:,:,k) = conv2(I(:,:,k),filtre,'same');
							end
						else
							I_filtre = I;
						end
						figure;
						imshow(I(:,:,1));
						figure;
						imshow(I_filtre(:,:,1));
						pause;
						close all;
					end
				end
				disp("------------------------------------");
			end
		end
	end
end
