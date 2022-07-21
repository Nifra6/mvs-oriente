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
valeur_bruitage = 2;
liste_surface = ["gaussienne_1_bruitee_" + int2str(valeur_bruitage), "gaussienne_1_pepper_bruitee_" ...
	+ int2str(valeur_bruitage), "gaussienne_2_bruitee_" + int2str(valeur_bruitage), "sinc_1_bruitee_" ...
	+ int2str(valeur_bruitage)];
liste_surface = ["gaussienne_1", "gaussienne_1_pepper", "gaussienne_2", "sinc_1"];
%liste_surface = ["gaussienne_1_bis", "gaussienne_1_pepper_bis", "gaussienne_2_bis", "sinc_1_bis"];
liste_surface = ["plan"];
liste_rayon_voisinage = 15;
nombre_iteration = 1;
%liste_ecart_type_I = [0:0.5:3];
liste_ecart_type_I = -2.5;
%liste_ecart_type_grad = [0:10];
liste_ecart_type_grad = -1;
filtrage = 0;
liste_nombre_vues = 2;
liste_nombre_profondeur_iteration = [5000];
utilisation_mediane_normale = 1;
i_pixel = 250;
j_pixel = 250;
i_pixel = 210;
j_pixel = 170;
i_pixel = 233;
j_pixel = 394;
i_pixel = 17;
j_pixel = 192;
%i_pixel = 350; %prb tracé ?
%j_pixel = 400;
i_pixel = 360;
j_pixel = 170;
i_pixel = 253;
j_pixel = 187;
%hehe
i_pixel = 91;
j_pixel = 341;
i_pixel = 401;
j_pixel = 371;
i_pixel = 131;
j_pixel = 471;
% hehe 2
i_pixel = 171;
j_pixel = 51;
i_pixel = 221;
j_pixel = 61;
i_pixel = 281;
j_pixel = 281;
i_pixel = 371;
j_pixel = 101;
i_pixel = 351;
j_pixel = 131;
%i_pixel = 111;
%j_pixel = 161;
%i_pixel = 141;
%j_pixel = 161;
%i_pixel = 231;
%j_pixel = 301;
i_pixel = 251;
j_pixel = 311;
% hehe 3
%i_pixel = 311;
%j_pixel = 201;
% oui
%i_pixel = 287;
%j_pixel = 227;

%% Variables
pixel_considere = sub2ind([500 500],i_pixel,j_pixel);
nb_surface = size(liste_surface,2);
nb_rayon_voisinage = size(liste_rayon_voisinage,2);
nb_nombre_profondeur = size(liste_nombre_profondeur_iteration,2);
nb_nombre_vues = size(liste_nombre_vues,2);
nb_ecart_type_grad = size(liste_ecart_type_grad,2);
nb_ecart_type_I = size(liste_ecart_type_I,2);

%% Affichage des données fixes
disp("=================");
if (nb_surface == 1)
	disp("Surface utilisée : " + liste_surface(1));
end
if (nb_nombre_vues == 1)
	disp("Nombre de vues utilisées : " + int2str(liste_nombre_vues(1)));
end
if (nb_ecart_type_grad == 1)
	disp("Écart type pour les gradients utilisé : " + num2str(liste_ecart_type_grad(1)));
end
if (nb_ecart_type_I == 1)
	disp("Écart type pour la cohérence photométrique utilisé : " + num2str(liste_ecart_type_I(1)));
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


%% Modifications noms de fichiers
if (utilisation_mediane_normale)
	fichier_mediane = "__normales_medianes";
else
	fichier_mediane = "";
end


%% Algorithme
for i_ecart_type_I = 1:nb_ecart_type_I
	ecart_type_I = liste_ecart_type_I(i_ecart_type_I);
	if (nb_ecart_type_I > 1)
		disp("-------------------- Écart type pour la cohérence photométrique utilisé : " + ecart_type_I);
	end
	for i_ecart_type_grad = 1:nb_ecart_type_grad
		ecart_type_grad = liste_ecart_type_grad(i_ecart_type_grad);
		if (nb_ecart_type_grad > 1)
			disp("----------------- Écart type pour les gradients utilisé : " + ecart_type_grad);
		end
		if (ecart_type_grad >= 0 && filtrage)
			fichier_bruite = "__bruite_" + int2str(valeur_bruitage) + "__filtre_I_" + num2str(ecart_type_I) + "__filtre_grad_" + num2str(ecart_type_grad);
		else
			fichier_bruite = "";
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
							disp("MVS modifié, valeur estimée :");
							utilisation_profondeur_GT = 1;
							utilisation_normale_GT = 0;
							[z_estime_mvsm,score,echantillons_z,erreur_z_mvsm,espace_z,normales_mvsm,erreur_angle_moy,erreur_angle_med,i_2,j_2, patch_I_1, patch_I_2] = mvs_modifie_bis(pixel_considere,premiere_iteration,surface,nombre_vues,rayon_voisinage,ecart_type_I,ecart_type_grad,nombre_z,z_estime_mvsm,espace_z,utilisation_profondeur_GT,utilisation_normale_GT,utilisation_mediane_normale);
							disp("MVS modifié, valeur VT :");
							utilisation_profondeur_GT = 1;
							utilisation_normale_GT = 1;
							[z_GT,score_GT,~,erreur_z_mvsm,~,normale_VT,~,~,i_2_GT,j_2_GT,patch_I_1_GT, patch_I_2_GT] = mvs_modifie_bis(pixel_considere,premiere_iteration,surface,nombre_vues,rayon_voisinage,ecart_type_I,ecart_type_grad,nombre_z,z_estime_mvsm,espace_z,utilisation_profondeur_GT,utilisation_normale_GT,utilisation_mediane_normale);
						end
						% Préparation des résultats
						[~,indice_nearest] = min(abs(echantillons_z - z_GT));
						[~,indice_mini] = min(score);
						[~,indice_mins] = mink(score,5);
						max_value = max([max(score) ; score_GT]);
						% Affichage des résultats
						figure;
						plot(echantillons_z',score,'b.');
						hold on;
						plot([echantillons_z(indice_mini) ; echantillons_z(indice_mini)],[0 ; max_value],'m','LineWidth',2);
						plot([echantillons_z(indice_nearest) ; echantillons_z(indice_nearest)],[0 ; max_value],'k');
						%plot([echantillons_z(indice_mins) ; echantillons_z(indice_mins)],[0 ; max_value],'r');
						plot(z_GT,score_GT,'g+','LineWidth',2);
						hold off;
						legend("Échantillonage MVSm", "Minimum trouvé", "Échantillon proche de la VT", "Profondeur VT");
						title("Surface " + surface);
						xlabel("Profondeur échantillonée");
						ylabel("Score de cohérence photométrique");

						% Affichage du pixel
						load("../../data/"+"simulateur_"+surface+"_formate.mat");
						figure;
						imshow(I(:,:,1));
						hold on;
						plot(j_pixel,i_pixel,'b+','LineWidth',3);
						title("Image de référence");
						hold off;

						figure;
						imshow(I(:,:,2));
						hold on;
						plot([j_2(1) j_2(end)],[i_2(1) i_2(end)],'b','LineWidth',2);
						plot(j_2(indice_mini),i_2(indice_mini),'mo','LineWidth',3,'MarkerSize',10);
						plot(j_2_GT,i_2_GT,'go','LineWidth',3,'MarkerSize',10);
						plot(j_2(indice_mins),i_2(indice_mins),'r+','LineWidth',3);
						plot(j_2(indice_nearest),i_2(indice_nearest),'k+','LineWidth',3);
						title("Image témoin");
						hold off;

						diff_VT = patch_I_1 - patch_I_2(:,:,indice_nearest);
						diff_estime = patch_I_1 - patch_I_2(:,:,indice_mini);
						limite_min = min([min(diff_VT(:)) min(diff_estime(:))]);
						limite_max = max([max(diff_VT(:)) max(diff_estime(:))]);

						figure;
						subplot(2,3,1);
						imshow(patch_I_1);
						title("Image de référence");
						subplot(2,3,2);
						imshow(patch_I_2(:,:,indice_nearest));
						title("Image témoin, VT");
						subplot(2,3,3);
						imshow(patch_I_2(:,:,indice_mini));
						title("Image témoin, estimé");
						subplot(2,3,4);
						imshow(patch_I_1);
						title("Image de référence");
						subplot(2,3,5);
						imagesc(diff_VT,[limite_min limite_max]);
						colorbar;
						title("Résidu témoin, VT");
						subplot(2,3,6);
						%imshow(abs(patch_I_1 - patch_I_2(:,:,indice_mini)));
						imagesc(diff_estime,[limite_min limite_max]);
						title("Résidu témoin, estimé");


						format long
						disp("Localisation du pixel reprojeté trouvé :");
						[j_2(indice_mini) i_2(indice_mini)]
						disp("Localisation du pixel reprojeté proche de la VT :");
						[j_2(indice_nearest) i_2(indice_nearest)]
						disp("Localisation du pixel reprojeté VT :");
						[j_2_GT i_2_GT]

						angle_normale(normale_VT,normales_mvsm)





						% Sauvegarde des résultats
						nom_fichier = "Résultat_" + surface + "__nb_vues_" + int2str(nombre_vues) ...
							+ "__patch_" + int2str(taille_patch) + "x" + int2str(taille_patch) ...
							+ "__nb_profondeur_" + int2str(nombre_profondeur_iteration) ...
							+ fichier_bruite + fichier_mediane + ".mat";
						path = "../../result/tests/";
						%save(path+nom_fichier,"surface","nombre_vues","taille_patch","nombre_profondeur_iteration","z_estime_mvs","z_estime_mvsm","erreur_z_mvs","erreur_z_mvsm","normales_mvs","normales_mvsm","erreur_angle_moy","erreur_angle_med");
						%disp("Enregistrement sous : " + nom_fichier);
					end
					disp("------------------------------------");
				end
			end
		end
	end
end
