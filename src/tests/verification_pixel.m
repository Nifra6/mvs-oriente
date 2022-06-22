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
liste_surface = ["gaussienne_1"];
liste_rayon_voisinage = 1;
nombre_iteration = 1;
%liste_ecart_type_I = [0:0.5:3];
liste_ecart_type_I = -2.5;
%liste_ecart_type_grad = [0:10];
liste_ecart_type_grad = -1;
filtrage = 0;
liste_nombre_vues = 2;
liste_nombre_profondeur_iteration = [5000];
utilisation_mediane_normale = 0;
i_pixel = 250;
j_pixel = 250;
%i_pixel = 210;
%j_pixel = 170;
%i_pixel = 233;
%j_pixel = 394;

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
							utilisation_profondeur_GT = 0;
							utilisation_normale_GT = 0;
							[z_estime_mvsm,score,echantillons_z,erreur_z_mvsm,espace_z,normales_mvsm,erreur_angle_moy,erreur_angle_med,i_2,j_2] = mvs_modifie_bis(pixel_considere,premiere_iteration,surface,nombre_vues,rayon_voisinage,ecart_type_I,ecart_type_grad,nombre_z,z_estime_mvsm,espace_z,utilisation_profondeur_GT,utilisation_normale_GT,utilisation_mediane_normale);
							disp("MVS modifié, valeur VT :");
							utilisation_profondeur_GT = 1;
							utilisation_normale_GT = 0;
							[z_GT,score_GT,~,erreur_z_mvsm,~,~,~,~,i_2_GT,j_2_GT] = mvs_modifie_bis(pixel_considere,premiere_iteration,surface,nombre_vues,rayon_voisinage,ecart_type_I,ecart_type_grad,nombre_z,z_estime_mvsm,espace_z,utilisation_profondeur_GT,utilisation_normale_GT,utilisation_mediane_normale);
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
						legend("Échantillonage MVSm", "Minimum trouvé", "Échantillon proche de la VT","5 minimums", "Profondeur VT");
						title("Surface " + surface);
						% Affichage du pixel
						load("../../data/"+"simulateur_"+surface+"_formate.mat");
						figure;
						imshow(I(:,:,1));
						hold on;
						plot(j_pixel,i_pixel,'b+','LineWidth',3);
						figure;
						imshow(I(:,:,2));
						hold on;
						plot([j_2(1) j_2(end)],[i_2(1) i_2(end)],'b','LineWidth',1);
						plot(j_2(indice_mini),i_2(indice_mini),'mo','LineWidth',3,'MarkerSize',10);
						plot(j_2_GT,i_2_GT,'go','LineWidth',3,'MarkerSize',10);
						plot(j_2(indice_mins),i_2(indice_mins),'r+','LineWidth',3);
						plot(j_2(indice_nearest),i_2(indice_nearest),'k+','LineWidth',3);
						hold off;
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
