% Lancer un ensemble de reconstruction sur des images prises avec une caméra perspective via
% MVS classique et MVS modifié.

%% Trucs de Matlab
% Clear
clear;
close all;
% Paramètres d'affichage
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

%% Paramètres
valeur_bruitage = 8;
liste_surface = ["gaussienne_1_bruitee_" + int2str(valeur_bruitage), "gaussienne_1_pepper_bruitee_" ...
	+ int2str(valeur_bruitage), "gaussienne_2_bruitee_" + int2str(valeur_bruitage), "sinc_1_bruitee_" ...
	+ int2str(valeur_bruitage)];
liste_surface = ["gaussienne_1_pepper_bruitee_" ...
	+ int2str(valeur_bruitage), "gaussienne_2_bruitee_" + int2str(valeur_bruitage), "sinc_1_bruitee_" ...
	+ int2str(valeur_bruitage)];
%liste_surface = ["gaussienne_1_bruitee_" + int2str(valeur_bruitage)];
liste_surface = ["gaussienne_1", "gaussienne_1_pepper", "gaussienne_2", "sinc_1"];
%liste_surface = ["gaussienne_1_bis", "gaussienne_1_pepper_bis", "gaussienne_2_bis", "sinc_1_bis"];
liste_surface = ["calotte"];
liste_rayon_voisinage = [4];
nombre_iteration = 1;
%liste_ecart_type_I = [0:0.5:3];
liste_ecart_type_I = -1;
%liste_ecart_type_grad = [0:10];
liste_ecart_type_grad = 0;
filtrage = 0;
grille_pixels = 5;
liste_nombre_vues = [3];
liste_nombre_profondeur_iteration = [5000];
utilisation_profondeur_GT = 0;
utilisation_normale_GT = 0;
utilisation_mediane_normale = 1;

%% Variables
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
if (utilisation_profondeur_GT)
	fichier_profondeur_GT = "__profondeurs_GT";
	liste_ecart_type_I = 0;
	nb_ecart_type_I = 1;
else
	fichier_profondeur_GT = "";
end
if (utilisation_normale_GT)
	fichier_normale_GT = "__normales_GT";
	liste_ecart_type_grad = 0;
	nb_ecart_type_grad = 1;
else
	fichier_normale_GT = "";
end
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
							disp("MVS :");
							[z_estime_mvs,erreur_z_mvs,angles_mvs,normales_mvs] = mvs(premiere_iteration,surface,nombre_vues,rayon_voisinage,ecart_type_I,ecart_type_grad,nombre_z,z_estime_mvs,espace_z,utilisation_profondeur_GT,grille_pixels);
							disp("MVS modifié :");
							[z_estime_mvsm,erreur_z_mvsm,espace_z,normales_mvsm,erreur_angle_moy,erreur_angle_med] = mvs_modifie(premiere_iteration,surface,nombre_vues,rayon_voisinage,ecart_type_I,ecart_type_grad,nombre_z,z_estime_mvsm,espace_z,utilisation_profondeur_GT,utilisation_normale_GT,utilisation_mediane_normale,grille_pixels);
						end
						% Sauvegarde des résultats
						nom_fichier = "Surface_" + surface + "__nb_vues_" + int2str(nombre_vues) ...
							+ "__patch_" + int2str(taille_patch) + "x" + int2str(taille_patch) ...
							+ "__nb_profondeur_" + int2str(nombre_profondeur_iteration) ...
							+ fichier_bruite + fichier_profondeur_GT + fichier_normale_GT ...
							+ fichier_mediane + ".mat";
						path = "../../result/tests/";
						save(path+nom_fichier,"surface","nombre_vues","taille_patch","nombre_profondeur_iteration","z_estime_mvs","z_estime_mvsm","erreur_z_mvs","erreur_z_mvsm","normales_mvs","normales_mvsm","erreur_angle_moy","erreur_angle_med","grille_pixels");
						disp("Enregistrement sous : " + nom_fichier);
					end
					disp("------------------------------------");
				end
			end
		end
	end
end
