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
surface = "gaussienne_1";
rayon_voisinage = 1;
nombre_iteration = 3;
liste_ecart_type = transpose(0:2:20);
nombre_profondeur_iteration = 10;

%% Variables
nombre_test = size(liste_ecart_type,1);
erreurs_angulaires_moyennes = zeros(nombre_test,1);
erreurs_angulaires_medianes = zeros(nombre_test,1);

%% Algorithme
for k = 1:nombre_test
	for i = 1:nombre_iteration
		disp("---")
		disp("Écart type " + int2str(k) + " / " + int2str(nombre_test) + ", itération " + int2str(i) + " / " + int2str(nombre_iteration));
		% Paramétrage
		premiere_iteration = (i == 1);
		if (premiere_iteration)
			z_estime = 0;
			nombre_z = nombre_profondeur_iteration + 1;
			espace_z = 0;
		else
			nombre_z = 2 * nombre_profondeur_iteration + 1;
		end
		% Exécution
		[z_estime,~,espace_z,~,erreur_angle_moy,erreur_angle_med] = mvs_modifie(premiere_iteration,surface,rayon_voisinage,liste_ecart_type(k),nombre_z,z_estime,espace_z);
		% Sauvegarde
		erreurs_angulaires_moyennes(k) = erreur_angle_moy;
		erreurs_angulaires_medianes(k) = erreur_angle_med;
	end
end

%% Affichage
figure
plot(liste_ecart_type,erreurs_angulaires_moyennes,'b-','LineWidth',2);
hold on
plot(liste_ecart_type,erreurs_angulaires_medianes,'g-','LineWidth',2);
hold off
legend('Erreurs moyennes','Erreurs médianes')
xlabel('Écart type')
ylabel('Erreurs angulaires')
