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
ecart_type = 2;
nombre_profondeur_iteration = 1000;

%% Variables

%% Algorithme
for i = 1:nombre_iteration
	disp("---")
	disp("Itération " + int2str(i) + " / " + int2str(nombre_iteration));
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
	[z_estime_mvs,erreur_z_mvs,~,normales_mvs] = mvs(premiere_iteration,surface,rayon_voisinage,ecart_type,nombre_z,z_estime_mvs,espace_z);
	[z_estime_mvsm,erreur_z_mvsm,espace_z,normales_mvsm,~,~] = mvs_modifie(premiere_iteration,surface,rayon_voisinage,ecart_type,nombre_z,z_estime_mvsm,espace_z);
end

%% Analyse des résultats
normales_fronto = zeros(size(normales_mvsm));
normales_fronto(3,:) = -1;
angles_mvs = angle_normale(normales_fronto,normales_mvs);
angles_mvsm = angle_normale(normales_fronto,normales_mvsm);

zones_angles = 0:10:180;
nombre_zones = size(zones_angles,2) - 1;
erreurs_mvs = zeros(1,nombre_zones);
erreurs_mvsm = zeros(1,nombre_zones);
label_zones = [];
for k = 1:size(zones_angles,2)-1
	indices_mvs = find(zones_angles(k) <= angles_mvs & angles_mvs < zones_angles(k+1));
	indices_mvsm = find(zones_angles(k) <= angles_mvsm & angles_mvsm < zones_angles(k+1));
	erreurs_mvs(k) = transpose(median(erreur_z_mvs(indices_mvsm)));
	erreurs_mvsm(k) = transpose(median(erreur_z_mvsm(indices_mvsm)));
	%label_zones = [ label_zones , int2str(zones_angles(k)) + "-" + int2str(zones_angles(k+1)) ];
	label_zones = [ label_zones , zones_angles(k+1) ];
end

%% Affichage
% Histogramme
figure
bar(label_zones,[erreurs_mvs ; erreurs_mvsm])
legend('MVS','MVS modifié')
xlabel('Angles des normales avec la direction de la caméra de référence')
ylabel('Erreurs de profondeurs médianes')

% Préparation de la reconstruction
load('../../data/simulateur_formate.mat','nombre_lignes','nombre_colonnes','facteur_k','u_0','v_0','s');
X = 1:nombre_colonnes;
X = (X - u_0) / facteur_k;
Y = 1:nombre_lignes;
Y = (Y - v_0) / facteur_k;
[X,Y] = meshgrid(X,Y);

% Affichage de la reconstruction
figure('Name','Relief MVS','Position',[0,0,0.33*L,0.5*H]);
sl = surfl(X,Y,-z_estime_mvs,s);
sl.EdgeColor = 'none';
grid off;
colormap gray;
axis equal;

% Affichage de la reconstruction
figure('Name','Relief MVS modifié','Position',[0,0,0.33*L,0.5*H]);
sl = surfl(X,Y,-z_estime_mvsm,s);
sl.EdgeColor = 'none';
grid off;
colormap gray;
axis equal;
