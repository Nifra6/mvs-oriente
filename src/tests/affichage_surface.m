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
taille_patch = 9;
nombre_profondeur_iteration = 10000;
grille_pixel = 10;

%% Variables
path = "../../data/";
fichier_surface = "simulateur_" + surface + "_formate.mat";
load(path+fichier_surface,'nombre_lignes','nombre_colonnes','facteur_k','u_0','v_0','s');
path = "../../result/tests/";
nom_fichier = "Surface_" + surface + "__patch_" + int2str(taille_patch) + "x" ...
		  	+ int2str(taille_patch) + "__nb_profondeur_" + int2str(nombre_profondeur_iteration) + ".mat";
fichier_surface = "simulateur_" + surface + "_formate.mat";
load(path+nom_fichier);

%% Algorithme
% Préparation
X = 1:grille_pixel:nombre_colonnes;
X = (X - u_0) / facteur_k;
Y = 1:grille_pixel:nombre_lignes;
Y = (Y - v_0) / facteur_k;
[X,Y] = meshgrid(X,Y);

% Affichage
figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
sl = surfl(X,Y,-z_estime_mvsm(1:grille_pixel:end,1:grille_pixel:end),s);
sl.EdgeColor = 'none';
grid off;
colormap gray;
axis equal;
