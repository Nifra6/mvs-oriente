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
surface = "gaussienne_decentree";

%% Variables
path = "../../data/perspectif/";
fichier_surface = "simulateur_" + surface + "_formate.mat";
load(path+fichier_surface);

%% Algorithme
% Préparation
X = 1:nb_colonnes;
X = (X - u_0) / 200;
Y = 1:nb_lignes;
Y = (Y - v_0) / 200;
[X,Y] = meshgrid(X,Y);
Z = z(:,:,1);
ind_z_0 = find(Z == 0);
Z(ind_z_0) = nan;

% Affichage surface
figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
sl = surfl(X,Y,-Z,s);
sl.EdgeColor = 'none';
grid off;
colormap gray;
axis equal;

% Affichage
figure;
imshow(I(:,:,1));
