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
surface = "gaussienne_2";

%% Variables
path = "../../data/";
fichier_surface = "simulateur_" + surface + "_formate.mat";
load(path+fichier_surface);

%% Algorithme
% Préparation
X = 1:nombre_colonnes;
X = (X - u_0) / facteur_k;
Y = 1:nombre_lignes;
Y = (Y - v_0) / facteur_k;
[X,Y] = meshgrid(X,Y);

% Affichage surface
figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
sl = surfl(X,Y,-z(:,:,1),s);
sl.EdgeColor = 'none';
grid off;
colormap gray;
axis equal;

% Affichage
figure;
imshow(I(:,:,1));
