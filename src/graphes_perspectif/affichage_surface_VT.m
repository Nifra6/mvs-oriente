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
Y = 1:nb_lignes;
[X,Y] = meshgrid(X,Y);
Z = z(:,:,1);
ind_z_0 = find(Z == 0);
Z(ind_z_0) = nan;
p = repmat(Z(:)',3,1) .* (inv(K) * [X(:)' ; Y(:)' ; ones(1,nb_lignes*nb_colonnes)]);
P = R(:,:,1)' * p - R(:,:,1)' * t(:,1);
X = reshape(P(1,:),nb_lignes,nb_colonnes);
Y = reshape(P(2,:),nb_lignes,nb_colonnes);
Z = reshape(P(3,:),nb_lignes,nb_colonnes);

% Affichage surface
figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
sl = surf(X,Y,Z,I(:,:,1));
sl.EdgeColor = 'none';
grid off;
%colormap gray;
axis equal;

% Affichage
figure;
imshow(I(:,:,1));
