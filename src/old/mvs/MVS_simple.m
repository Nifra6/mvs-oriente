%% Clear
clear;
close all;

%% Données
load ../../data/donnees_calotte;
[nombre_lignes, nombre_colonnes, nombres_images] = size(I);
[i_k, j_k]  = find(masque(:,:,1));
ind_1		= sub2ind([nombre_lignes nombre_colonnes], i_k, j_k);
ind			= ind_1;
P_k 		= zeros(3,size(i_k,1),nombres_images);
P_k(:,:,1) 	= [i_k - C_x, j_k - C_y, zeros(length(i_k), 1)].';

%% Paramètres
valeurs_z   = 60:.1:120;

%% Boucle de reconstruction
n		= length(valeurs_z);
erreurs	= 10*ones(length(i_k), n);

tic

for i = 1:n

	% Affichage de la progression des calculs
	if mod(i,round(n/25)) == 0
		disp("Progression à " + int2str(i/n*100) + "%");
	end

	% Sélection d'une profondeur
	valeur_z 	= valeurs_z(i);
	P_k(3,:,1) 	= valeur_z;

	% Changements de repère
	for k = 1:nombres_images-1
		P_k(:,:,k+1) = R(:,:,k) * P_k(:,:,1);
		i_k(:,k+1) = round(P_k(1,:,k+1) + C_x).';
		j_k(:,k+1) = round(P_k(2,:,k+1) + C_y).';
	end

	% Vérification des pixels hors images
	condition_image = ones(size(i_k(:,1)));
	for k = 1:nombres_images-1
		condition_image = condition_image & i_k(:,k+1) > 0 & i_k(:,k+1) <= size(masque,1) & j_k(:,k+1) > 0 & j_k(:,k+1) <= size(masque,2);
	end
	
	% Calcul des indices valides
	for k = 1:nombres_images-1
		i_k(:,k+1) = (ones(size(i_k,1),1) - condition_image) + condition_image .* i_k(:,k+1);
		j_k(:,k+1) = (ones(size(j_k,1),1) - condition_image) + condition_image .* j_k(:,k+1);
		ind(:,k+1) = sub2ind([nombre_lignes nombre_colonnes], i_k(:,k+1), j_k(:,k+1));
	end
	
	% Calcul de l'erreur
	erreur_k = zeros(size(ind,1), nombres_images-1);
	for k = 1:nombres_images-1
		erreur_k(:,k) = erreur_k(:,k) + condition_image .* (I(ind(:,1)) - I(ind(:,k+1) + k * nombre_lignes * nombre_colonnes));
	end
	erreurs(:,i) = (1 / (nombres_images - 1)) * sum(erreur_k.^2,2); % MSE
	%erreurs(:,i) = (1 / (nombres_images - 1)) * (1 - exp(-sum(erreur_k.^2,2)/0.2^2)); % E robuste
end

toc

%% Résultats
% Sélections des profondeurs avec l'erreur minimale
[~,indices_min] = min(erreurs,[],2);
z_in = transpose(valeurs_z(indices_min));
z = zeros(nombre_lignes, nombre_colonnes);
indice_masque  = find(masque(:,:,1));
z(indice_masque) = z_in;

% Mesures
disp("==============")
disp("Mesure relative de profondeur")
sum(abs(Z_1(indice_masque) - z(indice_masque)),'all') / size(z_in,1)
ecart_moyen = sum(Z_1(find(masque(:,:,1))) - z_in) / size(z_in,1);
disp("Mesure relative de forme")
sum(abs(Z_1(indice_masque) - (z(indice_masque) + ecart_moyen)),'all') / size(z_in,1)


% Affichage
figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
plot3(X,Y,z,'k.');
xlabel('$x$','Interpreter','Latex','FontSize',30);
ylabel('$y$','Interpreter','Latex','FontSize',30);
zlabel('$z$','Interpreter','Latex','FontSize',30);
axis equal;
rotate3d;
