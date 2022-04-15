%% Clear
clear;
close all;

%% Données
albedo = 1;
load ../../data/donnees_calotte;
eclairage = S;
[nombre_lignes, nombre_colonnes, nombres_images] = size(I);
[i_k, j_k]  = find(masque(:,:,1));
ind_1		= sub2ind([nombre_lignes nombre_colonnes], i_k, j_k);
ind			= ind_1;
P_k 		= zeros(3,size(i_k,1),nombres_images);
P_k(:,:,1) 	= [i_k - C_x, j_k - C_y, zeros(length(i_k), 1)].';

%% Paramètres
valeurs_z   = 60:.1:120;
lambda      = 1/(nombres_images-1);

%% Calcul des gradients
dx_I_k = zeros(size(I));
dy_I_k = zeros(size(I));
for k = 1:nombres_images
	[dy_I, dx_I] = gradient(I(:,:,k));
	dx_I_k(:,:,k) = dx_I;
	dy_I_k(:,:,k) = dy_I;
end
grad_I_1 	= [ dx_I_1(ind_1) , dy_I_1(ind_1) ].';
grad_I_x	= [dx_I_1(ind_1)'];
grad_I_y    = [dy_I_1(ind_1)'];

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
	
	% Calcul des gradients
	for k = 1:nombres_images-1
		i_k(:,k+1) = (ones(size(i_k,1),1) - condition_image) + condition_image .* i_k(:,k+1);
		j_k(:,k+1) = (ones(size(j_k,1),1) - condition_image) + condition_image .* j_k(:,k+1);
		ind(:,k+1) = sub2ind([nombre_lignes nombre_colonnes], i_k(:,k+1), j_k(:,k+1));
		grad_I_x(k+1,:) = dx_I_k(ind(:,k+1) + k*(size(dx_I_k,1) * size(dx_I_k,2)))';
		grad_I_y(k+1,:) = dy_I_k(ind(:,k+1) + k*(size(dy_I_k,1) * size(dy_I_k,2)))';
	end
	
	% Calcul des numérateurs et dénominateurs
	denominateur = [];
	numerateur_1 = [];
	numerateur_2 = [];
	for k = 1:nombres_images-1
		denominateur(k,:) = R(1:2,3,k)' * [grad_I_x(k+1,:); grad_I_y(k+1,:)];
		numerateur = [grad_I_x(1,:); grad_I_y(1,:)] - R(1:2,1:2,k)' * [grad_I_x(k+1,:); grad_I_y(k+1,:)];
		numerateur_p(k,:) = numerateur(1,:);
		numerateur_q(k,:) = numerateur(2,:);
	end

	% Calcul des coefficients p et q
	p_estim = zeros(nombres_images-1, size(denominateur,2));	
	q_estim = zeros(nombres_images-1, size(denominateur,2));	
	for k = 1:nombres_images-1
		p_estim(k,:) = numerateur_p(k,:) ./ (denominateur(k,:));
		q_estim(k,:) = numerateur_q(k,:) ./ (denominateur(k,:));
	end
	p_estim = p_estim.';
	q_estim = q_estim.';
	%condition_image(isnan(sum(p_estim,1))) = 0;
	%condition_image(isnan(sum(q_estim,1))) = 0;

	% Calcul de l'erreur
	erreur_k = zeros(size(ind,1), nombres_images-1);
	for k = 1:nombres_images-1
		%erreur_k(:,k) = (ones(size(ind,1),1) - condition_image) .* erreur_k(:,k) + condition_image .* (I(ind(:,1)) + ([p_estim(:,k) q_estim(:,k) -1*ones(size(p_estim,1),1)] * eclairage) ./ sqrt(p_estim(:,k).^2 + q_estim(:,k).^2 + 1));
		erreur_k(:,k) = (ones(size(ind,1),1) - condition_image) .* erreur_k(:,k) + condition_image .* (I(ind(:,1)) + -1 ./ sqrt(p_estim(:,k).^2 + q_estim(:,k).^2 + 1));
	end
	%erreurs(:,i) = erreurs(:,i) + lambda * condition_image .* (I(ind(:,1)) - I(ind(:,2) + nombre_lignes * nombre_colonnes)).^2;
	erreurs(:,i) = (1 / (nombres_images - 1)) * sum(erreur_k.^2,2); % MSE
	%erreurs(:,i) = (1 / (nombres_images - 1)) * (1 - exp(-sum(erreur_k.^2,2)/0.2^2)); % E robuste
end

toc

sum(isnan(erreurs),'all')
sum(isinf(erreurs),'all')
size(erreurs)

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
