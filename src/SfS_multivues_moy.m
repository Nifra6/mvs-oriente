%% Clear
clear;
close all;

%% Données
load ../data/donnees_calotte;
[nombre_lignes, nombre_colonnes, nombres_images] = size(I);
[i_k, j_k]  = find(masque(:,:,1));
ind_1		= sub2ind([nombre_lignes nombre_colonnes], i_k, j_k);
ind			= ind_1;
P_k 		= zeros(3,size(i_k,1),nombres_images);
P_k(:,:,1) 	= [i_k - u_0, j_k - v_0, zeros(length(i_k), 1)].';

%% Paramètres
valeurs_z   	= 60:.1:120;
lambda      	= 1/(nombres_images-1);
interpolation 	= 'cubic';
ertimateur		= 'MSE'

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
		P_k(:,:,k+1) = R(:,:,k).' * P_k(:,:,1);
		i_k(:,k+1) = (P_k(1,:,k+1) + u_0).';
		j_k(:,k+1) = (P_k(2,:,k+1) + v_0).';
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
		grad_I_x(k+1,:) = interp2(dx_I_k(:,:,k+1),j_k(:,k+1),i_k(:,k+1),interpolation)';
		grad_I_y(k+1,:) = interp2(dy_I_k(:,:,k+1),j_k(:,k+1),i_k(:,k+1),interpolation)';
		i_k(:,k+1) = round(i_k(:,k+1));
		j_k(:,k+1) = round(j_k(:,k+1));
		ind(:,k+1) = sub2ind([nombre_lignes nombre_colonnes], i_k(:,k+1), j_k(:,k+1));
	end

	% Calcul des numérateurs et dénominateurs
	A 	= [];
	B_1 = [];
	B_2 = [];
	for k = 1:nombres_images-1
		A(k,:) = R(3,1,k) * grad_I_x(k+1,:) + R(3,2,k) * grad_I_y(k+1,:);
		b = [grad_I_x(1,:); grad_I_y(1,:)] - R(1:2,1:2,k) * [grad_I_x(k+1,:); grad_I_y(k+1,:)];
		B_1(k,:) = b(1,:);
		B_2(k,:) = b(2,:);
	end

	% Calcul des coefficients p et q
	p_q = 0;	
	for k = 1:nombres_images-1
		p_q = p_q + A(k,:) .* [B_1(k,:); B_2(k,:)];
	end
	p_q 	= p_q ./ sum(A.^2, 1);
	p_estim = p_q(1, :).';
	q_estim = p_q(2, :).';

	% Calcul de l'erreur
	erreur_k = zeros(size(ind,1), nombres_images);
	for k = 1:nombres_images
		erreur_k(:,k) = (ones(size(ind,1),1) - condition_image) .* erreur_k(:,k) + condition_image .* (I(ind(:,k) + (k-1) * nombre_lignes * nombre_colonnes) + -1 ./ sqrt(p_estim.^2 + q_estim.^2 + 1));
	end
	%erreurs(:,i) = erreurs(:,i) + lambda * condition_image .* (I(ind(:,1)) - I(ind(:,2) + nombre_lignes * nombre_colonnes)).^2;
	switch (estimateur)
		case 'MSE'
			erreurs(:,i) = (1 / nombres_images) * sum(erreur_k.^2,2);
		case 'Robuste'
			erreurs(:,i) = (1 / nombres_images) * (1 - exp(-sum(erreur_k.^2,2)/0.2^2));
	end
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
