%% Clear
clear;
close all;

%% Données
load ../data/donnees_calotte;

%% Paramètres
valeurs_z   = 60:.1:120;
lambda      = 1;


[i_1, j_1]  = find(masque_1);
ind_1       = sub2ind(size(I_1), i_1, j_1);
I_1_in      = I_1(ind_1);
p_1         = [i_1 - C_x, j_1 - C_y, zeros(length(i_1), 1)].';

%% Calcul des gradients
[dy_I_1, dx_I_1] = gradient(I_1);
[dy_I_2, dx_I_2] = gradient(I_2);
[dy_I_3, dx_I_3] = gradient(I_3);

grad_I_1 = [ dx_I_1(ind_1) , dy_I_1(ind_1) ].';

%% Boucle de reconstruction
n       = length(valeurs_z);
erreurs = 10*ones(length(i_1), n);

for i = 1:n
	% Sélection d'une profondeur
	valeur_z = valeurs_z(i);
	p_1(3, :) = valeur_z;

	% Changement de repère image 2 
	p_2 = R_2.' * p_1;
	i_2 = round(p_2(1, :) + C_x).';
	j_2 = round(p_2(2, :) + C_y).';

	% Changement de repère image 3
	p_3 = R_3.' * p_1;
	i_3 = round(p_3(1, :) + C_x).';
	j_3 = round(p_3(2, :) + C_y).';

	% Vérification des pixels hors images
	condition_image = i_2 > 0 & i_2 <= size(masque_2,1) & i_3 > 0 & i_3 <= size(masque_3,1) & j_2 > 0 & j_2 <= size(masque_2,2) & j_3 > 0 & j_3 <= size(masque_3,2); 

	% Calcul du gradient dans l'image 2
	i_2 = (ones(size(i_2,1),1) - condition_image) + condition_image .* i_2;
	j_2 = (ones(size(j_2,1),1) - condition_image) + condition_image .* j_2;
	ind_2 = sub2ind(size(I_2), i_2, j_2);
	grad_I_2 = [ dx_I_2(ind_2) , dy_I_2(ind_2) ].';

	% Calcul du gradient dans l'image 3
	i_3 = (ones(size(i_3,1),1) - condition_image) + condition_image .* i_3;
	j_3 = (ones(size(j_3,1),1) - condition_image) + condition_image .* j_3;

	ind_3 = sub2ind(size(I_3), i_3, j_3);
	grad_I_3 = [ dx_I_3(ind_3) , dy_I_3(ind_3) ].';

	% Calcul des numérateurs et dénominateurs
	a_1 = R_2(3,1) * grad_I_2(1,:) + R_2(3,2) * grad_I_2(2,:);
	b_1 = grad_I_1 - R_2(1:2,1:2) * grad_I_2;
	a_2 = R_3(3,1) * grad_I_3(1,:) + R_3(3,2) * grad_I_3(2,:);
	b_2 = grad_I_1 - R_3(1:2,1:2) * grad_I_3;

	% Calcul des coefficients p et q
	p_q = (a_1 .* b_1 + a_2 .* b_2) ./ (a_1.^2 + a_2.^2);
	p_estim = p_q(1, :).';
	q_estim = p_q(2, :).';

	% Calcul de l'erreur
	erreurs(:, i) = (ones(size(condition_image,1),1) - condition_image).* erreurs(:,i) + ... 
	condition_image .* (I_1(ind_1) - 1 ./ sqrt(p_estim.^2 + q_estim.^2 + 1)).^2 + ...
		(I_2(ind_2) - 1 ./ sqrt(p_estim.^2 + q_estim.^2 + 1)).^2 + ...
		(I_3(ind_3) - 1 ./ sqrt(p_estim.^2 + q_estim.^2 + 1)).^2 + ...
		lambda * (I_1(ind_1) - I_2(ind_2)).^2;

end

%% Résultats
% Sélections des profondeurs avec l'erreur minimale
[~,indices_min] = min(erreurs,[],2);
z_in = transpose(valeurs_z(indices_min));
z = zeros(256, 256);
z(find(masque_1)) = z_in;

% Affichage
figure('Name','Relief','Position',[0,0,0.33*L,0.5*H]);
plot3(X,Y,z,'k.');
xlabel('$x$','Interpreter','Latex','FontSize',30);
ylabel('$y$','Interpreter','Latex','FontSize',30);
zlabel('$z$','Interpreter','Latex','FontSize',30);
axis equal;
rotate3d;
