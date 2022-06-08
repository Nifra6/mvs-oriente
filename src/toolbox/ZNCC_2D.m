function score = ZNCC_2D(I_1,I_2)
	% Calcul de l'espérance
	E_I1 = mean(I_1,2);
	E_I2 = mean(I_2,2);
	% Centrage des données
	I_1_centre = I_1 - E_I1;
	I_2_centre = I_2 - E_I2;
	% Calcul de la ZNCC
	numerateur = dot(I_1_centre,I_2_centre,2);
	%denominateur = norm(I_1_centre) * norm(I_2_centre);
	denominateur = std(I_1,0,2) .* std(I_2,0,2);
	score = numerateur ./ denominateur;
end
