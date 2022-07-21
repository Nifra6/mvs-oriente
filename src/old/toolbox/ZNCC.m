function score = ZNCC(I_1,I_2)
	% Vectorisation des images
	I_1_vec = I_1(:);
	I_2_vec = I_2(:);
	% Calcul de l'espérance
	E_I1 = mean(I_1_vec);
	E_I2 = mean(I_2_vec);
	% Centrage des données
	I_1_centre = I_1_vec - E_I1;
	I_2_centre = I_2_vec - E_I2;
	% Calcul de la ZNCC
	numerateur = I_1_centre' * I_2_centre;
	%denominateur = norm(I_1_centre) * norm(I_2_centre);
	denominateur = std(I_1_vec) * std(I_2_vec);
	score = numerateur / denominateur;
end
