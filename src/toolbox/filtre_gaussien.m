function I_filtre = filtre_gaussien(I,taille_filtre)
	sigma_filtre = taille_filtre / 4;
	filtre_I = fspecial('gauss',taille_filtre,sigma_filtre);
	filtre_I = filtre_I / sum(filtre_I(:));
	I_filtre = zeros(size(I));
	for k = 1:size(I,3)
		I_filtre(:,:,k) = conv2(I(:,:,k),filtre_I,'same');
	end
