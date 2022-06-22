function normales = normales_medianes(p_estimes,q_estimes)
	% Données
	[nb_paire_image,nb_pixels] = size(p_estimes);
	normales = zeros(3,nb_pixels);

	% Calcul d'une normale médiane
	for pixel = 1:nb_pixels

		% Calcul des normales du pixel considéré
		liste_normales = [p_estimes(:,pixel)' ; q_estimes(:,pixel)' ; -ones(1,nb_paire_image)];
		liste_normales = liste_normales ./ vecnorm(liste_normales);
		coordones_zero = zeros(1,size(liste_normales,2));
		quiver3(coordones_zero,coordones_zero,coordones_zero,liste_normales(1,:),liste_normales(2,:),liste_normales(3,:));

		% Calcul des erreurs géodésiques
		matrice_angles = zeros(nb_paire_image);
		for image_temoin = 1:nb_paire_image-1
			normale_image_temoin = liste_normales(:,image_temoin);
			liste_normales_autres = liste_normales(:,image_temoin+1:end);
			liste_normale_temoin = repmat(normale_image_temoin,1,size(liste_normales_autres,2));
			angles_calcules = angle_normale(liste_normale_temoin,liste_normales_autres);
			matrice_angles(image_temoin+1:end,image_temoin) = angles_calcules;
			matrice_angles(image_temoin,image_temoin+1:end) = angles_calcules';
		end

		% Sélection de la normale médiane
		scores = sum(matrice_angles,2);
		[~,indice_normale_mediane] = min(scores);
		normales(:,pixel) = liste_normales(:,indice_normale_mediane);
	end
end
