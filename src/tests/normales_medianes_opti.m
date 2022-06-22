function normales = normales_medianes_opti(p_estimes,q_estimes)
	% Mise en forme des données
	[nb_paire_image,nb_pixels] = size(p_estimes);
	normales = zeros(3,nb_pixels);
	p = reshape(p_estimes,[1 nb_paire_image nb_pixels]);
	q = reshape(q_estimes,[1 nb_paire_image nb_pixels]);

	% Calcul des normales
	liste_normales = [p ; q ; -ones(1,nb_paire_image, nb_pixels)];
	liste_normales = liste_normales ./ vecnorm(liste_normales);

	% Calcul des erreurs géodésiques
	matrice_angles = zeros(nb_paire_image,nb_paire_image,nb_pixels);
	for image_temoin = 1:nb_paire_image-1
		normale_image_temoin = liste_normales(:,image_temoin,:);
		liste_normales_autres = liste_normales(:,image_temoin+1:end,:);
		liste_normale_temoin = repmat(normale_image_temoin,1,size(liste_normales_autres,2),1);
		angles_calcules = angle_normale(liste_normale_temoin,liste_normales_autres);
		matrice_angles(image_temoin+1:end,image_temoin,:) = angles_calcules;
		matrice_angles(image_temoin,image_temoin+1:end,:) = permute(angles_calcules,[2 1 3]);
	end
	
	% Sélection de la normale médiane
	scores = sum(matrice_angles,2);
	[~,indice_normale_mediane] = min(scores,[],1);
	indice_normale_mediane = reshape(indice_normale_mediane,[1,nb_pixels]);	
	for k = 1:3
		ind = sub2ind([3, nb_paire_image, nb_pixels], k*ones(1,nb_pixels),indice_normale_mediane,1:nb_pixels);
		normales(k,:) = liste_normales(ind);
	end
end
