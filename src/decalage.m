function images_decalees = decalage(image)
	[nombre_lignes, nombre_colonnes] = size(image);
	images_decalees = zeros(nombre_lignes, nombre_colonnes, 9);
	directions = ['S-'; 'SE'; 'E-'; 'NE'; 'N-'; 'NO'; 'O-'; 'SO'];
	images_decalees(:,:,1) = image;
	for k = 1:8
		images_decalees(:,:,k+1) = decalage_unique(image, directions(k,:));
	end
end

function image_decalee = decalage_unique(image, direction)
	image_decalee = zeros(size(image));
	switch direction
		case 'S-'
			image_decalee(2:end,:) = image(1:end-1,:);
		case 'E-'
			image_decalee(:,2:end) = image(:,1:end-1);
		case 'N-'
			image_decalee(1:end-1,:) = image(2:end,:);
		case 'O-'
			image_decalee(:,1:end-1) = image(:,2:end);
		case 'SE'
			image_decalee(2:end,2:end) = image(1:end-1,1:end-1);
		case 'NE'
			image_decalee(1:end-1,2:end) = image(2:end,1:end-1);
		case 'NO'
			image_decalee(1:end-1,1:end-1) = image(2:end,2:end);
		case 'SO'
			image_decalee(2:end,1:end-1) = image(1:end-1,2:end);
	end
end
