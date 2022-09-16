% Calculer l'angle entre deux normales unitaires.
function angle = angle_normale(normale_1,normale_2)
	if ~isreal(normale_1)
		normale_1 = real(normale_1);
	end
	if ~isreal(normale_2)
		normale_2 = real(normale_2);
	end
	angle = abs((180/pi) * atan2(vecnorm(cross(normale_1,normale_2)),dot(normale_1,normale_2)));
end
