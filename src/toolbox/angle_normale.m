% Calculer l'angle entre deux normales unitaires.
function angle = angle_normale(normale_1,normale_2)
	angle = abs((180/pi) * atan2(vecnorm(cross(normale_1,normale_2)),dot(normale_2,normale_2)));
end
