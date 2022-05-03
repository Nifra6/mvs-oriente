function [dx_V, dy_V] = gradient_correct(V,masque,facteur)
	addpath(genpath('../../../Développement/Ortho/Toolbox/'));
	[G,vMasque] = make_gradient(masque);
	Dx_main = G(1:2:end-1,:);
	Dy_main = G(2:2:end,:);
	clear G;
	dx_V_calcule = facteur * Dx_main * V(vMasque);
	dy_V_calcule = facteur * Dy_main * V(vMasque);
	dx_V = zeros(size(V));
	dy_V = zeros(size(V));
	dx_V(vMasque) = dx_V_calcule;
	dy_V(vMasque) = dy_V_calcule;
	clear dx_V_calcule dy_V_calcule vMasque;
end
