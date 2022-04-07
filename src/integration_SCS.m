function [U,rmse] = integration_SCS(p,q,gt)
% Ref : Direct Analytical Methods for Solving Poisson Equation in 
% Computer Vision problems - PAMI 1990
%
% example :
% p = ones(100,100); 
% q = ones(100,100); 
% u = simchony(p,q); 
%
% This performs the least square solution to \nabla u = [p,q], i.e. :
% min \int_\Omega \| \nablua U - [p,q] \|^2
% where \Omega is square and the natural Neumann boundary condition 
% \mu \cdotp (\nabla u -[p,q]) = 0 is used (2nd order derivatives 
% are neglected), where \mu is the outer 
% normal to \Omega. 
%
% % Axis : O->y
%          |
%          x
%
% The problem comes to solve a linear system Ax=b, where A is a bloc 
% Toplitz matrix. Fast solution is provided by Discrete Cosine Transform

% Compute div(p,q):
px = 0.5*(p([2:end end],:)-p([1 1:end-1],:));
qy = 0.5*(q(:,[2:end end])-q(:,[1 1:end-1]));

% Div(p,q) + Boundary Condition:
f = px+qy;
f(1,2:end-1) = 0.5*(p(1,2:end-1)+p(2,2:end-1));
f(end,2:end-1) = 0.5*(-p(end,2:end-1)-p(end-1,2:end-1));
f(2:end-1,1) = 0.5*(q(2:end-1,1)+q(2:end-1,2));
f(2:end-1,end) = 0.5*(-q(2:end-1,end)-q(2:end-1,end-1));

f(1,1) = 0.5*(p(1,1)+p(2,1)+q(1,1)+q(1,2));
f(end,1) = 0.5*(-p(end,1)-p(end-1,1)+q(end,1)+q(end,2));
f(1,end) = 0.5*(p(1,end)+p(2,end)-q(1,end)-q(1,end-1));
f(end,end) = 0.5*(-p(end,end)-p(end-1,end)-q(end,end)-q(end,end-1));

% Cosine Transform of f:
fsin = dct2(f);

% Denominator:
[x,y] = meshgrid(0:size(p,2)-1,0:size(p,1)-1);
denom = (2*cos(pi*x/(size(p,2)))-2) + (2*cos(pi*y/(size(p,1))) - 2);
Z = fsin./(denom);
Z(1,1) = 0.5*Z(1,2)+0.5*Z(2,1);		% Or whatever...

% Inverse Cosine Transform:
U = idct2(Z);

% Ground truth available:
if (nargin>2)
	moyenne_ecarts = mean(U(:)-gt(:));
	U = U-moyenne_ecarts;
	npix = size(p,1)*size(p,2);
	rmse = sqrt((sum((U(:)-gt(:)).^2))/npix);
else
	U = U-min(U(:));
end

end
