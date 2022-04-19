function XYZ = shapeFromDmOrtho(dm, imMask)

imask1 = find(imMask > 0);
imask2 = find(dm > 0);

imask = intersect(imask1,imask2);

maskIm = zeros(size(imMask));
maskIm(imask) = 1;

[nrows, ncols] = size(dm);

XYZ = NaN*ones(nrows*ncols,3);

[xx,yy] = meshgrid(1:ncols,1:nrows);	

xx = xx(imask);
yy = yy(imask);

XYZ(imask,1) = xx;
XYZ(imask,2) = yy;
XYZ(imask,3) = dm(imask);

XYZ = reshape(XYZ,[nrows,ncols,3]);
