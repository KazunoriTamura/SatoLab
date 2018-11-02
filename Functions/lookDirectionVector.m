function [ L ] = lookDirectionVector( theta, phi )
% ‹ü•ûŒüƒxƒNƒgƒ‹

x = sin(theta) * cos(phi);
y = sin(theta) * sin(phi);
z = cos(theta);

L = [x,y,z].';

end

