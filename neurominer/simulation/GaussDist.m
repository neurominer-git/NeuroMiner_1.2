function output = GaussDist(xo,sigma,PA)
r = sqrt(-2.0.*(sigma^2).*log(rand(PA,1)));
phi = 2.0.*pi.*rand(PA,1);
output = xo+r.*cos(phi);