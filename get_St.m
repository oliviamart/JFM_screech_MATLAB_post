function [St] = get_St(f, M, h, b)

De = 2 * sqrt(h*b/pi);     % effective nozzle diameter
c = 1;                     % non-dimensionalized speed of sound
U = M*(1+0.2*M^2)^(-1/2);  % U/c_infinity
f = f*c/h;
St = f * De/U;

end
