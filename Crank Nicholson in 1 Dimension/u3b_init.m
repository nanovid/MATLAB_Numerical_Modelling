% The initial conditions for Laboratory Problem 3.8b.

function [u] = u3b_init(x)
	u = 4/pi*sin( x ) + 4/3/pi*sin( 3*x ) + 4/5/pi*( sin( 5*x ) );
end
