% The analytic solution for Laboratory Problem 3.8b.

function [u] = u3b(x, t)
	u = 4/pi*exp( -t ).*sin( x ) + 4/3/pi*exp( -9*t ).*sin( 3*x ) + 4/5/pi*exp( -25*t ).*( sin( 5*x ) );
end
