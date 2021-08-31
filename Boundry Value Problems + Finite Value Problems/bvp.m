% bvp
% I am trying to implement a solution that is able  
% to solve a boundary-value problem by creating a matrix of values that are points      
% approximated for the function. 
%
% Parameters
% ==========
%    c = Vector of values that represent the coefficients of the ODE
%        x_int = Vector of values representing the boundaries of x
%    u_int = Vector of values representing the boundaries of u
%
%    g = Forcing function
%
%    n = Scalar value representing the number of values we want returned
%
% Return Values
% =============
%    x = The vector xout is a column vector of n x-values going from a to b that 
%    are equally spaced
%    u = The vector uout is is a column vector of n u-values that are            
%   approximations of the x-values

function [x_out, u_out] = bvp( c, g, x_bndry, u_bndry, n )

% Argument Checking

    if ~all( size( c ) == [3, 1] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument c is not a 3-dimensional column vector' ) );
    end

    if ~isscalar( n ) || ( n ~= round( n ) )  
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument n is not an integer' ) );
    end
    
    if ~isa( g, 'function_handle' )
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument g is not a function handle' ) );
    end

    if ~all( size( x_int ) == [2, 1] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument x_int is not a 2-dimensional column vector' ) );
    end

    if ~all( size( u_int ) == [2, 1] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument u_int is not a 2-dimensional column vector' ) );
    end

%Defining the h value based off the boundaries of x
h = (x_bndry(end) - x_bndry(1))/(n-1);

%Defined the equation for the subdiagonal
d_sub = 2*c(1) - h*c(2);

%Defined the equation for the center diagonal
d = (2)*(h^2)*c(3) - 4*c(1);

%Defined the equation for the super diagonal
d_super = 2*c(1) + h*c(2);

%Defined an array of x values and took out the first and last value
x_out = linspace(x_bndry(1), x_bndry(2), (n))';
x_out(1) = [];
x_out(end) = [];

%Defined the error term 
b = (2*(h^2).*g(x_out))';
b(1) = b(1) - d_sub*u_bndry(1);
b(end) = b(end) - d_super*u_bndry(end);

%Initialized a (n-2)x(n-2) vector filled with zeroes
Matrix = zeros((n-2), (n-2));

%took the transpose of b so it can be used with the other matrix
b = b';

%Created a vector of the d-values repeated the number of times it will
%appear on the diagonal
d_sub_vector = (d_sub.*ones(n-3, 1));
d_super_vector = (d_super.*ones(n-3, 1));
d_vector = (d.*ones(n-2, 1));

%placed the respective d-values in their diagonals
diag_sub = diag(d_sub_vector, -1);
diag_super = diag(d_super_vector, 1);
diag_d = diag(d_vector);

%Added the diagonal values to the initialized zeros matrix
Matrix = Matrix + diag_sub + diag_super + diag_d;

u_out = Matrix\b;

%plotting the points of the approximation
plot(x_out, u_out, 'or')
hold on
plot(x_bndry(end), u_bndry(end), 'or')
hold on 
plot(x_bndry(1), u_bndry(1), 'or')
hold on
title('sshakim')