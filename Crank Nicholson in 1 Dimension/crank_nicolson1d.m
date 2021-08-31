% Using the heat-conduction/diffusion equation, we are able to determine the temperature at a particular time for a specific point in
%space. This will allow us to know the temperature and how it changes over time until an equilibrium is obtained. 
%Parameters
% ==========
%    kappa = The diffusivity coefficient 
%    x_rng = The bounds of the space range that is being considered
%    t_rng = The range of time that is being considered
%
%    u_init = The function handle that is giving the initial state 
%    u_bndry = The function handle giving the boundary conditions
%
%    nx = The number of intervals the x-range should be divided by
%    nt = The number of intervals the time should be divided into 
%
% Return Values
% =============
%    x_out = The vector x is the boundary of values for the U_out
%    t_out = The vector t is the top and bottom boundaries for the matrix
%    U_out =  The matrix U is the corresponding temperature values at a given time for each point in that space

function [x_out, t_out, U_out] = crank_nicolson1d1( kappa, x_rng, nx, t_rng, nt, u_init, u_bndry )




%Argument checking

% if ~isscalar( kappa ) 
%         throw( MException( 'MATLAB:invalid_argument', ...
%         'the argument kappa is not a scalar' ) );
%     end
% 
% if ~all( size( x_rng ) == [1, 2] ) 
%         throw( MException( 'MATLAB:invalid_argument', ...
%         'the argument x_rng is not a 2-dimensional row vector' ) );
%     end
%   
%     
% 
% if ~isscalar( nx ) || ( nx ~= round( nx ) )  
%     throw( MException( 'MATLAB:invalid_argument', ...
%     'the argument nx is not an integer' ) );
% end
% 
% if ~all( size( t_rng ) == [1, 2] ) 
%         throw( MException( 'MATLAB:invalid_argument', ...
%         'the argument t_rng is not a 2-dimensional row vector' ) );
% end
% 
% if ~isscalar( tx ) || ( tx ~= round( tx ) )  
%     throw( MException( 'MATLAB:invalid_argument', ...
%     'the argument tx is not an integer' ) );
% end
% 
% if ~isa( u_init, 'function_handle' )
%         throw( MException( 'MATLAB:invalid_argument', ...
%         'the argument u_init is not a function handle' ) );
% end
% if ~isa( u_bndry, 'function_handle' )
%     throw( MException( 'MATLAB:invalid_argument', ...
%     'the argument u_bndry is not a function handle' ) );
% end


%Determining the number of discrete points we want for a certain space
%interval
h = (x_rng(2) - x_rng(1))/(nx - 1);

%Determining the change in time for a specific interval 
dt = (t_rng(2)- t_rng(1))/(nt-1);

%Determining the constant for the heat/diffusion equation
r = (kappa*dt)/h^2
%   Check to make sure the kappa*dt/h^2 value is less than 0.5 to ensure small error. If it is not less than 0.5, we must alert the user to change their nt value and stop the function.
% if const <0.5
    
%Determining the x-values depending on the boundaries
x_out = linspace(x_rng(1), x_rng(2), nx)'

%Determining the t-values depending on the boundaries
t_out = linspace(t_rng(1), t_rng(2), nt);

%Creating a matrix of zeros that is going to be used for the values of the
%resulting heat/diffusion equation
U_out = zeros(nx, nt);

%Intializing the values for the first column 
col_1 = u_init(x_out);

%Assigning the values for the first column of the resulting matrix
U_out(:,1) = col_1;

%Initializing the top and bottom boundaries of the matrix
boundaries= u_bndry(t_out(1, 2:end));
a_bndry = boundaries(1,:);
b_bndry = boundaries(2,:);

% Putting the boundaries into the matrix
U_out(1,2:end) = a_bndry;
U_out(end,2:end) = b_bndry;
a = isnan(a_bndry(1));
b = isnan(b_bndry(1));
if a == 1
    
    %Creating a matrix of zeros for the crank nicolson method:
    A_crank = zeros(nx-2, nx-2);
    
    %Defining the values for the middle diagonal, the sub diagonal and the
    %super diagonal
    top_and_bottom = -1.*r;
    middle = 2.*(1+r);
    
    
    %Creating matrix for the diagonal
    crank_middle = ones(1, nx-2);
    crank_middle = crank_middle.*middle;
    crank_middle(1) = 2+(2/3)*r;
    crank_middle = diag(crank_middle, 0);
    
    %Creating matrix for the super diagonal
    crank_super = ones(1, nx-3);
    crank_super = crank_super.*top_and_bottom;
    crank_super(1) = (-2/3)*r;
    
    crank_super = diag(crank_super, 1);
    
    %Creating a matriz for the sub diagonal
    crank_sub = ones(1, nx-3);
    crank_sub = crank_sub.*top_and_bottom;
    crank_sub = diag(crank_sub, -1);

    % inserting the values to the crank nicolson matrix
    A_crank = A_crank + crank_middle + crank_super + crank_sub;


    %Calculating and assigning values for the solution matrix
    for x = 1:nt - 1
        if(x < 2)
            difference = diff(U_out(:, x), 2);
            ui = U_out(2:end-1, x);
            z = (ui + difference.*r);

            U_out(2:end-1, x+1) = z;
            U_out(1, x+1) = (4/3).*U_out(2, x+1) - (1/3).*U_out(3, x+1);
        end
        if(x >= 2)
            z = U_out(2:end-1, x);
            y = U_out(1:end, x); 
            crank_known_1 = 2.*z;
            crank_known_2 = r.*diff(y, 2);
            crank_known_3 = zeros(nx-2, 1);
%             crank_known_3(1) = r.*U_out(1, x);
            crank_known_3(end) = r.*U_out(end, x);

            crank_known = crank_known_1 + crank_known_2 + crank_known_3;
            A_crank;

            new_z = A_crank\crank_known;
            U_out(2:end-1, x+1) =  new_z;
            U_out(1, x+1) = (4/3).*U_out(2, x+1) - (1/3).*U_out(3, x+1);

        end 
    end
end

%If the bar is insulated
if b == 1
    
    %Creating a matrix of zeros for the crank nicolson method
    A_crank = zeros(nx-2, nx-2);
    
    %Defining the values for the middle diagonal, the sub diagonal and the
    %super diagonal
    
    top_and_bottom = -1.*r;
    middle = 2.*(1+r);
    
     
    %Creating matrix for the diagonal
    crank_middle = ones(1, nx-2);
    crank_middle = crank_middle.*middle;
    crank_middle(end) = 2+(2/3)*r;
    crank_middle = diag(crank_middle, 0);

    %Creating matrix for the super diagonal
    crank_super = ones(1, nx-3);
    crank_super = crank_super.*top_and_bottom;
    
    crank_super = diag(crank_super, 1);

    %Creating matrix for the sub diagonal
    crank_sub = ones(1, nx-3);
    crank_sub = crank_sub.*top_and_bottom;
    crank_sub(end) = (-2/3)*r;
    crank_sub = diag(crank_sub, -1);

    % inserting the values to the crank nicolson matrix
    A_crank = A_crank + crank_middle + crank_super + crank_sub;


    %Calculating and assigning values for the solution matrix
    for x = 1:nt - 1
        if(x < 2)
            difference = diff(U_out(:, x), 2);
            ui = U_out(2:end-1, x);
            z = (ui + difference.*r);

            U_out(2:end-1, x+1) = z;
            U_out(end, x+1) = (4/3).*U_out(end-1, x+1) - (1/3).*U_out(end-2, x+1);
        end
        if(x >= 2)
            z = U_out(2:end-1, x);
            y = U_out(1:end, x); 
            crank_known_1 = 2.*z;
            crank_known_2 = r.*diff(y, 2);
            crank_known_3 = zeros(nx-2, 1);
            crank_known_3(1) = r.*U_out(1, x);
%             crank_known_3(end) = r.*U_out(end, x);

            crank_known = crank_known_1 + crank_known_2 + crank_known_3;
%             crank_matrix_A;

            new_z = A_crank\crank_known;
            U_out(2:end-1, x+1) =  new_z;
            U_out(end, x+1) = (4/3).*U_out(end-1, x+1) - (1/3).*U_out(end-2, x+1);

        end
end
end 

if (a == 0 && b == 0) 
%Looping through all the points in the matrix with a zero value and
%determining the value using the diffusion equation 

    A_crank = zeros(nx-2, nx-2);
    middle = 2.*(1+r);
    top_and_bottom = -1.*r;
    crank_middle = ones(1, nx-2);
    crank_middle = crank_middle.*middle;
    crank_middle = diag(crank_middle, 0);

    crank_super = ones(1, nx-3);
    crank_super = crank_super.*top_and_bottom;
    crank_super = diag(crank_super, 1);

    crank_sub = ones(1, nx-3);
    crank_sub = crank_sub.*top_and_bottom;
    crank_sub = diag(crank_sub, -1);


    A_crank = A_crank + crank_middle + crank_super + crank_sub;



    for x = 1:nt - 1
        if(x < 2)
            difference = diff(U_out(:, x), 2);
            ui = U_out(2:end-1, x);
            z = (ui + difference.*r);

            U_out(2:end-1, x+1) = z;
        end
        if(x >= 2)
            z = U_out(2:end-1, x);
            y = U_out(1:end, x); 
            crank_known_1 = 2.*z;
            crank_known_2 = r.*diff(y, 2);
            crank_known_3 = zeros(nx-2, 1);
            crank_known_3(1) = r.*U_out(1, x);
            crank_known_3(end) = r.*U_out(end, x);

            crank_known = crank_known_1 + crank_known_2 + crank_known_3;
            A_crank;

            new_z = A_crank\crank_known;
            U_out(2:end-1, x+1) =  new_z;

        end 

end
end
U_out;
t_out;
mesh( t_out, x_out, U_out )
title ('m3muir and sshakim')
% else
%     rec_nt=ceil((2*(kappa*(t_rng(2)- t_rng(1)))/(h^2))+0.5);
%     disp("Your value of nt is too small for your inputted parameters. Choose a value greater than ")
%     disp(rec_nt)
% end

    
%Plotting the solution