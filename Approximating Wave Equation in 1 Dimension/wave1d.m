% Solving the wave equation 
%
% Parameters
% ==========
%    c = speed of the wave
%    x_int = x interval
%    n_x = number of values to split x interval into
%    t_int = t interval
%    n_t = number of values to split t interval into
%    u_init =  initial value function
%    du_init = intital rate of change of function.
%    u_bndry = boundary condition function
%
% Return Values
% =============
%    x_out = The vector x is the boundary of values 
%    t_out = The vector t is the boundaries for the matrix
%    U_out =  matrix giving ghe temperature values for each point in space
%    in a certain time

function [x_out, t_out, U_out] = wave1d( c, x_int,  n_x, t_int, n_t, u_init, du_init, u_bndry )

    % ARGUMENT CHECKING
    
    
    
    %initializing the variables for the boundaries:
    a = x_int(1);
    b = x_int(2);   
    
    %initializing the boundaries for 
    t_initial = t_int(1);
    t_final = t_int(2);
    
    
    dt = (t_final - t_initial)/(n_t - 1);
    
    %determining the h value
    h = (b - a)/(n_x - 1);
    
    %determining the value of r
    r = (c*dt/h)^2;
    
      %initializing the matrix of zeroes
    U_mat = zeros(n_x, n_t);
    
    
    % ERROR CHECKING
    % Calculate the error ratio (c * del_t / h)^2 to see if it is greater or equal to 1
        err = (c*dt/h)^2;
        
        if err>= 1
            
            noft = ceil((c*(t_final-t_initial)/h) + 1);
            
            throw( MException( 'MATLAB:invalid_argument', ...
                'The ratio (c * del_t / h)^2 = %f >= 1, consider using n_t = %d', ...
                err, noft ) );
        end
        
    %creating a t vector using n_t
    t_vec = linspace(t_initial, t_final, n_t);
    
    %creating an x vector using n_x
    x_1 = linspace(a, b, n_x)';
    x_vec = u_init(x_1);
    
    % Assign intial values using x vector
    U_mat(:,1) = x_vec;
    
    % putting the boundaries into the vector
    bound = u_bndry(t_vec);
    top_bound = bound(1);
    bot_bound = bound(2);
    
    
    % Assign boundary space values over time
    U_mat(1, 2:end) = top_bound;
    U_mat(end, 2:end) = bot_bound;
    
    % SOLVING
    % used the finite difference fomula to find every value in the matrix
    for t = 2:n_t
        
        % Special case for t = t_2
        % Use Euler's approximation.
        if t == 2
            u_x2 = U_mat(2:end - 1, 1) + dt*(du_init(x_vec(2:end - 1)));
            U_mat(2:end - 1, 2) = u_x2;
            
            
            
            if isnan(U_mat(1, 2))
                % Solve using forward and backward divided difference approximation
                U_mat(1, 2) = (4/3)*U_mat(2, 2) - (1/3)*U_mat(3, 2);
            end
            if isnan(U_mat(end, 2))
                % Solve using forward and backward divided difference approximation
                U_mat(end, 2) = (4/3)*U_mat(end - 1, 2) - (1/3)*U_mat(end - 2, 2);
            end
            
      
        else
            u_i_k = U_mat(:, t-1);
            u_i_k1 = 2*(U_mat(2:end - 1, t-1)) - (U_mat(2:end - 1, t-2)) + r.*(diff(u_i_k, 2));
            
            U_mat(2:end - 1, t) = u_i_k1;
            
            
            if isnan(U_mat(1, t))
                %using forward and backward divided difference approximation
                U_mat(1, t) = (4/3)*U_mat(2, t) - (1/3)*U_mat(3, t);
            end
            if isnan(U_mat(end, t))
                %using forward and backward divided difference approximation
                U_mat(end, t) = (4/3)*U_mat(end - 1, t) - (1/3)*U_mat(end - 2, t);
            end
        end
    end
  
    x_out = x_1;
    t_out = t_vec;
    U_out = U_mat;

end
