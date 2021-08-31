function [U_soln] = laplace3d( U )
% Error and Warning Checking
% ==========================
%
% Checking to see if the input is a 2D array. If it is a higher  
% or lower order array, it will output an error and a warning    
% message because the code can not use it to solve for the       
% unknowns.


if ~all( size( U ) ~= [2, 2] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument U is not a 2-dimensional array' ) );    

sz=size (U);
for i=1:(sz(1,2))
    if U(1,i)==-Inf
        throw( MException( 'MATLAB:invalid_argument', ...
        'There are -Inf values in the matrix in positions which are not to be solved' ) );
    end
    if U((sz(1,1)),i)==-Inf
        throw( MException( 'MATLAB:invalid_argument', ...
        'There are -Inf values in the matrix in positions which are not to be solved' ) );
    end
end
for j=1:(sz(1,1))
    if U(i,1)==-Inf
        throw( MException( 'MATLAB:invalid_argument', ...
        'There are -Inf values in the matrix in positions which are not to be solved' ) );
    end
    if U(i,sz(1,2))==-Inf
        throw( MException( 'MATLAB:invalid_argument', ...
        'There are -Inf values in the matrix in positions which are not to be solved' ) );
    end
end


% Initialization
% ==============
%
%   Creating a 2D matrix given the size of the input array that  
% can be used to map in the next step.

    [n_x, n_y, n_z] = size( U );
    U_soln = U;

% Mapping the unknown points to a unique number from 1 to m
% =========================================================
%
%   Mapping each of the “m” unsolved points to a unique number   
% from 1 to “m”. This allows us to know where the unknowns are as  
% well as the numbers around them to solve.

    u_to_w = zeros( n_x, n_y, n_z );
    w_to_u = zeros( 3, n_x * n_y * n_z );
    m = 0;

    for i = 1:n_x
        for j = 1:n_y
            for k = 1:n_z
                if U(i, j, k) == -Inf
                    m = m + 1;
                    u_to_w(i, j, k) = m;
                    w_to_u(:, m) = [i, j, k]';
                end
            end
        end
    end

    % Create the sparse system of linear equations
    M = spalloc( m, m, 50*m );
    b = zeros( m, 1 );
    

% Creating and solving a system of linear equations
% =================================================
%
%   Using the information from step 3, creating a appropriately  
% size matrix of zeros and vector of zeros. To update the matrix  
% and the vector, check the adjacent points to the ith unknown  
% point. If it is an insulated boundary, do nothing. If it is a  
% Dirichlet boundary condition, subtract 1 from the ith diagonal  
% entry of the matrix and subtract the value from the ith entry  
% of the vector. If there is another unknown at the position j,  
% subtract 1 from the ith diagonal entry of the matrix and add 1 
% to the (i, j)th entry of the matrix.


    for k = 1:m% Get the coordinates of the kth point
        c = w_to_u(:,k);
        % Determing the 6 adjacent points
        top_point = c + [-1 0 0]';
        bottom_point = c + [1 0 0]';
        left_point = c + [0 -1 0]';
        right_point = c + [0 1 0]';
        forward_point = c + [0 0 1]';
        backward_point = c + [0 0 -1]';
        
        %put the values of the adjacent points in an array so it can be
        %acccesed in the loop
        points = [top_point, bottom_point, left_point, right_point, forward_point, backward_point];
        for i = 1:6
            p_1 = points(1, i);
            p_2 = points(2, i);
            p_3 = points(3, i);
            
            
            % Checking to see if it will follow the Dirichlet condition
            if(isnan(U(p_1, p_2, p_3)) == 0 && isinf(U(p_1, p_2, p_3)) == 0) 
                %subtract 1 from the jth diagonal entry of M
                M(k, k) = M(k, k) - 1;
                
                %subtract the value from the jth entry of the vector b
                b(k, 1) = b(k, 1) - U(p_1, p_2, p_3);
            
            %Checking to see if is a point is an unknown
            elseif (U(p_1, p_2, p_3) == -Inf)
                %Subtract 1 from the ith diagonal entry of M
                %and Add 1 to the (i, j)th entry of M

                unknowns = u_to_w(p_1, p_2, p_3);
                M(k, k) = M(k, k) - 1;
                M(k, unknowns) = M(k, unknowns) + 1;
            else 
                continue  
            end 
        end
        %w_to_u
        %four_points
        % the point is an insluated boundary point, a Dirichlet    end
    end 
    
% Substituting the values back into the matrix U_out
% ===================================================
%
%   Solve for Uout by dividing the “b” vector by the “M” matrix.

    w = M \ b;
    %A = full(M);
    for k = 1:m
        % put the value of w into the solution matrix
        assign_1 = w_to_u(1, k);
        assign_2 = w_to_u(2, k);
        assign_3 = w_to_u(3, k);
        U_soln(assign_1, assign_2, assign_3) = w(k);
    end

end
