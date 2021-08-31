function [t, U_out] = wave2d( c, h, U_init, dU_init, U_bndry, t_int, n_t )


    if ~isa( dU_init, 'function_handle' )
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument u_bndry is not a function handle' ) );
    end

    if ~isa( U_init, 'function_handle' )
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument u_init is not a function handle' ) );
    end

    if ~isa( U_bndry, 'function_handle' )
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument u_bndry is not a function handle' ) );
    end

    if ~isscalar( c ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument c is not a scalar' ) );
    end

    if ~all( size( t_int ) == [1, 2] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument t_int is not a 2-dimensional row vector' ) );
    end

    if ~isscalar( n_t ) || ( n_t  ~= round( n_t ) )  
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument n_t is not an integer' ) );
    end

    if ~isscalar( h ) || ( h ~= round( h ) )  
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument h is not an integer' ) );
    end
    
    
   %Initialize the constant that will be used in the wave equation as well as the matrix to input   
   % the final values into. This will give a matrix of answers. 
   % Also initialized the boundaries of the matrix that will help to find 
   % the solution of the interior points
    ti = t_int(1);
    tf = t_int(2);
    dt = (tf - ti)/(n_t - 1);
    t = linspace( ti, tf, n_t );
 
    [nx, ny] = size( U_init );
 
    U_out = zeros( nx, ny, n_t );
    U_out(:, :, 1) = U_init;
    
    r = (c*dt/h)^2;
    
        
    %   Check to make sure the c*dt/h^2 value is less than ¼ for 
    % 2D. This ensures that there is only a 
    % small error. If it is not less than those values, we must alert 
    % the user to change their nt value and stop the function.
        err = (c*dt/h)^2;
        
        if err >= 0.25
            % calculate minimum integer nt, so that err_rat < 1
            say = ceil((2*c*(tf-t0)/h) + 1);
            
            throw( MException( 'MATLAB:invalid_argument', ...
                'The ratio (c * del_t / h)^2 = %f >= 1/4, consider using n_t = %d', ...
                err, say ) );
        end

    %   Took the value of the constant, the matrix with the initial  
    % boundaries and solved for the heat conduction and diffusion     
    % equation as well as the wave equation. Then we inputted the    
    % values for the respective point in the matrix.
    U_out(:, :, 2) = U_out(:, :, 1) + dt*dU_init;

    for it = 3:n_t
        U_out(:, :, it) = U_bndry( t(it), nx, ny );
        
        for ix = 1:nx
            for iy = 1:ny
                if U_out(ix, iy, it) == -Inf
                    Utmp = U_out(ix, iy, it - 1);
                    U_out(ix, iy, it) = 2*Utmp - U_out(ix, iy, it - 2);
                    
                    for dxy = [-1 1 0 0; 0 0 -1 1]
                        dix = ix + dxy(1);
                        diy = iy + dxy(2);

                        if ~isnan( U_out(dix, diy, it - 1) )
                            U_out(ix, iy, it) = U_out(ix, iy, it) + ...
                                r*( U_out(dix, diy, it - 1) - Utmp );
                        end
                    end
                end
            end
        end
    end

end