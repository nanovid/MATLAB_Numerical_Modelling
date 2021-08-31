function [t, U_soln] = diffusion3d( kappa, h, U_init, U_bndry, t_int, nt )
    ti = t_int(1);
    tf = t_int(2);
    dt = (tf - ti)/(nt - 1);
    t = linspace( ti, tf, nt );
 
    [nx, ny, nz] = size( U_init );
 
    U_soln = zeros( nx, ny, nz, nt );
    U_soln(:, :, :, 1) = U_init;
    
    r = kappa*dt/h^2;
    
        
    %
    % Calculate (c * del_t / h)^2 to ensure that it is the proper value
    % if it is >= 0.25
    % Provide a suitable value of n_t 
        err_rat = (c*dt/h)^2;
        
        if err_rat>= 0.25
            % calculate minimum integer nt, so that err_rat < 1
            noft = ceil((2*c*(tf-t0)/h) + 1);
            
            throw( MException( 'MATLAB:invalid_argument', ...
                'The ratio (c * del_t / h)^2 = %f >= 0.25, consider using n_t = %d', ...
                err_rat, noft ) );
        end

    for it = 2:nt
        U_soln(:, :, :, it) = U_bndry( t(it), nx, ny, nz );
        
        for ix = 1:nx
            for iy = 1:ny
                for iz = 1:nz
                    if U_soln(ix, iy, iz, it) == -Inf
                        Utmp = U_soln(ix, iy, iz, it - 1);
                        U_soln(ix, iy, iz, it) = Utmp;

                        for dxyz = [-1 1 0 0 0 0; 0 0 -1 1 0 0; 0 0 0 0 -1 1]
                            dix = ix + dxyz(1);
                            diy = iy + dxyz(2);
                            diz = iz + dxyz(3);

                            if ~isnan( U_soln(dix, diy, diz, it - 1) )
                                U_soln(ix, iy, iz, it) = U_soln(ix, iy, iz, it) + ...
                                    r*( U_soln(dix, diy, diz, it - 1) - Utmp );
                        
                            end
                        end
                    end
                end
            end
        end
    end

end