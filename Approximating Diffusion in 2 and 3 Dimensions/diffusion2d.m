function [t, U_soln] = diffusion2d( kappa, h, U_init, U_bndry, t_int, nt )
    ti = t_int(1);
    tf = t_int(2);
    dt = (tf - ti)/(nt - 1);
    t = linspace( ti, tf, nt );
 
    [nx, ny] = size( U_init );
 
    U_soln = zeros( nx, ny, nt );
    U_soln(:, :, 1) = U_init;
    
    r = kappa*dt/h^2;

    for it = 2:nt
        U_soln(:, :, it) = U_bndry( t(it), nx, ny );
        
        for ix = 1:nx
            for iy = 1:ny
                if U_soln(ix, iy, it) == -Inf
                    Utmp = U_soln(ix, iy, it - 1);
                    U_soln(ix, iy, it) = Utmp;
                    
                    for dxy = [-1 1 0 0; 0 0 -1 1]
                        dix = ix + dxy(1);
                        diy = iy + dxy(2);

                        if ~isnan( U_soln(dix, diy, it - 1) )
                            U_soln(ix, iy, it) = U_soln(ix, iy, it) + ...
                                r*( U_soln(dix, diy, it - 1) - Utmp );
                        end
                    end
                end
            end
        end
    end

end