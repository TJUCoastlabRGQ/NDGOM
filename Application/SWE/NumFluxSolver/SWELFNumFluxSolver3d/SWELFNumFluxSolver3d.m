classdef SWELFNumFluxSolver3d
    %SWELFNUMFLUXSOLVER3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    methods
        function [ fluxS ] = evaluate( obj, hmin, gra, nx, ny, nz, fm, fp )
            %             fluxS = mxEvaluate( hmin, gra, nx, ny, fm, fp );
            fluxS = obj.LFNumFluxTerm( hmin, gra, fm, fp, nx, ny, nz );
        end
    end
    
    methods( Access = protected )
        
        %         function E = NodalFluxTerm( obj, gra, h, hu, hv )
        %             E(:, :, 1) = hu;
        %             E(:, :, 2) = 0.5 * gra * h.^2 + hu.^2 ./ h;
        %             E(:, :, 3) = hu .* hv ./ h;
        %         end
        
        function [ FluxS ] = LFNumFluxTerm( obj, hmin, gra, fm, fp, nx, ny, nz )
            
            FluxM(:, :, 1) = ( fm(:,:,1).^2./fm(:,:,4)  + 1/2 * gra * ( fm(:,:,4).^2 - obj.H0^2 ) ) .* nx +...
                ( fm(:,:,1) .* fm(:,:,2) ./ fm(:,:,4) ) .* ny + ( fm(:,:,1) .* fm(:,:,3) ./ fm(:,:,4) ) .* nz;
            %> $(hu^2+\frac{1}{2}g(H^2 - h^2) )* nx + huv * ny + u\omega * nz$
            FluxP(:, :, 1) = ( fp(:,:,1).^2./fp(:,:,4) + 1/2 * obj.gra * ( fp(:,:,4).^2 - obj.H0^2 )) .* nx +...
                ( fp(:,:,1) .* fp(:,:,2) ./ fp(:,:,4) ) .* ny+ ( fp(:,:,1) .* fp(:,:,3) ./ fp(:,:,4) ) .* nz;
            %> $huv * nx + ( hv^2+\frac{1}{2}g(H^2 - h^2) ) * ny + v\omega * nz$
            FluxM(:, :, 2) = ( fm(:,:,1) .* fm(:,:,2) ./ fm(:,:,4) ) .* nx +...
                ( fm(:,:,2).^2./fm(:,:,4) + 1/2 * obj.gra * ( fm(:,:,4).^2 - obj.H0^2 ) ) .* ny +...
                 ( fm(:,:,2) .* fm(:,:,3) ./ fm(:,:,4) ) .* nz;
            %> $huv * nx + ( hv^2+\frac{1}{2}g(H^2 - h^2) ) * ny + v\omega * nz$
            FluxP(:, :, 2) = ( fp(:,:,1) .* fp(:,:,2) ./ fp(:,:,4) ) .* nx +...
                ( fp(:,:,2).^2./fp(:,:,4) + 1/2 * obj.gra * ( fp(:,:,4).^2 - obj.H0^2 ) ) .* ny+ ...
                ( fp(:,:,2) .* fp(:,:,3) ./ fp(:,:,4) ) .* nz;
            
            %> $\lambda = abs( max( u*nx + v*ny + \frac{\omega}{H}*nz, ( u+\sqrt(gH) )*nx + ( v+\sqrt(gH) )*ny + \frac{\omega}{H}*nz,( u-\sqrt(gH) )*nx + ( v-\sqrt(gH) )*ny + \frac{\omega}{H}*nz )  )$
            lambda = abs( max( vertcat (max ( max ( fm(:, :, 1) ./ fm(:, :, 4) .* nx +  fm(:, :, 2) ./ fm(:, :, 4) .* ny + fm(:,:,3) ./ fm(:,:,4) .* nz,...
                - fp(:, :, 1) ./ fp(:, :, 4) .* nx -  fp(:, :, 2) ./ fp(:, :, 4) .* ny - fp(:,:,3) ./ fp(:,:,4) .* nz )), ...
                max( max (( fm(:, :, 1) ./ fm(:, :, 4) + sqrt(obj.gra* fm(:, :, 4))) .* nx...
                + ( fm(:, :, 2) ./ fm(:, :, 4) + sqrt(obj.gra* fm(:, :, 4))) .* ny + fm(:,:,3) ./ fm(:,:,4) .* nz ,...
                ( -fp(:, :, 1) ./ fp(:, :, 4) - sqrt(obj.gra* fp(:, :, 4))) .* nx...
                + ( -fp(:, :, 2) ./ fp(:, :, 4) - sqrt(obj.gra* fp(:, :, 4))) .* ny - fp(:,:,3) ./ fp(:,:,4) .* nz )),...
                max( max (( fm(:, :, 1) ./ fm(:, :, 4) - sqrt(obj.gra* fm(:, :, 4))) .* nx...
                + ( fm(:, :, 2) ./ fm(:, :, 4) - sqrt(obj.gra* fm(:, :, 4))) .* ny + fm(:,:,3) ./ fm(:,:,4) .* nz,...
                ( -fp(:, :, 1) ./ fp(:, :, 4) +  sqrt(obj.gra* fp(:, :, 4))) .* nx...
                + ( -fp(:, :, 2) ./ fp(:, :, 4) + sqrt(obj.gra* fp(:, :, 4))) .* ny - fp(:,:,3) ./ fp(:,:,4) .* nz )))));
            
            %> $\mathbf n\cdot\mathbf {F^*} = \mathbf n\cdot\frac{\mathbf{F^{(+)}}+\mathbf{F^{(-)}}}{2} - \frac{\lambda}{2}(\bold U^+ - \bold U^-)$
            FluxS(:, :, 1) = 0.5 .* ( FluxM(:, :, 1) + FluxP(:, :, 1) - ...
                bsxfun( @times, lambda , ( fp( :, :, 1 ) - fm( :, :, 1 ) ) ) );
            FluxS(:, :, 2) = 0.5 .* ( FluxM(:, :, 2) + FluxP(:, :, 2) - ...
                bsxfun( @times, lambda , ( fp( :, :, 2 ) - fm( :, :, 2 ) ) ) );
            
            
            %             % normal flux
            %             qnM = + fm(:, :, 2) .* nx + fm(:, :, 3) .* ny;
            %             qvM = - fm(:, :, 2) .* ny + fm(:, :, 3) .* nx;
            %             qnP = + fp(:, :, 2) .* nx + fp(:, :, 3) .* ny;
            %             qvP = - fp(:, :, 2) .* ny + fp(:, :, 3) .* nx;
            %
            %             % evaluate local eigenvalues of the flux Jacobian
            %             lamda = max( ...
            %                 sqrt( gra * fm(:, :, 1) ) + sqrt( qnM.^2 + qvM.^2 )./fm(:, :, 1), ...
            %                 sqrt( gra * fp(:, :, 1) ) + sqrt( qnP.^2 + qvP.^2 )./fp(:, :, 1) );
            %
            %             EM = obj.NodalFluxTerm( gra, fm(:, :, 1), qnM, qvM );
            %             EP = obj.NodalFluxTerm( gra, fp(:, :, 1), qnP, qvP );
            %
            %             FluxS(:, :, 1) = ( EM(:, :, 1) + EP(:, :, 1) ) .* 0.5 ...
            %                 - lamda .* ( fp(:, :, 1) - fm(:, :, 1) );
            %             Fqn = ( EM(:, :, 2) + EP(:, :, 2) ) .* 0.5 - lamda .* ( qnP - qnM );
            %             Fqv = ( EM(:, :, 3) + EP(:, :, 3) ) .* 0.5 - lamda .* ( qvP - qvM );
            %
            %             FluxS(:, :, 2) = (Fqn .* nx - Fqv .* ny);
            %             FluxS(:, :, 3) = (Fqn .* ny + Fqv .* nx);
        end
    end
    
end

