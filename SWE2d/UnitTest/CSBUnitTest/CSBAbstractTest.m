classdef CSBAbstractTest < SWEPreBlanaced2d
    
    methods

        function checkVolumeTerm( obj, m, k )

            mesh = obj.meshUnion(m);
            % volume integral
            fprintf('\nCheck Volume Flux Term\n');
            % obj.advectionSolver.evaluateAdvectionRHS( obj.fphys );
            [ E, G ] = obj.matEvaluateFlux( mesh, obj.fphys{m} );
            cell = mesh.cell;
            fprintf('-------- volume flux --------\n');
            fprintf('|  Field  |  Node Id  |  E  |  G  |  F  |\n');
            for fld = 1:3
                for n = 1:cell.Np
                    fprintf('|  fld = %d  |  %d  |  E = %e  |  G = %e  |  F = %e  |\n', fld, n, ...
                        E(n, k, fld), G(n, k, fld), obj.frhs{m}(n, k, fld)  );
                end
            end
            % volume integral
            obj.frhs{m} = zeros( obj.meshUnion.cell.Np, obj.meshUnion.K );
            obj.matEvaluateSourceTerm( obj.fphys );
            fprintf('-------- topography source --------\n');
            for fld = 1:3
                for n = 1:cell.Np
                    fprintf('|  fld = %d  | %d | S = %e  |\n', fld, n, obj.frhs{m}(n, k, fld) );
                end
            end
            
            fprintf('\nExact Volume Flux Term\n');
            Eh = ones( Np, 1 ); Ehu = ones( Np, 1 ); Ehv = ones( Np, 1 );
            Gh = ones( Np, 1 ); Ghu = ones( Np, 1 ); Ghv = ones( Np, 1 );
            for n = 1:Np
                [ Eh(n), Ehu(n), Ehv(n), Gh(n), Ghu(n), Ghv(n) ] = obj.evaluateVolumeFluxTerm...
                    ( h(n), hu(n), hv(n), z(n) );
            end
            
            for n = 1:Np
                fprintf('n = %d, E(h) = %f, E(hu) = %f, E(hv) = %f\n', ...
                    n, Eh(n), Ehu(n), Ehv(n) );
            end
            for n = 1:Np
                fprintf('n = %d, G(h) = %f, G(hu) = %f, G(hv) = %f\n', ...
                    n, Gh(n), Ghu(n), Ghv(n) );
            end
        end% func
        
        function checkSourceTerm( obj, m, k )
            mesh = obj.meshUnion(m);
            Np = mesh.cell.Np;
            
            h = obj.fphys{m}(:, k, 1); 
            z = obj.fphys{m}(:, k, 4);
            zx = obj.zGrad{m}(:, k, 1);
            zy = obj.zGrad{m}(:, k, 2);
            fprintf('\nConservative Variable\n');
            for n = 1:Np
                fprintf('n = %d, h = %f,  z = %f, zx = %f, zy = %f, eta = %f\n', ...
                    n, h(n), z(n), zx(n), zy(n), h(n) + z(n));
            end
            
            obj.matEvaluateTopographySourceTerm( obj.fphys );
            Shu = obj.frhs{m}(:, :, 2);
            Shv = obj.frhs{m}(:, :, 3);
            fprintf('Source Term\n');
            for n = 1:Np
                fprintf('n = %d, Shu = %f, Shv = %f\n', n, Shu(n, k), Shv(n, k));
            end
            
            fprintf('\nExact Source Term\n');
            
            Shu = zeros(Np, 1); Shv = zeros(Np, 1);
            for n = 1:Np
                [ Shu(n), Shv(n) ] = evaluateSourceTerm( obj, h(n), z(n), zx(n), zy(n) );
            end
            
            for n = 1:Np
                fprintf('n = %d, Shu = %f, Shv = %f\n', n, Shu(n), Shv(n));
            end
        end
    end
    
    methods( Access = protected )
        function [ Shu, Shv ] = evaluateSourceTerm( obj, h, b, bx, by )
            Shu = - obj.gra .* ( h +b ) .* bx;
            Shv = - obj.gra .* ( h +b ) .* by;
        end
        
        function [ Eh, Ehu, Ehv, Gh, Ghu, Ghv ] = evaluateVolumeFluxTerm( obj, h, hu, hv, z )
            Eh = hu;
            Gh = hv;
            Ehu = ( hu.^2 / h + 0.5 * obj.gra * (h ^ 2 - z ^ 2) );
            Ehv = hu * hv / h;
            Ghu = hu * hv / h;
            Ghv = ( hv.^2 / h + 0.5 * obj.gra * (h ^ 2 - z ^ 2) );
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 0.2;
            outputIntervalNum = 2;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('temporalDiscreteType') = NdgTemporalIntervalType.DeltaTime;
            option('obcType') = NdgBCType.None;
            option('outputIntervalType') = NdgIOIntervalType.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = NdgTemporalDiscreteType.RK45;
            option('integralType') = NdgDiscreteIntegralType.QuadratureFree;
            option('equationType') = NdgDiscreteEquationType.Strong;
            option('WellBlancedType') = true;
            option('limiterType') = NdgLimiterType.Vert;
        end
    end
end


