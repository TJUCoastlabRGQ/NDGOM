classdef SWENonhydroHLLNumFluxSolver1d < SWEHLLNumFluxSolver1d
    %SWENONHYDROHLLNUMFLUXSOLVER1D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    methods
        function [ fluxS ] = evaluate( obj, hmin, gra, nx, fm, fp, mesh, edge )
            fluxS = evaluate@SWEHLLNumFluxSolver1d( obj, hmin, gra, nx, fm, fp, mesh, edge );
            % upwind flux for the HW term
            tempfluxS = ( fm(:,:,2) .* fm(:,:,5) ./ fm(:,:,1) + ...
                fp(:,:,2) .* fp(:,:,5) ./ fp(:,:,1) ) .* nx ./ 2 ;
            temphum = fm(:,:,2);  temphm = fm(:,:,1);temphwm = fm(:,:,5);
            temphup = fp(:,:,2);  temphp = fp(:,:,1);temphwp = fp(:,:,5);
            Index = ( temphum .* nx  > 0 & - temphup .* nx  <= 0 );
            tempfluxS( Index ) =  ( temphum(Index) .* temphwm(Index) ./ temphm(Index) ) .* nx( Index ) ;
            Index = ( temphum .* nx  <= 0 & - temphup .* nx > 0 );
            tempfluxS( Index ) =   ( temphup(Index) .* temphwp(Index) ./ temphp(Index) ) .* nx( Index );
            fluxS(:,:,3) =  tempfluxS;    
            
            % upwind flux for the HW term
%             fluxS(:,:,3) = mxEvaluateUpwindNumFlux( int8(mesh.status), edge.FToE, fm(:,:,1), fm(:,:,2),...
%                  fm(:,:,5), fp(:,:,1), fp(:,:,2), fp(:,:,5), nx );            
        end        
        
    end
    
end

