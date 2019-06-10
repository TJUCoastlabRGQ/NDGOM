classdef SWENonhydroHLLNumFluxSolver2d < SWEHLLNumFluxSolver2d
    %SWENONHYDROHLLNUMFLUXSOLVER2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    methods
        function [ fluxS ] = evaluate( obj, hmin, gra, nx, ny, fm, fp )
            fluxS = evaluate@SWEHLLNumFluxSolver2d( obj, hmin, gra, nx, ny, fm, fp );
            % upwind flux for the HW term
            tempfluxS = ( fm(:,:,2) .* fm(:,:,6) ./ fm(:,:,1) + ...
                fp(:,:,2) .* fp(:,:,6) ./ fp(:,:,1) ) .* nx ./ 2 + ...
                ( fm(:,:,3) .* fm(:,:,6) ./ fm(:,:,1) + ...
                fp(:,:,3) .* fp(:,:,6) ./ fp(:,:,1) ) .* ny ./ 2;
            temphum = fm(:,:,2); temphvm = fm(:,:,3); temphm = fm(:,:,1);temphwm = fm(:,:,6);
            temphup = fp(:,:,2); temphvp = fp(:,:,3); temphp = fp(:,:,1);temphwp = fp(:,:,6);
            Index = ( temphum .* nx + temphvm .* ny > 0 & - temphup .* nx - temphvp .* ny <= 0 );
            tempfluxS( Index ) =  ( temphum(Index) .* temphwm(Index) ./ temphm(Index) ) .* nx( Index ) + ...
            ( temphvm(Index) .* temphwm(Index) ./ temphm(Index) ) .* ny( Index ) ;
            Index = ( temphum .* nx + temphvm .* ny <= 0 & - temphup .* nx - temphvp .* ny > 0 );
            tempfluxS( Index ) =   ( temphup(Index) .* temphwp(Index) ./ temphp(Index) ) .* nx( Index ) +...
            ( temphvp(Index) .* temphwp(Index) ./ temphp(Index) ) .* ny( Index )  ;
            fluxS(:,:,4) =  tempfluxS;            
        end
    end
    
end

