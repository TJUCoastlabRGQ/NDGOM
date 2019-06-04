classdef SWEPreBlanaced2d < SWEConventional2d
    
    methods( Hidden )
        function [ E, G ] = matEvaluateFlux( obj, mesh, fphys )
              [ E, G ] = obj.volumefluxSolver.evaluate( obj.hmin, obj.gra, mesh, fphys );
        end
    end
    
    methods
        function matUpdateWetDryState(obj, fphys)
            for m = 1:obj.Nmesh
                wetflag = all( fphys{m}(:,:,1) > obj.hmin );
                obj.meshUnion(m).status( ~wetflag ) = int8( enumSWERegion.Dry );
                obj.meshUnion(m).status(  wetflag ) = int8( enumSWERegion.Wet );
            end
        end
    end
    
    methods( Access = protected )
        
        function matEvaluateTopographySourceTerm( obj, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                obj.frhs{m}(:,:,1:3) = obj.frhs{m}(:,:,1:3) + mxEvaluateSourceTopography2d...
                    ( obj.gra, mesh.status, fphys{m}, obj.zGrad{m} );
            end
        end
    end
    
end

