classdef SWEConventional1d < SWEAbstract1d
    
    methods( Hidden )
        function [ E ] = matEvaluateFlux( obj, mesh, fphys )
%             [ E ] = mxEvaluateFlux1d( obj.hmin, obj.gra, mesh.status, fphys );
              [ E ] = obj.volumefluxSolver.evaluate(obj.hmin, obj.gra, mesh, fphys);
        end
    end
    
    methods( Access = protected )
        function matUpdateWetDryState(obj, fphys)
            for m = 1:obj.Nmesh
                wetflag = all( fphys{m}(:,:,1) > obj.hmin );
                obj.meshUnion(m).status( ~wetflag ) = int8( enumSWERegion.Dry );
                obj.meshUnion(m).status(  wetflag ) = int8( enumSWERegion.Wet );
            end
        end
        
        function [ fphys ] = matEvaluatePostFunc(obj, fphys)
            obj.matUpdateWetDryState( fphys );
            for m = 1:obj.Nmesh
                hc = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,1) );
                qxc = obj.meshUnion(m).GetMeshAverageValue( fphys{m}(:,:,2) );
                fphys{m}(:,:,1:2) = mxEvaluatePostFunc1d( obj.hmin, fphys{m}, hc, qxc );
            end
            obj.matUpdateWetDryState( fphys );
        end
        
        function matEvaluateTopographySourceTerm( obj, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                obj.frhs{m}(:,:,1:2) = obj.frhs{m}(:,:,1:2) + mxEvaluateSourceTopography1d...
                    ( obj.gra, mesh.EToR, fphys{m}, obj.zGrad{m} );
            end
        end
    end
    
end

