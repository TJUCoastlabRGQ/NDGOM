classdef NdgAbstractVisSolver < handle
    
    properties( SetAccess = protected )
        % num of mesh
        Nmesh
        % physical object
        phys
        % auxiliary variabel, the second index mean the partial derivative with recpect to the corresponding direction, 
        % while the third index the variable in the corresponding direction
        pxx, pyy, pzz, pzx, pzy, pxy, pyx
        % variable index, and its corresponding rhs index
        varId, rhsId
        % num of field
        Nfield
    end
    
    properties( SetAccess = public )
        % viscosity
        mx, my, mz
    end
    
    methods
        function obj = NdgAbstractVisSolver( phys, varId, rhsId )
            obj.phys = phys;
            obj.varId = varId;
            obj.rhsId = rhsId;
            obj.Nfield = numel( varId );
            
            obj.Nmesh = phys.Nmesh;
            obj.pxx = cell( obj.Nmesh, 1 );
            obj.pyy = cell( obj.Nmesh, 1 );
            obj.pzz = cell( obj.Nmesh, 1 );

            obj.pzx = cell( obj.Nmesh, 1 );
            obj.pzy = cell( obj.Nmesh, 1 );
            obj.pxy = cell( obj.Nmesh, 1 );
            obj.pyx = cell( obj.Nmesh, 1 );
%             obj.mx = cell( obj.Nmesh, 1 );
%             obj.my = cell( obj.Nmesh, 1 );
%             obj.mz = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                obj.pxx{m} = zeros( phys.meshUnion(m).cell.Np, phys.meshUnion(m).K );
                obj.pyy{m} = zeros( phys.meshUnion(m).cell.Np, phys.meshUnion(m).K );
                obj.pzz{m} = zeros( phys.meshUnion(m).cell.Np, phys.meshUnion(m).K );
                obj.pzx{m} = zeros( phys.meshUnion(m).cell.Np, phys.meshUnion(m).K );
                obj.pzy{m} = zeros( phys.meshUnion(m).cell.Np, phys.meshUnion(m).K );
                obj.pxy{m} = zeros( phys.meshUnion(m).cell.Np, phys.meshUnion(m).K );          
                obj.pyx{m} = zeros( phys.meshUnion(m).cell.Np, phys.meshUnion(m).K );       
            end
        end% func
    end% methods
    
    methods(Abstract)
        matEvaluateRHS( obj, fphys, frhs );
    end% methods
    
end
