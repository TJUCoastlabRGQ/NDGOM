classdef NdgFilter3d < NdgAbstractFilter
    %NDGFILTER3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgFilter3d( mesh )
            obj = obj@NdgAbstractFilter( mesh );
            obj.matAssembleFilterVandMatirx( mesh );
            obj.matSetSmoothRequireParameter( mesh );
            obj.matSetFilterStrengthMatrix( mesh );
        end
    end
    
    methods( Hidden )
        function matAssembleFilterVandMatirx( obj, mesh )
            obj.VF = mesh.cell.V;
            if mesh.cell.type == enumStdCell.PrismQuad
                % x direction
                obj.VF(:,mod( 1:mesh.cell.Np, mesh.cell.N+1 ) == 0 ) = 0;
                % y direction
                obj.VF(:,repmat(fix((0:(mesh.cell.N+1)*(mesh.cell.N+1)-1)/(mesh.cell.N + 1)),1,(mesh.cell.Nz+1)) == 1 ) = 0;
                %> z direction
                obj.VF(:,  fix(((1:mesh.cell.Np)-1)./(mesh.cell.N+1)/(mesh.cell.N+1)) == 1 ) = 0;
            end
        end
        
        function  matSetSmoothRequireParameter( obj, mesh )
            obj.SR = -2 - 6*log10(mesh.cell.N);
        end
        
        function matSetFilterStrengthMatrix( obj, mesh )
           CoeX = rem((1:mesh.cell.Np)-1,mesh.cell.N+1)/mesh.cell.N;
           CoeY = repmat(fix((0:(mesh.cell.N+1)*(mesh.cell.N+1)-1)/(mesh.cell.N + 1))./mesh.cell.N,1,(mesh.cell.Nz+1));
           CoeZ = fix(((1:mesh.cell.Np)-1)./(mesh.cell.N+1)/(mesh.cell.N+1))/mesh.cell.Nz;
           obj.FS = exp( -obj.alpha * CoeX.^obj.r ) .* exp( -obj.alpha * CoeY.^obj.r ) .* ...
               exp( -obj.alpha * CoeZ.^obj.r );
        end
    end
    
end

