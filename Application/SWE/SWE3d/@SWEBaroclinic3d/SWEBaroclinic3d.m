classdef SWEBaroclinic3d < SWEBarotropic3d
    %SWEBAROCLINIC3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        %> num of 3d physical field
        Nfield = 15
        %> field name of each field
        fieldName3d = {'hu','hv','omega', 'h','nv','z','eta','zx','zy','w', 'hw','hc','rho','hT', 'hS'};
        %> number of variable field
        Nvar = 4
        %> index of variable in physical field
        varFieldIndex = [ 1, 2, 14, 15 ]
    end
    
    properties
        %> the 3d field to be put in the output file
        outputFieldOrder3d = [1 2 3 12 13]
    end
    
    methods
        function obj = SWEBaroclinic3d(  )
            clear mxCalculateDensityField;
            clear mxCalculateBaroclinicTerm;
        end
    end
    
    methods( Access = protected )
        function rho = matCalculateDensityField( obj, fphys )
            rho = mxCalculateDensityField( fphys(:,:,4), fphys(:,:,12), fphys(:,:,13), obj.meshUnion.z,  obj.hcrit, ...
                obj.meshUnion.cell.Np, obj.meshUnion.K);
        end
        
        function matEvaluateBaroclinicTerm( obj, fphys )
            warning('off');
            [ obj.frhs{1}(:,:,1), obj.frhs{1}(:,:,2) ] = mxCalculateBaroclinicTerm( obj.frhs{1}(:,:,1), obj.frhs{1}(:,:,2), ...
                struct(obj.meshUnion.InnerEdge), struct(obj.meshUnion.BoundaryEdge), struct(obj.meshUnion.BottomEdge), ...
                struct(obj.meshUnion.cell), struct(obj.meshUnioin), fphys{1}(:,:,4), fphys{1}(:,:,14), obj.gra );
            warning('on');
        end
        
    end    
    
end

