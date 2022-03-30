classdef SWEBaroclinic3d < SWEBarotropic3d
    %SWEBAROCLINIC3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        
    end
    
    properties
        %> The following parameter is used for calculation of density when the linear
        %> EOS is adopted, the linear EOS is taken as $\rho = \rho_0 - \alpha_T(T-T_0) + \beta_S(S-S_0) $
        T0 = 25
        
        S0 = 35
         
        alphaT = 0.2
        
        betaS = 0.1
    end
    
    methods
        function obj = SWEBaroclinic3d(  )
            %> num of 3d physical field
            obj.Nfield = 15;
            %> field name of each field
            obj.fieldName3d = {'hu','hv','omega', 'h','nv','z','eta','zx','zy','w', 'hw','hc','rho','hT', 'hS'};
            %> number of variable field
            obj.Nvar = 4;
            %> index of variable in physical field
            obj.varFieldIndex = [ 1, 2, 14, 15 ];
            %> the 3d field to be put in the output file
            obj.outputFieldOrder3d = [1 2 3 13 14 15];
            
            clear mxCalculateDensityField;
            
            clear mxCalculateBaroclinicTerm;
        end
    end
    
    methods( Access = protected )
        function rho = matCalculateDensityField( obj, fphys )
            
            rho = mxCalculateDensityField( fphys(:,:,4), fphys(:,:,14), fphys(:,:,15), obj.meshUnion.z,  obj.hcrit, ...
                obj.meshUnion.cell.Np, obj.meshUnion.K, obj.rho0, obj.alphaT, obj.betaS, obj.T0, obj.S0, char(obj.getOption('EosType')));
            
        end
        
        function matEvaluateBaroclinicTerm( obj, fphys )
            warning('off');
            [ obj.frhs{1}(:,:,1), obj.frhs{1}(:,:,2) ] = mxCalculateBaroclinicTerm( obj.frhs{1}(:,:,1), obj.frhs{1}(:,:,2), ...
                struct(obj.meshUnion.InnerEdge), struct(obj.meshUnion.BoundaryEdge), struct(obj.meshUnion.BottomEdge), ...
                struct(obj.meshUnion.cell), struct(obj.meshUnion), fphys{1}(:,:,4), fphys{1}(:,:,13), obj.gra );
            warning('on');
        end
        
    end
    
end

