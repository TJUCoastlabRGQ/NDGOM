classdef SWEBarotropic3d < SWEAbstract3d
    %SWEBAROTROPIC3D 此处显示有关此类的摘要
    %   此处显示详细说明
        
%     properties ( Constant )
%         %> number of physical field
%         Nfield2d = 7 %[eta HU HV H Z Zx Zy]
%         %> num of 3d physical field
%         Nfield3d = 11 % [ Hu, Hv, w, tau_x, tau_y, H, eta, Z, Zx, Zy, MUy ]
%         %> number of variable field
%         Nvar3d = 2
%         %> index of variable in physical field
%         varFieldIndex3d = [ 1, 2 ]
%         %> number of variable field
%         Nvar2d = 1
%         %> index of variable in physical field
%         varFieldIndex2d = 1        
%     end
    properties
        %> number of physical field
        Nfield2d = 7 %[eta HU HV H Z Zx Zy]
        %> num of 3d physical field
        Nfield = 11 % [ Hu, Hv, w, tau_x, tau_y, H, eta, Z, Zx, Zy, MUy ]
        %> number of variable field
        Nvar = 2
        %> index of variable in physical field
        varFieldIndex = [ 1, 2 ]
        %> number of variable field
        Nvar2d = 1
        %> index of variable in physical field
        varFieldIndex2d = 1        
    end
    
    properties
        %> the 2d field to be put in the output file
        outputFieldOrder2d =  1 
        %> the 3d field to be put in the output file 
        outputFieldOrder = [1 2 3]
    end
    
    methods
        
%         function matSolve( obj )
%             matEvaluateRK45( obj );
%         end
%         
%         matEvaluateRK45( obj );
    end
    
    methods ( Access = protected )
          matEvaluateRK45( obj );
    end
    
end

