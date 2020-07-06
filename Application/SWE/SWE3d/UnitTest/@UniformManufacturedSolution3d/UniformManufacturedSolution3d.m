classdef UniformManufacturedSolution3d < ManufacturedSolution3d
    %MANUFACTUREDSOLUTION3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    methods
        function obj = UniformManufacturedSolution3d( N, Nz, M, Mz )
            obj = obj@ManufacturedSolution3d( N, Nz, M, Mz );
                        obj.ExactValue = cell(1);
%             [obj.ExactValue{1}(:,:,1), obj.ExactValue{1}(:,:,2), obj.ExactValue{1}(:,:,3), ~] = ...
%                obj.matGetExactSolution( obj.meshUnion.x, obj.meshUnion.y, obj.meshUnion.z, obj.ftime); 
           
        end
    end
    
    methods ( Access = protected )
        
        function matEvaluateError( obj, fphys, time)
%             [ h, hu, hv ] = obj.matGetExactSolution( obj.mesh3d, time);
%             fext = cell(1);
%             fext{1}(:,:,1) = hu;
%             fext{1}(:,:,2) = hv;
%             fext{1}(:,:,3) = h;
%             Tempfphys = cell(1);
%             Tempfphys{1}(:,:,1) = fphys(:,:,1);
%             Tempfphys{1}(:,:,2) = fphys(:,:,2);
%             Tempfphys{1}(:,:,3) = fphys(:,:,4);
%             Err2 = obj.PostProcess.evaluateNormErr2( Tempfphys, fext );
%             obj.timePoint(obj.Index) = time;
%             obj.ErrNorm2{1}(obj.Index)  = Err2(1);
%             obj.ErrNorm2{2}(obj.Index)  = Err2(2);
%             obj.ErrNorm2{3}(obj.Index)  = Err2(3);
%             obj.Index = obj.Index + 1;
        end        
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh2d = obj.mesh2d(m);
                mesh3d = obj.mesh3d(m);
                x = mesh3d.x; y = mesh3d.y;
                fphys2d{m} = zeros( mesh2d.cell.Np, mesh2d.K, obj.Nfield2d );
                fphys{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield );
                % bottom elevation
                fphys2d{m}(:, :, 4) = -(2-0.005*(mesh2d.x+mesh2d.y));
                
                fphys2d{m}(:,:,5) = ( mesh2d.rx .* ( mesh2d.cell.Dr * fphys2d{m}(:,:,4) ) ) + ...
                    ( mesh2d.sx .* ( mesh2d.cell.Ds * fphys2d{m}(:,:,4) ) );
                
                fphys2d{m}(:,:,6) = ( mesh2d.ry .* ( mesh2d.cell.Dr * fphys2d{m}(:,:,4) ) ) + ...
                    ( mesh2d.sy .* ( mesh2d.cell.Ds * fphys2d{m}(:,:,4) ) );
                %water depth
                fphys2d{m}(:,:,1) = obj.e*(sin(obj.w*(mesh2d.x+0)) + sin(obj.w*(mesh2d.y+0))) - ...
                    fphys2d{m}(:, :, 4);
                
                %> Z in extended three dimensional fields
                fphys{m}(:,:,6) = mesh3d.Extend2dField( fphys2d{m}(:, :, 4) );
                %> '$\eta$' in extended three dimensional fields
                fphys{m}(:,:,7) = fphys{m}(:,:,4) + fphys{m}(:,:,6);
                
                [ fphys{m}(:,:,4), fphys{m}(:,:,1), fphys{m}(:,:,2), ~ ] = obj.matGetExactSolution( mesh3d.x, mesh3d.y, mesh3d.z, 0);
                
                [ ~, fphys2d{m}(:,:,2), fphys2d{m}(:,:,3), ~ ] = obj.matGetExactSolution( mesh2d.x, mesh2d.y, mesh2d.x, 0);
%                 fphys2d{m}(:, :, 2) = mesh3d.VerticalColumnIntegralField( fphys{m}(:, :, 1) );
%                 fphys2d{m}(:, :, 3) = mesh3d.VerticalColumnIntegralField( fphys{m}(:, :, 2) );
                
            end
        end
        
        function  matEvaluateSourceTerm(obj, fphys, time)
            mesh = obj.meshUnion(1);
            mesh2d = obj.mesh2d(1);
            matEvaluateSourceTerm@SWEAbstract3d( obj, fphys);
            
            h2 = @(x,y,t)obj.e * ( sin(obj.w*(x+t)) + sin(obj.w*(y+t)) ) + 2 - 0.005*( x + y );
            h2d = h2(mesh2d.x, mesh2d.y, time);
            
%             h2d =  obj.e * ( sin(obj.w*(mesh2d.x+time)) + sin(obj.w*(mesh2d.y+time)) ) + 2 - 0.005*( mesh2d.x + mesh2d.y );
            h2dt = obj.e * obj.w .* cos(obj.w*(mesh2d.x+time)) + obj.e * obj.w .* cos(obj.w*(mesh2d.y+time));
            u2d = sin(obj.w*(mesh2d.x+time));
            u2dx = obj.w * cos(obj.w*(mesh2d.x + time));
            v2dy = obj.w * cos(obj.w*(mesh2d.y + time));
            v2d = sin(obj.w*(mesh2d.y+time));
            h2dx = obj.e * obj.w *( cos(obj.w*(mesh2d.x+time)) ) - 0.005;
            h2dy = obj.e * obj.w *( cos(obj.w*(mesh2d.y+time)) ) - 0.005;
            obj.frhs2d{1}(:,:,1) = obj.frhs2d{1}(:,:,1) +  h2dt  + ...
                h2d .* u2dx + u2d.*h2dx + h2d .* v2dy + v2d .* h2dy;
            
            h3d = obj.e*(sin(obj.w.*(mesh.x+time)) + sin(obj.w*(mesh.y+time))) + (2-0.005.*(mesh.x + mesh.y));
            eta = obj.e*(sin(obj.w.*(mesh.x+time)) + sin(obj.w*(mesh.y+time)));
            h3dt = obj.e * obj.w * ( cos(obj.w.*(mesh.x+time)) + cos(obj.w.*(mesh.y+time)) );
            u3d = sin(obj.w.*(mesh.x+time));
            u3dt = obj.w .* cos(obj.w.*(mesh.x+time));
            v3d = sin(obj.w.*(mesh.y+time));
            v3dt = obj.w .* cos(obj.w.*(mesh.y+time));
            
            h3dx = obj.e * obj.w .* cos(obj.w.*(mesh.x + time)) - 0.005;
            h3dy = obj.e * obj.w .* cos(obj.w.*(mesh.y + time)) - 0.005;
            u3dx = obj.w .* cos(obj.w.*(mesh.x + time));
            v3dy = obj.w .* cos(obj.w.*(mesh.y + time));
            
            
            obj.frhs{1}(:,:,1) = obj.frhs{1}(:,:,1) + u3d.*h3dt + h3d.*u3dt + u3d.*(u3d.*h3dx+h3d.*u3dx) + h3d.*u3d.*u3dx + obj.gra.*h3d.*h3dx - ...
                obj.gra .* fphys{1}(:,:,6).*fphys{1}(:,:,8) + u3d .* (v3d.*h3dy + h3d.*v3dy) + obj.gra .* eta .* fphys{1}(:,:,8);
            obj.frhs{1}(:,:,2) = obj.frhs{1}(:,:,2) + v3d.*h3dt + h3d.*v3dt + v3d.*(u3d.*h3dx+h3d.*u3dx) + v3d.*(v3d.*h3dy + h3d.*v3dy) + h3d.*v3d.*v3dy + ...
                obj.gra.*h3d.*h3dy - obj.gra .* fphys{1}(:,:,6).*fphys{1}(:,:,9) + obj.gra .* eta .* fphys{1}(:,:,9);
        end
        
        function matUpdateExternalField( obj, time, ~, ~ )

             [obj.fext3d{1}( :, :, 4 ), obj.fext3d{1}( :, :, 1 ), obj.fext3d{1}( :, :, 2 ), ~] = ...
                obj.matGetExactSolution( obj.mesh3d.BoundaryEdge.xb, obj.mesh3d.BoundaryEdge.yb, obj.mesh3d.BoundaryEdge.zb, time);
            
           [ obj.fext2d{1}( :, :, 1 ), obj.fext2d{1}( :, :, 2 ), obj.fext2d{1}( :, :, 3 ), ~ ] =...
               obj.matGetExactSolution( obj.mesh2d.BoundaryEdge.xb, obj.mesh2d.BoundaryEdge.yb, obj.mesh2d.BoundaryEdge.xb, time);
            
%             obj.fext2d{1}( :, :, 2 ) = obj.meshUnion.BoundaryEdge.VerticalColumnIntegralField(  obj.fext3d{1}( :, :, 1 ) );      
%             obj.fext2d{1}( :, :, 3 ) = obj.meshUnion.BoundaryEdge.VerticalColumnIntegralField(  obj.fext3d{1}( :, :, 2 ) );            
            
        end
        
        function [ h, hu, hv, Omega ] = matGetExactSolution( obj, x, y, ~, time)
            h = obj.e*(sin(obj.w.*(x+time)) + sin(obj.w*(y+time))) + (2-0.005.*(x + y));
            hu = sin(obj.w.*(x+time)) .* h;
            hv = sin(obj.w.*(y+time)) .* h;
            Omega = [];
        end
        
        
        function [ option ] = setOption( obj, option )
            ftime = 100;
            outputIntervalNum = 100;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 5;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK343;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.None;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('ConstantVerticalEddyViscosityValue') = 0.03;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 100;
        end
        
        
        
    end
    
end