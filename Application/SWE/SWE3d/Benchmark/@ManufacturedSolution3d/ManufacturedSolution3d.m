classdef ManufacturedSolution3d < SWEBarotropic3d
    %MANUFACTUREDSOLUTION3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        PostProcess
        timePoint
        ErrNorm2
        Index
    end
    
    properties ( Constant )
        %> channel length
        ChLength = 100;
        ChWidth = 100;
        e = 0.001;
        w = 0.01;
        d = 0.1;
        hcrit = 0.001;
    end
    
    methods
        function obj = ManufacturedSolution3d( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder3d = [ 1 2 3 10];
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            %> time interval
            obj.PostProcess = NdgPostProcess(obj.meshUnion(1),...
                strcat('WindDrivenFlow/3d','/','WindDrivenFlow'));
            obj.ErrNorm2 = cell(3);
            obj.Index = 1;
        end
        
        function matPlotErrorNorm2(obj)
            close all
            figure(1);
            hold on;
            color = {'r-','k-','b-','r-'};
            LWidth = 1.5;
            FontSize = 15;
            for i = 1:4
                plot(obj.timePoint, obj.ErrNorm2{i},color{i},'LineWidth',LWidth);
            end
            lendstr = {'$hu$',...
                '$hv$','$h$','$\omega$'};
            legend(lendstr,'Interpreter','Latex','FontSize',FontSize);
            xlabel('$time(s)$','Interpreter','Latex');
            ylabel('$L_2 error$','Interpreter','Latex');
            box on;
        end
    end
    
    methods ( Access = protected )
        
        matCalculateExplicitRHSTerm( obj, fphys2d, fphys, Stage, RKIndex, time);
        
        function matEvaluateError( obj, fphys, time)
            [ h, hu, hv, Omega ] = obj.matGetExactSolution( obj.mesh3d, time);
            fext{1}(:,:,1) = hu;
            fext{1}(:,:,2) = hv;
            fext{1}(:,:,3) = h;
            fext{1}(:,:,4) = Omega;
            Tempfphys = cell(1);
            Tempfphys{1}(:,:,1) = fphys(:,:,1);
            Tempfphys{1}(:,:,2) = fphys(:,:,2);
            Tempfphys{1}(:,:,3) = fphys(:,:,4);
            Tempfphys{1}(:,:,4) = fphys(:,:,3);
            Err2 = obj.PostProcess.evaluateNormErr2( Tempfphys, fext );
            obj.timePoint(obj.Index) = time;
            obj.ErrNorm2{1}(obj.Index)  = Err2(1);
            obj.ErrNorm2{2}(obj.Index)  = Err2(2);
            obj.ErrNorm2{3}(obj.Index)  = Err2(3);
            obj.ErrNorm2{4}(obj.Index)  = Err2(4);
            obj.Index = obj.Index + 1;
        end
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh2d = obj.mesh2d(m);
                mesh3d = obj.mesh3d(m);
                z = mesh3d.z; x = mesh3d.x; y = mesh3d.y;
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
                
                fphys{m}(:,:,4) = mesh3d.Extend2dField(fphys2d{m}(:,:,1));
                %> Z in extended three dimensional fields
                fphys{m}(:,:,6) = mesh3d.Extend2dField( fphys2d{m}(:, :, 4) );
                %> '$\eta$' in extended three dimensional fields
                fphys{m}(:,:,7) = fphys{m}(:,:,4) + fphys{m}(:,:,6);
                
                fphys{m}(:,:,1) = -fphys{m}(:,:,4).*obj.d.*z.*sin(obj.w*(x+0)) .* fphys{m}(:,:,4);
                fphys{m}(:,:,2) = -fphys{m}(:,:,4).*obj.d.*z.*sin(obj.w*(y+0)) .* fphys{m}(:,:,4);
                
                fphys2d{m}(:, :, 2) = mesh3d.VerticalColumnIntegralField( fphys{m}(:, :, 1) );
                fphys2d{m}(:, :, 3) = mesh3d.VerticalColumnIntegralField( fphys{m}(:, :, 2) );
                
                %Omega
                fphys{m}(:,:,3) = -fphys{m}(:,:,4).^2 .* obj.d./2.*obj.w.*z.*cos(obj.w.*(x+0)) - ...
                    fphys{m}(:,:,4) * obj.d ./2 .* sin(obj.w.*(x+0)).*z.*(obj.e*obj.w.*cos(obj.w.*(x+0)) - 0.005) + ...
                    fphys{m}(:,:,4).^2 .* obj.d .* z.^2./2*obj.w.*cos(obj.w.*(x+0)) + ...
                    fphys{m}(:,:,4) .* obj.d .* z.^2./2 .* sin(obj.w.*(x+0)).*(obj.e*obj.w.*cos(obj.w.*(x+0)) - 0.005) - ...
                    fphys{m}(:,:,4).^2 .* obj.d./2.*obj.w.*z.*cos(obj.w.*(y+0)) - ...
                    fphys{m}(:,:,4) * obj.d ./2 .* sin(obj.w.*(y+0)).*z.*(obj.e*obj.w.*cos(obj.w.*(y+0)) - 0.005) + ...
                    fphys{m}(:,:,4).^2 .* obj.d .* z.^2./2*obj.w.*cos(obj.w.*(y+0)) + ...
                    fphys{m}(:,:,4) .* obj.d .* z.^2./2 .* sin(obj.w.*(y+0)).*(obj.e*obj.w.*cos(obj.w.*(y+0)) - 0.005);
                
            end
        end
        
        function  matEvaluateSourceTerm(obj, fphys, time)
            mesh = obj.meshUnion(1);
            mesh2d = obj.mesh2d(1);
            matEvaluateSourceTerm@SWEAbstract3d( obj, fphys);
            h2d =  obj.e * ( sin(obj.w*(mesh2d.x+time)) + sin(obj.w*(mesh2d.y+time)) ) + 2 - 0.005*( mesh2d.x + mesh2d.y );
            %[Problematic
            obj.frhs2d{1}(:,:,1) = obj.frhs2d{1}(:,:,1) + obj.e*(obj.w * cos(obj.w.*(mesh2d.x+time)) + obj.w * cos(obj.w.*(mesh2d.y+time))) - ...
                h2d .* h2d * obj.d *obj.w./2.*cos(obj.w*(mesh2d.x+time)) - h2d.*obj.d./2.*sin(obj.w*(mesh2d.x+time)).*( obj.e * obj.w * cos(obj.w.*(mesh2d.x+time)) - 0.005 ) - ...
                h2d .* h2d * obj.d *obj.w./2.*cos(obj.w*(mesh2d.y+time)) - h2d.*obj.d./2.*sin(obj.w*(mesh2d.y+time)).*( obj.e * obj.w * cos(obj.w.*(mesh2d.y+time)) - 0.005 );
            
            h = obj.e*(sin(obj.w.*(mesh.x+time)) + sin(obj.w*(mesh.y+time))) + (2-0.005.*(mesh.x + mesh.y));
            eta = obj.e*(sin(obj.w.*(mesh.x+time)) + sin(obj.w*(mesh.y+time)));
            ht = obj.e*(obj.w*cos(obj.w.*(mesh.x+time)) + obj.w*cos(obj.w.*(mesh.y+time)));
            u = -h.*obj.d.*mesh.z.*sin(obj.w.*(mesh.x+time));
            ut = -h .* obj.d .* mesh.z .* obj.w .* cos(obj.w.*(mesh.x+time));
            v = -h.*obj.d.*mesh.z.*sin(obj.w.*(mesh.y+time));
            vt = -h .* obj.d .* mesh.z .* obj.w .* cos(obj.w.*(mesh.y+time));
            
            hx = obj.e * obj.w .* cos(obj.w.*(mesh.x + time)) - 0.005;
            hy = obj.e * obj.w .* cos(obj.w.*(mesh.y + time)) - 0.005;
            ux = -h .* obj.d .* mesh.z .*  obj.w .* cos(obj.w.*(mesh.x + time));
            vy = -h .* obj.d .* mesh.z .*  obj.w .* cos(obj.w.*(mesh.y + time));
            
            us = -h .* obj.d .* sin(obj.w.*(mesh.x + time));
            vs = -h .* obj.d .* sin(obj.w.*(mesh.y + time));
            
            Omega = -h.^2 .* obj.d./2.*obj.w.*mesh.z.*cos(obj.w.*(mesh.x+time)) - ...
                h * obj.d ./2 .* sin(obj.w.*(mesh.x+time)).*mesh.z.*(obj.e*obj.w.*cos(obj.w.*(mesh.x+time)) - 0.005) + ...
                h.^2 .* obj.d .* mesh.z.^2./2*obj.w.*cos(obj.w.*(mesh.x+time)) + ...
                h .* obj.d .* mesh.z.^2./2 .* sin(obj.w.*(mesh.x+time)).*(obj.e*obj.w.*cos(obj.w.*(mesh.x+time)) - 0.005) - ...
                h.^2 .* obj.d./2.*obj.w.*mesh.z.*cos(obj.w.*(mesh.y+time)) - ...
                h * obj.d ./2 .* sin(obj.w.*(mesh.y+time)).*mesh.z.*(obj.e*obj.w.*cos(obj.w.*(mesh.y+time)) - 0.005) + ...
                h.^2 .* obj.d .* mesh.z.^2./2*obj.w.*cos(obj.w.*(mesh.y+time)) + ...
                h .* obj.d .* mesh.z.^2./2 .* sin(obj.w.*(mesh.y+time)).*(obj.e*obj.w.*cos(obj.w.*(mesh.y+time)) - 0.005);
            
            Os = - h.^2 .* obj.d./2.*obj.w.*cos(obj.w.*(mesh.x+time)) - ...
                h * obj.d ./2 .* sin(obj.w.*(mesh.x+time)).*(obj.e*obj.w.*cos(obj.w.*(mesh.x+time)) - 0.005) + ...
                h.^2 .* obj.d .* mesh.z*obj.w.*cos(obj.w.*(mesh.x+time)) + ...
                h .* obj.d .* mesh.z .* sin(obj.w.*(mesh.x+time)).*(obj.e*obj.w.*cos(obj.w.*(mesh.x+time)) - 0.005) - ...
                h.^2 .* obj.d./2.*obj.w.*cos(obj.w.*(mesh.y+time)) - ...
                h * obj.d ./2 .* sin(obj.w.*(mesh.y+time)).*(obj.e*obj.w.*cos(obj.w.*(mesh.y+time)) - 0.005) + ...
                h.^2 .* obj.d .* mesh.z .* obj.w.*cos(obj.w.*(mesh.y+time)) + ...
                h .* obj.d .* mesh.z .* sin(obj.w.*(mesh.y+time)).*(obj.e*obj.w.*cos(obj.w.*(mesh.y+time)) - 0.005);
            
            obj.frhs{1}(:,:,1) = obj.frhs{1}(:,:,1) + u.*ht + h.*ut + u.*(u.*hx+h.*ux) + h.*u.*ux + obj.gra.*h.*hx - ...
                obj.gra .* fphys{1}(:,:,6).*fphys{1}(:,:,8) + u .* (v.*hy + h.*vy) + u.*Os + Omega .* us + obj.gra .* eta .* fphys{1}(:,:,8);
            obj.frhs{1}(:,:,2) = obj.frhs{1}(:,:,2) + v.*ht + h.*vt + v.*(u.*hx+h.*ux) + v.*(v.*hy + h.*vy) + h.*v.*vy + ...
                obj.gra.*h.*hy - obj.gra .* fphys{1}(:,:,6).*fphys{1}(:,:,9) + v.*Os + Omega .* vs + obj.gra .* eta .* fphys{1}(:,:,9);
        end
        
        function matUpdateExternalField( obj, time, ~, ~ ) 
            h = obj.e * ( sin(obj.w*(obj.meshUnion.BoundaryEdge.xb+time)) + sin(obj.w*(obj.meshUnion.BoundaryEdge.yb+time)) )...
                + 2 - 0.005*( obj.meshUnion.BoundaryEdge.xb + obj.meshUnion.BoundaryEdge.yb );
            obj.fext3d{1}( :, :, 4 ) = h;
            obj.fext3d{1}( :, :, 1 ) = obj.fext3d{1}( :, :, 4 ) .* (-h.*obj.d.*obj.meshUnion.BoundaryEdge.zb.*sin(obj.w.*(obj.meshUnion.BoundaryEdge.xb+time)));
            obj.fext3d{1}( :, :, 2 ) = obj.fext3d{1}( :, :, 4 ) .* (-h.*obj.d.*obj.meshUnion.BoundaryEdge.zb.*sin(obj.w.*(obj.meshUnion.BoundaryEdge.yb+time)));
            
            obj.fext2d{1}( :, :, 1 ) = obj.e * ( sin(obj.w*(obj.mesh2d.BoundaryEdge.xb+time)) + sin(obj.w*(obj.mesh2d.BoundaryEdge.yb+time)) )...
                + 2 - 0.005*( obj.mesh2d.BoundaryEdge.xb + obj.mesh2d.BoundaryEdge.yb );
            
            obj.fext2d{1}( :, :, 2 ) = obj.meshUnion.BoundaryEdge.VerticalColumnIntegralField(  obj.fext3d{1}( :, :, 1 ) );      
            obj.fext2d{1}( :, :, 3 ) = obj.meshUnion.BoundaryEdge.VerticalColumnIntegralField(  obj.fext3d{1}( :, :, 2 ) );            
            
        end        
        
        
        function [ h, hu, hv, Omega ] = matGetExactSolution( obj, mesh, time)
            h = obj.e*(sin(obj.w.*(mesh.x+time)) + sin(obj.w*(mesh.y+time))) + (2-0.005.*(mesh.x + mesh.y));
            
            hu = -h.*obj.d.*mesh.z.*sin(obj.w.*(mesh.x+time)) .* h;
            hv = -h.*obj.d.*mesh.z.*sin(obj.w.*(mesh.y+time)) .* h;
            
            Omega = -h.^2 .* obj.d./2.*obj.w.*mesh.z.*cos(obj.w.*(mesh.x+time)) - ...
                h * obj.d ./2 .* sin(obj.w.*(mesh.x+time)).*mesh.z.*(obj.e*obj.w.*cos(obj.w.*(mesh.x+time)) - 0.005) + ...
                h.^2 .* obj.d .* mesh.z.^2./2*obj.w.*cos(obj.w.*(mesh.x+time)) + ...
                h .* obj.d .* mesh.z.^2./2 .* sin(obj.w.*(mesh.x+time)).*(obj.e*obj.w.*cos(obj.w.*(mesh.x+time)) - 0.005) - ...
                h.^2 .* obj.d./2.*obj.w.*mesh.z.*cos(obj.w.*(mesh.y+time)) - ...
                h * obj.d ./2 .* sin(obj.w.*(mesh.y+time)).*mesh.z.*(obj.e*obj.w.*cos(obj.w.*(mesh.y+time)) - 0.005) + ...
                h.^2 .* obj.d .* mesh.z.^2./2*obj.w.*cos(obj.w.*(mesh.y+time)) + ...
                h .* obj.d .* mesh.z.^2./2 .* sin(obj.w.*(mesh.y+time)).*(obj.e*obj.w.*cos(obj.w.*(mesh.y+time)) - 0.005);            
           
        end
        
        
        function [ option ] = setOption( obj, option )
            ftime = 200;
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

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

% bctype = [ ...
%     enumBoundaryCondition.SlipWall, ...
%     enumBoundaryCondition.SlipWall, ...
%     enumBoundaryCondition.ZeroGrad, ...
%     enumBoundaryCondition.ZeroGrad ];

bctype = [ ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped ];

mesh2d = makeUniformQuadMesh( N, ...
    [ 0, obj.ChLength ], [ 0, obj.ChWidth ], ceil(obj.ChLength/M), ceil(obj.ChWidth/M), bctype);

cell = StdPrismQuad( N, Nz );
zs = ones(mesh2d.Nv, 1); zb = zeros(mesh2d.Nv, 1);
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
end