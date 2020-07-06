classdef ManufacturedSolution3d < SWEBarotropic3d
    %MANUFACTUREDSOLUTION3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        PostProcess
        timePoint
        ErrNorm2
        Index
        ExactValue
    end
    
    properties
        h
        ht
        hx
        hy
        eta
        u
        ux
        uy
        uz
        ut
        v
        vx
        vy
        vz
        vt
        u2d
        u2dx
        u2dy
        u2dt
        v2d
        v2dx
        v2dy
        v2dt
        Ou2d
        Ov2d
        OH
        Ou
        Ov
        Oux
        Ovy
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
            obj.matGetFunction();
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            
            
            %> time interval
            obj.PostProcess = NdgPostProcess(obj.meshUnion(1),...
                strcat('WindDrivenFlow/3d','/','WindDrivenFlow'));
            obj.ErrNorm2 = cell(3);
            obj.Index = 1;
                        obj.ExactValue = cell(1);
            [obj.ExactValue{1}(:,:,1), obj.ExactValue{1}(:,:,2), obj.ExactValue{1}(:,:,3), obj.ExactValue{1}(:,:,4)] = ...
               obj.matGetExactSolution( obj.meshUnion.x, obj.meshUnion.y, obj.meshUnion.z, obj.ftime);             
        end
        
        function matPlotErrorNorm2(obj)
            close all
            figure(1);
            hold on;
            color = {'r-','k-','b-','r-'};
            LWidth = 1.5;
            FontSize = 15;
            for i = 1:3
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
        
        [ Omega , W ] = matEvaluateVerticalVelocity( obj, mesh3d, fphys2d, fphys3d, time );
        matCalculateExplicitRHSTerm( obj, fphys2d, fphys, Stage, RKIndex, time);
        
        function matGetFunction(obj)
            syms x y z t;
           obj.h = @(x,y,t)obj.e * ( sin(obj.w*(x+t)) + sin(obj.w*(y+t)) ) + 2 - 0.005*( x + y );
           obj.ht = @(x,y,t)cos(t/100 + x/100)/100000 + cos(t/100 + y/100)/100000;
           obj.hx = @(x,y,t)cos(t/100 + x/100)/100000 - 1/200;
           obj.hy = @(x,y,t)cos(t/100 + y/100)/100000 - 1/200;
           obj.eta = @(x,y,t)obj.e * ( sin(obj.w*(x+t)) + sin(obj.w*(y+t)) );
           obj.u = @(x,y,z,t)-(obj.e * ( sin(obj.w*(x+t)) + sin(obj.w*(y+t)) ) + 2 - 0.005*( x + y ))*obj.d.*( z + 1 ).*sin(obj.w.*(x+t));
           obj.uz = @(x,y,t)-sin(t/100 + x/100).*(sin(t/100 + x/100)/10000 - y/2000 - x/2000 + sin(t/100 + y/100)/10000 + 1/5);
           obj.ut = @(x,y,z,t)- (cos(t/100 + x/100).*(z + 1).*(sin(t/100 + x/100)/10000 - y/2000 - x/2000 + sin(t/100 + y/100)/10000 + 1/5))/100 ...
           - sin(t/100 + x/100).*(z + 1).*(cos(t/100 + x/100)/1000000 + cos(t/100 + y/100)/1000000);
           obj.ux = @(x,y,z,t)- (cos(t/100 + x/100).*(z + 1).*(sin(t/100 + x/100)/10000 - y/2000 - x/2000 + ...
               sin(t/100 + y/100)/10000 + 1/5))/100 - sin(t/100 + x/100).*(cos(t/100 + x/100)/1000000 - 1/2000).*(z + 1);
           obj.uy = @(x,y,z,t)-sin(t/100 + x/100).*(cos(t/100 + y/100)/1000000 - 1/2000).*(z + 1);
           obj.v = @(x,y,z,t)-(obj.e * ( sin(obj.w*(x+t)) + sin(obj.w*(y+t)) ) + 2 - 0.005*( x + y ))*obj.d.*( z + 1 ).*sin(obj.w.*(y+t)); 
           obj.vz = @(x,y,t)-sin(t/100 + y/100) .* (sin(t/100 + x/100)/10000 - y/2000 - x/2000 + sin(t/100 + y/100)/10000 + 1/5);
           obj.vt = @(x,y,z,t)- (cos(t/100 + y/100).*(z + 1).*(sin(t/100 + x/100)/10000 - y/2000 - x/2000 + sin(t/100 + y/100)/10000 + 1/5))/100 ...
           - sin(t/100 + y/100).*(z + 1).*(cos(t/100 + x/100)/1000000 + cos(t/100 + y/100)/1000000);
           obj.vx = @(x,y,z,t)-sin(t/100 + y/100).*(cos(t/100 + x/100)/1000000 - 1/2000).*(z + 1);
           obj.vy = @(x,y,z,t)- (cos(t/100 + y/100).*(z + 1).*(sin(t/100 + x/100)/10000 - y/2000 - x/2000 + ...
               sin(t/100 + y/100)/10000 + 1/5))/100 - sin(t/100 + y/100).*(cos(t/100 + y/100)/1000000 - 1/2000).*(z + 1);
           obj.u2d = @(x,y,t)-(sin(t/100 + x/100).*(sin(t/100 + x/100)/10000 - y/2000 - x/2000 + sin(t/100 + y/100)/10000 + 1/5))/2;
           obj.u2dt = @(x,y,t)- (cos(t/100 + x/100).*(sin(t/100 + x/100)/10000 - y/2000 - x/2000 + sin(t/100 + y/100)/10000 + 1/5))/200 ...
           - (sin(t/100 + x/100).*(cos(t/100 + x/100)/1000000 + cos(t/100 + y/100)/1000000))/2;
           obj.u2dx = @(x,y,t)- (sin(t/100 + x/100).*(cos(t/100 + x/100)/1000000 - 1/2000))/2 - (cos(t/100 + x/100).*(sin(t/100 + x/100)/10000 -...
               y/2000 - x/2000 + sin(t/100 + y/100)/10000 + 1/5))/200;
           obj.u2dy = @(x,y,t)-(sin(t/100 + x/100).*(cos(t/100 + y/100)/1000000 - 1/2000))/2;
           obj.v2d = @(x,y,t)-(sin(t/100 + y/100).*(sin(t/100 + x/100)/10000 - y/2000 - x/2000 + sin(t/100 + y/100)/10000 + 1/5))/2;
           obj.v2dt = @(x,y,t)- (cos(t/100 + y/100).*(sin(t/100 + x/100)/10000 - y/2000 - x/2000 + sin(t/100 + y/100)/10000 + 1/5))/200 ...
               - (sin(t/100 + y/100).*(cos(t/100 + x/100)/1000000 + cos(t/100 + y/100)/1000000))/2;
           obj.v2dy = @(x,y,t)- (sin(t/100 + y/100).*(cos(t/100 + y/100)/1000000 - 1/2000))/2 - (cos(t/100 + y/100).*(sin(t/100 + x/100)/10000 - ...
               y/2000 - x/2000 + sin(t/100 + y/100)/10000 + 1/5))/200;
           obj.v2dx = @(x,y,t)-(sin(t/100 + y/100).*(cos(t/100 + x/100)/1000000 - 1/2000))/2;
           obj.Ou2d = @(x,y,z,t)-(sin(t/100 + x/100).*(z + 1).*(sin(t/100 + x/100)/10000 - y/2000 - x/2000 + sin(t/100 + y/100)/10000 + 1/5))/2;
           obj.OH = @(x,y,z,t)((z + 1).*(sin(t/100 + x/100) - 5*y - 5*x + sin(t/100 + y/100) + 2000))/1000;
           obj.Ou = @(x,y,z,t)-(sin(t/100 + x/100).*(z + 1).^2.*(sin(t/100 + x/100) - 5*y - 5*x + sin(t/100 + y/100) + 2000))/20000;
           obj.Oux = @(x,y,z,t)-((z + 1).^2.*(2000*cos(t/100 + x/100) - 500*sin(t/100 + x/100) - 5*x.*cos(t/100 + x/100) - ...
               5*y.*cos(t/100 + x/100) + 2*cos(t/100 + x/100).*sin(t/100 + x/100) + cos(t/100 + x/100).*sin(t/100 + y/100)))/2000000;
           obj.Ov2d = @(x,y,z,t)-(sin(t/100 + y/100).*(z + 1).*(sin(t/100 + x/100)/10000 - y/2000 - x/2000 + sin(t/100 + y/100)/10000 + 1/5))/2;
           obj.Ov = @(x,y,z,t)-(sin(t/100 + y/100).*(z + 1).^2.*(sin(t/100 + x/100) - 5*y - 5*x + sin(t/100 + y/100) + 2000))/20000;
           obj.Ovy = @(x,y,z,t)-((z + 1).^2.*(2000*cos(t/100 + y/100) - 500*sin(t/100 + y/100) - 5.*x.*cos(t/100 + y/100) - ...
               5.*y.*cos(t/100 + y/100) + cos(t/100 + y/100).*sin(t/100 + x/100) + 2*cos(t/100 + y/100).*sin(t/100 + y/100)))/2000000;
        end
        
        function matEvaluateError( obj, fphys, time)
            fext{1}(:,:,1) = obj.h(obj.meshUnion.x, obj.meshUnion.y, time) ...
                .* obj.u(obj.meshUnion.x, obj.meshUnion.y, obj.meshUnion.z, time);
            fext{1}(:,:,2) = obj.h(obj.meshUnion.x, obj.meshUnion.y, time) ...
                .* obj.v(obj.meshUnion.x, obj.meshUnion.y, obj.meshUnion.z, time);
            fext{1}(:,:,3) = obj.h(obj.meshUnion.x, obj.meshUnion.y, time);
            fext{1}(:,:,4) = obj.Ou2d(obj.meshUnion.x, obj.meshUnion.y, obj.meshUnion.z, time) .* obj.hx(obj.meshUnion.x, obj.meshUnion.y, time) + ...
                obj.OH(obj.meshUnion.x, obj.meshUnion.y, obj.meshUnion.z, time) .* obj.u2dx(obj.meshUnion.x, obj.meshUnion.y, time) - ...
                obj.Ou(obj.meshUnion.x, obj.meshUnion.y, obj.meshUnion.z, time) .* obj.hx(obj.meshUnion.x, obj.meshUnion.y, time) - ...
                obj.h(obj.meshUnion.x, obj.meshUnion.y, time) .* obj.Oux(obj.meshUnion.x, obj.meshUnion.y, obj.meshUnion.z, time) + ...
            obj.Ov2d(obj.meshUnion.x, obj.meshUnion.y, obj.meshUnion.z, time) .* obj.hy(obj.meshUnion.x, obj.meshUnion.y, time) + ...
            obj.OH(obj.meshUnion.x, obj.meshUnion.y, obj.meshUnion.z, time) .* obj.v2dy(obj.meshUnion.x, obj.meshUnion.y, time) - ...
            obj.Ov(obj.meshUnion.x, obj.meshUnion.y, obj.meshUnion.z, time) .* obj.hy(obj.meshUnion.x, obj.meshUnion.y, time) - ...
            obj.h(obj.meshUnion.x, obj.meshUnion.y, time) .* obj.Ovy(obj.meshUnion.x, obj.meshUnion.y, obj.meshUnion.z, time);
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
                fphys2d{m} = zeros( mesh2d.cell.Np, mesh2d.K, obj.Nfield2d );
                fphys{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield );
                % bottom elevation
                fphys2d{m}(:, :, 4) = -(2-0.005*(mesh2d.x+mesh2d.y));
                
                fphys2d{m}(:,:,5) = ( mesh2d.rx .* ( mesh2d.cell.Dr * fphys2d{m}(:,:,4) ) ) + ...
                    ( mesh2d.sx .* ( mesh2d.cell.Ds * fphys2d{m}(:,:,4) ) );
                
                fphys2d{m}(:,:,6) = ( mesh2d.ry .* ( mesh2d.cell.Dr * fphys2d{m}(:,:,4) ) ) + ...
                    ( mesh2d.sy .* ( mesh2d.cell.Ds * fphys2d{m}(:,:,4) ) );
                %water depth
                fphys2d{m}(:,:,1) = obj.h( mesh2d.x, mesh2d.y, 0 );
                fphys{m}(:,:,4) = obj.h( mesh3d.x, mesh3d.y, 0 );
                fphys{m}(:,:,1) = obj.h( mesh3d.x, mesh3d.y, 0 ) .* obj.u( mesh3d.x, mesh3d.y, mesh3d.z, 0 );
                fphys{m}(:,:,2) = obj.h( mesh3d.x, mesh3d.y, 0 ) .* obj.v( mesh3d.x, mesh3d.y, mesh3d.z, 0 );
                fphys{m}(:,:,3) = obj.Ou2d( mesh3d.x, mesh3d.y, mesh3d.z, zeros(size(mesh3d.z))) .* obj.hx( mesh3d.x, mesh3d.y, 0) + ...
                    obj.OH( mesh3d.x, mesh3d.y, mesh3d.z, 0) .* obj.u2dx( mesh3d.x, mesh3d.y, 0) - ...
                    obj.Ou( mesh3d.x, mesh3d.y, mesh3d.z, 0) .* obj.hx( mesh3d.x, mesh3d.y, 0) - ...
                    obj.h( mesh3d.x, mesh3d.y, 0) .* obj.Oux( mesh3d.x, mesh3d.y, mesh3d.z, 0) + ...
                    obj.Ov2d( mesh3d.x, mesh3d.y, mesh3d.z, 0) .* obj.hy( mesh3d.x, mesh3d.y, 0) + ...
                    obj.OH( mesh3d.x, mesh3d.y, mesh3d.z, 0) .* obj.v2dy( mesh3d.x, mesh3d.y, 0) - ...
                    obj.Ov( mesh3d.x, mesh3d.y, mesh3d.z, 0) .* obj.hy( mesh3d.x, mesh3d.y, 0) - ...
                    obj.h( mesh3d.x, mesh3d.y, 0) .* obj.Ovy( mesh3d.x, mesh3d.y, mesh3d.z, 0);
%                 [ fphys{m}(:,:,4), fphys{m}(:,:,1), fphys{m}(:,:,2), fphys{m}(:,:,3), ~ ] = obj.matGetExactSolution( mesh3d.x, mesh3d.y, mesh3d.z, 0);
                
                %> Z in extended three dimensional fields
                fphys{m}(:,:,6) = mesh3d.Extend2dField( fphys2d{m}(:, :, 4) );
                %> '$\eta$' in extended three dimensional fields
                fphys{m}(:,:,7) = fphys{m}(:,:,4) + fphys{m}(:,:,6);
                
                fphys2d{m}(:, :, 2) = obj.u2d(mesh2d.x,mesh2d.y,0).*obj.h(mesh2d.x,mesh2d.y,0);
                fphys2d{m}(:, :, 3) = obj.v2d(mesh2d.x,mesh2d.y,0).*obj.h(mesh2d.x,mesh2d.y,0);
                
            end
        end
        
        function  matEvaluateSourceTerm(obj, fphys, time)
            mesh = obj.meshUnion(1);
            mesh2d = obj.mesh2d(1);
            matEvaluateSourceTerm@SWEAbstract3d( obj, fphys);
            %[Problematic
            obj.frhs2d{1}(:,:,1) = obj.frhs2d{1}(:,:,1) + obj.ht(mesh2d.x, mesh2d.y, time) + ...
                obj.h(mesh2d.x, mesh2d.y, time) .* obj.u2dx(mesh2d.x, mesh2d.y, time) + ...
                obj.u2d(mesh2d.x, mesh2d.y, time) .* obj.hx(mesh2d.x, mesh2d.y, time) + ...
                obj.h(mesh2d.x, mesh2d.y, time) .* obj.v2dy(mesh2d.x, mesh2d.y, time) + ...
                obj.v2d(mesh2d.x, mesh2d.y, time) .* obj.hy(mesh2d.x, mesh2d.y, time) ;
            
            
            Omega = obj.Ou2d( mesh.x, mesh.y, mesh.z, time) .* obj.hx( mesh.x, mesh.y, time) + ...
                    obj.OH( mesh.x, mesh.y, mesh.z, time) .* obj.u2dx( mesh.x, mesh.y, time) - ...
                    obj.Ou( mesh.x, mesh.y, mesh.z, time) .* obj.hx( mesh.x, mesh.y, time) - ...
                    obj.h( mesh.x, mesh.y, time) .* obj.Oux( mesh.x, mesh.y, mesh.z, time) + ...
                    obj.Ov2d( mesh.x, mesh.y, mesh.z, time) .* obj.hy( mesh.x, mesh.y, time) + ...
                    obj.OH( mesh.x, mesh.y, mesh.z, time) .* obj.v2dy( mesh.x, mesh.y, time)- ...
                    obj.Ov( mesh.x, mesh.y, mesh.z, time) .* obj.hy( mesh.x, mesh.y, time) - ...
                    obj.h( mesh.x, mesh.y, time) .* obj.Ovy( mesh.x, mesh.y, mesh.z, time);
            
            Os = obj.u2d( mesh.x, mesh.y, time) .* obj.hx( mesh.x, mesh.y, time) + ...
                obj.h( mesh.x, mesh.y, time).* obj.u2dx( mesh.x, mesh.y, time)  - ...
                obj.u( mesh.x, mesh.y, mesh.z, time) .* obj.hx( mesh.x, mesh.y, time) - ...
                obj.h( mesh.x, mesh.y, time).* obj.ux( mesh.x, mesh.y, mesh.z, time)  + ...
                obj.v2d( mesh.x, mesh.y, time) .* obj.hy( mesh.x, mesh.y, time) + ...
                obj.h( mesh.x, mesh.y, time).* obj.v2dy( mesh.x, mesh.y, time)  - ...
                obj.v( mesh.x, mesh.y, mesh.z, time) .* obj.hy( mesh.x, mesh.y, time) - ...
                obj.h( mesh.x, mesh.y, time).* obj.vy( mesh.x, mesh.y, mesh.z, time);
            
            
            obj.frhs{1}(:,:,1) = obj.frhs{1}(:,:,1) + obj.u( mesh.x, mesh.y, mesh.z, time ) .* obj.ht( mesh.x, mesh.y, time ) +...
                obj.h( mesh.x, mesh.y, time ) .* obj.ut( mesh.x, mesh.y, mesh.z, time ) + obj.u( mesh.x, mesh.y, mesh.z, time ).*...
                (obj.u( mesh.x, mesh.y, mesh.z, time ).*obj.hx( mesh.x, mesh.y, time ) + obj.h( mesh.x, mesh.y, time ).*obj.ux( mesh.x, mesh.y, mesh.z, time ))...
                + obj.h( mesh.x, mesh.y, time ).*obj.u( mesh.x, mesh.y, mesh.z, time ).*obj.ux( mesh.x, mesh.y, mesh.z, time ) + ...
                obj.gra .* obj.h( mesh.x, mesh.y, time ) .* obj.hx( mesh.x, mesh.y, time ) - obj.gra .* fphys{1}(:,:,6).*fphys{1}(:,:,8) +...
                obj.u( mesh.x, mesh.y, mesh.z, time ) .* (obj.v( mesh.x, mesh.y, mesh.z, time ).*obj.hy( mesh.x, mesh.y, time ) + ...
                obj.h( mesh.x, mesh.y, time ) .* obj.vy( mesh.x, mesh.y, mesh.z, time )) + obj.u( mesh.x, mesh.y, mesh.z, time ).*Os + ...
                Omega .* obj.uz( mesh.x, mesh.y, time ) + obj.gra .* obj.eta( mesh.x, mesh.y, time ) .* fphys{1}(:,:,8);
            

            obj.frhs{1}(:,:,2) = obj.frhs{1}(:,:,2) + obj.v( mesh.x, mesh.y, mesh.z, time ).*obj.ht( mesh.x, mesh.y, time ) + ...
                obj.h( mesh.x, mesh.y, time ).*obj.vt( mesh.x, mesh.y, mesh.z, time ) + obj.v( mesh.x, mesh.y, mesh.z, time ).* ...
                (obj.u( mesh.x, mesh.y, mesh.z, time ).*obj.hx( mesh.x, mesh.y, time ) + obj.h( mesh.x, mesh.y, time ).*obj.ux( mesh.x, mesh.y, mesh.z, time ))...
                + obj.v( mesh.x, mesh.y, mesh.z, time ).*(obj.v( mesh.x, mesh.y, mesh.z, time ).*obj.hy( mesh.x, mesh.y, time ) + obj.h( mesh.x, mesh.y, time ).* ...
                obj.vy( mesh.x, mesh.y, mesh.z, time )) + obj.h( mesh.x, mesh.y, time ).*obj.v( mesh.x, mesh.y, mesh.z, time ).*obj.vy( mesh.x, mesh.y, mesh.z, time ) + ...
                obj.gra.*obj.h( mesh.x, mesh.y, time ).*obj.hy( mesh.x, mesh.y, time ) - obj.gra .* fphys{1}(:,:,6).*fphys{1}(:,:,9) + obj.v( mesh.x, mesh.y, mesh.z, time ).*Os + ...
                Omega .* obj.vz( mesh.x, mesh.y, time ) + obj.gra .* obj.eta( mesh.x, mesh.y, time ) .* fphys{1}(:,:,9);
        end
        
        function matUpdateExternalField( obj, time, ~, ~ ) 
            
            obj.fext3d{1}( :, :, 1 ) = obj.h(obj.mesh3d.BoundaryEdge.xb, obj.mesh3d.BoundaryEdge.yb, time) .* ...
                obj.u(obj.mesh3d.BoundaryEdge.xb, obj.mesh3d.BoundaryEdge.yb, obj.mesh3d.BoundaryEdge.zb, time);
            obj.fext3d{1}( :, :, 2 ) = obj.h(obj.mesh3d.BoundaryEdge.xb, obj.mesh3d.BoundaryEdge.yb, time) .* ...
                obj.v(obj.mesh3d.BoundaryEdge.xb, obj.mesh3d.BoundaryEdge.yb, obj.mesh3d.BoundaryEdge.zb, time);
            obj.fext3d{1}( :, :, 4 ) = obj.h(obj.mesh3d.BoundaryEdge.xb, obj.mesh3d.BoundaryEdge.yb, time);
            obj.fext2d{1}( :, :, 1 ) = obj.h(obj.mesh2d.BoundaryEdge.xb, obj.mesh2d.BoundaryEdge.yb, time);
                        
            obj.fext2d{1}( :, :, 2 ) = obj.meshUnion.BoundaryEdge.VerticalColumnIntegralField(  obj.fext3d{1}( :, :, 1 ) );      
            obj.fext2d{1}( :, :, 3 ) = obj.meshUnion.BoundaryEdge.VerticalColumnIntegralField(  obj.fext3d{1}( :, :, 2 ) );        
            
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
zs = zeros(mesh2d.Nv, 1); zb = zeros(mesh2d.Nv, 1) - ones(mesh2d.Nv, 1);
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
end