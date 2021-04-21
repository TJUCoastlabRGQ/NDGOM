classdef UniformManufacturedSolution3d < ManufacturedSolution3d
    %MANUFACTUREDSOLUTION3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    methods
        function obj = UniformManufacturedSolution3d( N, Nz, M, Mz )
            obj = obj@ManufacturedSolution3d( N, Nz, M, Mz );
        end
    end
    
    methods ( Access = protected )
        
        function matGetFunction(obj)
            syms x y z t;
            obj.eta = ( obj.e * ( sin(obj.w*(x+t)) + sin(obj.w*(y+t)) ) );
            obj.b = - ( 2 - 0.005*( x + y ));
            obj.h = obj.eta - obj.b;
            obj.u = sin(obj.w.*(x+t));
            obj.v = sin(obj.w.*(y+t));
            obj.ht = diff(obj.h,t);
            obj.u2d = int( obj.u, z, [-1,0] );
            obj.v2d = int( obj.v, z, [-1,0] );
            obj.Omega = int(-obj.ht - diff(obj.h*obj.u, x) - diff(obj.h*obj.v, y), z, [-1,z]);
            obj.Source2d = obj.ht + diff(obj.h*obj.u2d, x) + diff(obj.h*obj.v2d, y);

            obj.hut = diff( obj.h* obj.u, t);
            obj.mhux = diff( obj.h * obj.u * obj.u + 0.5 * obj.gra * ( obj.h * obj.h - obj.b * obj.b), x);
            obj.mhuy = diff( obj.h * obj.u * obj.v, y);
            obj.mhuz = diff( obj.u * obj.Omega, z);
            obj.hvt = diff( obj.h* obj.v, t);
            obj.mhvx = diff( obj.h * obj.u * obj.v, x);
            obj.mhvy = diff( obj.h * obj.v * obj.v + 0.5 * obj.gra * ( obj.h * obj.h - obj.b * obj.b), y);
            obj.mhvz = diff( obj.v * obj.Omega, z);
            obj.mh2dx = diff( obj.h * obj.u2d, x);
            obj.mh2dy = diff( obj.h * obj.v2d, y);
            obj.mhuzz = diff( obj.miu/obj.h/obj.h * diff( obj.h*obj.u , z ), z );
            obj.mhvzz = diff( obj.miu/obj.h/obj.h * diff( obj.h*obj.v , z ), z );
            obj.hufz = diff( obj.h*obj.u , z );
            obj.hvfz = diff( obj.h*obj.v , z );
        end
        
        function PostInit(obj)
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder3d = [ 1 2 3];
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            obj.PostProcess = NdgPostProcess(obj.meshUnion(1),...
                strcat('UniformManufacturedSolution3d/','3d/','UniformManufacturedSolution3d'));
            obj.lendstr = {'$hu$',...
                '$hv$','$h$'};
            obj.ErrNorm2 = cell(3,1);
            obj.Index = 1;
            obj.ExactValue = cell(1);
            [obj.ExactValue{1}(:,:,1), obj.ExactValue{1}(:,:,2), ~, obj.ExactValue{1}(:,:,3)] = ...
                obj.matGetExactSolution( obj.meshUnion.x, obj.meshUnion.y, obj.meshUnion.z, obj.ftime);
        end
        
        function matEvaluateError( obj, fphys, time)
            [ hu, hv, ~, h ] = obj.matGetExactSolution( obj.mesh3d.x, obj.mesh3d.y, obj.mesh3d.z, time);
            fext = cell(1);
            fext{1}(:,:,1) = hu;
            fext{1}(:,:,2) = hv;
            fext{1}(:,:,3) = h;
            Tempfphys = cell(1);
            Tempfphys{1}(:,:,1) = fphys(:,:,1);
            Tempfphys{1}(:,:,2) = fphys(:,:,2);
            Tempfphys{1}(:,:,3) = fphys(:,:,4);
            Err2 = obj.PostProcess.evaluateNormErr2( Tempfphys, fext );
            obj.timePoint(obj.Index) = time;
            obj.ErrNorm2{1}(obj.Index)  = Err2(1);
            obj.ErrNorm2{2}(obj.Index)  = Err2(2);
            obj.ErrNorm2{3}(obj.Index)  = Err2(3);
            obj.Index = obj.Index + 1;
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