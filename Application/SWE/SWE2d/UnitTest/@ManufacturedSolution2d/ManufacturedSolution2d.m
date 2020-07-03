classdef ManufacturedSolution2d < SWEPreBlanaced2d
    %MANUFACTUREDSOLUTION2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        PostProcess
        timePoint
        ErrNorm2
        Index
        ExactValue
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
        function obj = ManufacturedSolution2d(N, deltax, cellType)
            [ mesh ] = makeUniformMesh(N, deltax, cellType);
            obj = obj@SWEPreBlanaced2d();
            obj.hmin = 1e-3;
            obj.initPhysFromOptions( mesh );
            %> time interval
            obj.PostProcess = NdgPostProcess(obj.meshUnion(1),...
                strcat('WindDrivenFlow/3d','/','WindDrivenFlow'));
            obj.ErrNorm2 = cell(3);
            obj.Index = 1; 
            obj.ExactValue = cell(1);
            [obj.ExactValue{1}(:,:,1), obj.ExactValue{1}(:,:,2), obj.ExactValue{1}(:,:,3)] = ...
               obj.matGetExactSolution( obj.meshUnion.x, obj.meshUnion.y, obj.ftime); 
        end
        
        function matPlotErrorNorm2(obj)
            close all
            figure(1);
            hold on;
            color = {'r-','k-','b-'};
            LWidth = 1.5;
            FontSize = 15;
            for i = 1:3
                plot(obj.timePoint, obj.ErrNorm2{i},color{i},'LineWidth',LWidth);
            end
            lendstr = {'$h$',...
                '$hu$','$hv$'};
            legend(lendstr,'Interpreter','Latex','FontSize',FontSize);
            xlabel('$time(s)$','Interpreter','Latex');
            ylabel('$L_2 error$','Interpreter','Latex');
            box on;
        end
        
    end
    methods ( Access = protected )
        
        matCalculateExplicitRHSTerm( obj, fphys2d, fphys, Stage, RKIndex, time);
        
        function matEvaluateError( obj, fphys, time)
%             [ h, hu, hv ] = obj.matGetExactSolution( obj.meshUnion.x, obj.meshUnion.y, time);
%             fext = cell(1);
%             fext{1}(:,:,1) = h;
%             fext{1}(:,:,2) = hu;
%             fext{1}(:,:,3) = hv;
%             Tempfphys = cell(1);
%             Tempfphys{1}(:,:,1) = fphys(:,:,1);
%             Tempfphys{1}(:,:,2) = fphys(:,:,2);
%             Tempfphys{1}(:,:,3) = fphys(:,:,3);
%             Err2 = obj.PostProcess.evaluateNormErr2( Tempfphys, fext );
%             obj.timePoint(obj.Index) = time;
%             obj.ErrNorm2{1}(obj.Index)  = Err2(1);
%             obj.ErrNorm2{2}(obj.Index)  = Err2(2);
%             obj.ErrNorm2{3}(obj.Index)  = Err2(3);
%             obj.Index = obj.Index + 1;
        end
        
        %> set initial function
        function [ fphys] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                % bottom elevation
                fphys{m}(:, :, 4) = -(2-0.005*(mesh.x+mesh.y));
                
                fphys{m}(:,:,5) = ( mesh.rx .* ( mesh.cell.Dr * fphys{m}(:,:,4) ) ) + ...
                    ( mesh.sx .* ( mesh.cell.Ds * fphys{m}(:,:,4) ) );
                
                fphys{m}(:,:,6) = ( mesh.ry .* ( mesh.cell.Dr * fphys{m}(:,:,4) ) ) + ...
                    ( mesh.sy .* ( mesh.cell.Ds * fphys{m}(:,:,4) ) );
                %water depth
                fphys{m}(:,:,1) = obj.e*(sin(obj.w*(mesh.x+0)) + sin(obj.w*(mesh.y+0))) - ...
                    fphys{m}(:, :, 4);
                fphys{m}(:,:,2) = sin(obj.w*(mesh.x+0)) .* fphys{m}(:, :, 1);        
                fphys{m}(:,:,3) = sin(obj.w*(mesh.y+0)) .* fphys{m}(:, :, 1);                  
            end
        end
        
        
        function matUpdateExternalField( obj, time, ~, ~ )
            
            [obj.fext{1}( :, :, 1 ), obj.fext{1}( :, :, 2 ), obj.fext{1}( :, :, 3 )] = obj.matGetExactSolution( obj.meshUnion.BoundaryEdge.xb, obj.meshUnion.BoundaryEdge.yb, time);
            
        end
        
        function [ h, hu, hv ] = matGetExactSolution( obj, x, y, time)
            h = obj.e*(sin(obj.w.*(x+time)) + sin(obj.w*(y+time))) + (2-0.005.*(x + y));
            hu = sin(obj.w.*(x+time)) .* h;
            hv = sin(obj.w.*(y+time)) .* h;
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
    
    methods ( Access = protected )
        
        function  matEvaluateSourceTerm(obj, fphys, time)
            mesh = obj.meshUnion(1);
            matEvaluateSourceTerm@SWEAbstract2d( obj, fphys);
            h2d =  obj.e * ( sin(obj.w*(mesh.x+time)) + sin(obj.w*(mesh.y+time)) ) + 2 - 0.005*( mesh.x + mesh.y );
            h2dt = obj.e * obj.w .* cos(obj.w*(mesh.x+time)) + obj.e * obj.w .* cos(obj.w*(mesh.y+time));
            u2d = sin(obj.w*(mesh.x+time));
            u2dx = obj.w * cos(obj.w*(mesh.x + time));
            u2dt = obj.w * cos(obj.w*(mesh.x + time));
            
            v2d = sin(obj.w*(mesh.y+time));
            v2dy = obj.w * cos(obj.w*(mesh.y + time));
            v2dt = obj.w * cos(obj.w*(mesh.y + time));
            h2dx = obj.e * obj.w *( cos(obj.w*(mesh.x+time)) ) - 0.005;
            h2dy = obj.e * obj.w *( cos(obj.w*(mesh.y+time)) ) - 0.005;
            eta = obj.e*(sin(obj.w.*(mesh.x+time)) + sin(obj.w*(mesh.y+time)));
            
            obj.frhs{1}(:,:,1) = obj.frhs{1}(:,:,1) +  h2dt  + ...
                h2d .* u2dx + u2d.*h2dx + h2d .* v2dy + v2d .* h2dy;
            
            
            obj.frhs{1}(:,:,2) = obj.frhs{1}(:,:,2) + u2d.*h2dt + h2d.*u2dt + u2d.*(u2d.*h2dx+h2d.*u2dx) + h2d.*u2d.*u2dx + obj.gra.*h2d.*h2dx - ...
                obj.gra .* fphys{1}(:,:,4).*fphys{1}(:,:,5) + u2d .* (v2d.*h2dy + h2d.*v2dy) + obj.gra .* eta .* fphys{1}(:,:,5);
            obj.frhs{1}(:,:,3) = obj.frhs{1}(:,:,3) + v2d.*h2dt + h2d.*v2dt + v2d.*(u2d.*h2dx+h2d.*u2dx) + v2d.*(v2d.*h2dy + h2d.*v2dy) + h2d.*v2d.*v2dy + ...
                obj.gra.*h2d.*h2dy - obj.gra .* fphys{1}(:,:,4).*fphys{1}(:,:,6) + obj.gra .* eta .* fphys{1}(:,:,6);
        end
    end
    
end

function [ mesh ] = makeUniformMesh(N, deltax, type)
bctype = [...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped];

if (type == enumStdCell.Tri)
    mesh = makeUniformTriMesh(N, [0, 100], [0, 100], 100/deltax, 100/deltax, bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N,[0, 100], [0, 100], 100/deltax, 100/deltax, bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func