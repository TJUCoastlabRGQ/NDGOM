classdef NdgNonhydrostaticAbstractTest < SWEConventional2d
    %NDGNONHYDROSTATICABSTRACTTEST 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
        gra = 9.81
        hmin = 1e-3
        rho = 1000
        tol = 10^(-9)
    end
    
    methods
        function obj = NdgNonhydrostaticAbstractTest(N, cellType)
            obj = obj@SWEConventional2d();
            [ mesh ] = obj.makeUniformMesh(N, cellType);
            obj.initPhysFromOptions( mesh );
        end
        
        function testmatEvaluateNonhydrostaticSurfaceValue(obj)
            [ExactFm, ExactFp] = obj.getExactNonhydroSurfaceValue;
            [fm, fp] = obj.NonhydrostaticSolver.matEvaluateNonhydrostaticSurfaceValue(obj.meshUnion(1), ones(size(obj.meshUnion(1).x)));
            obj.Assert(ExactFm, fm);
            obj.Assert(ExactFp, fp);
        end
        
%         function testmatEvaluatePenaltyTerm(obj)
%             [ExactPenaltyX, ExactPenaltyY] = obj.getExactPenaltyTerm;
%             [fm, fp] = obj.NonhydrostaticSolver.matEvaluateNonhydrostaticSurfaceValue(obj.meshUnion(1), ones(size(obj.meshUnion(1).x)));
%             PenaltyX = obj.NonhydrostaticSolver.matEvaluatePenaltyTerm(obj.meshUnion(1), fm, fp, obj.meshUnion(1).nx);
%             PenaltyY = obj.NonhydrostaticSolver.matEvaluatePenaltyTerm(obj.meshUnion(1), fm, fp, obj.meshUnion(1).ny);
%             obj.Assert(ExactPenaltyX, PenaltyX);
%             obj.Assert(ExactPenaltyY, PenaltyY);
%         end
        
%         function testmatEvaluateVelocityLikeTerms(obj)
%             ExactharmonicTerm = obj.getExactHarmonicTerm;
%             [fm, fp] = obj.NonhydrostaticSolver.matEvaluateNonhydrostaticSurfaceValue(obj.meshUnion(1), ones(size(obj.meshUnion(1).x)));
%             harmonicTerm = obj.NonhydrostaticSolver.matEvaluateVelocityLikeTerms( obj.meshUnion(1), fm, fp );
%             obj.Assert(ExactharmonicTerm, harmonicTerm);
%         end
        
%         function testmatEvaluateNonconservativeFlux(obj)
%             [ ExactQx, ExactQy ]= obj.getExactNonconservativeFlux;
%             [fm, fp] = obj.NonhydrostaticSolver.matEvaluateNonhydrostaticSurfaceValue(obj.meshUnion(1), ones(size(obj.meshUnion(1).x)));
%             harmonicH = obj.NonhydrostaticSolver.matEvaluateVelocityLikeTerms( obj.meshUnion(1), fm, fp );
%             PenaltyX = obj.NonhydrostaticSolver.matEvaluatePenaltyTerm(obj.meshUnion(1), fm, fp, obj.meshUnion(1).nx);
%             PenaltyY = obj.NonhydrostaticSolver.matEvaluatePenaltyTerm(obj.meshUnion(1), fm, fp, obj.meshUnion(1).ny);
%             qx = obj.NonhydrostaticSolver.matEvaluateNonconservativeFlux( harmonicH, fm, fp, PenaltyX);
%             qy = obj.NonhydrostaticSolver.matEvaluateNonconservativeFlux( harmonicH, fm, fp, PenaltyY);
%             obj.Assert(ExactQx, qx);
%             obj.Assert(ExactQy, qy);
%         end
        
%         function testmatEvaluateDeltaSurfaceFlux(obj)
%             ExactDeltaFlux = obj.getExactDeltaFlux;
%             [fm, fp] = obj.NonhydrostaticSolver.matEvaluateNonhydrostaticSurfaceValue(obj.meshUnion(1), ones(size(obj.meshUnion(1).x)));
%             harmonicH = obj.NonhydrostaticSolver.matEvaluateVelocityLikeTerms( obj.meshUnion(1), fm, fp );
%             PenaltyX = obj.NonhydrostaticSolver.matEvaluatePenaltyTerm(obj.meshUnion(1), fm, fp, obj.meshUnion(1).nx);
%             PenaltyY = obj.NonhydrostaticSolver.matEvaluatePenaltyTerm(obj.meshUnion(1), fm, fp, obj.meshUnion(1).ny);
%             qx = obj.NonhydrostaticSolver.matEvaluateNonconservativeFlux( harmonicH, fm, fp, PenaltyX);
%             qy = obj.NonhydrostaticSolver.matEvaluateNonconservativeFlux( harmonicH, fm, fp, PenaltyY);
%             DeltaFlux = obj.NonhydrostaticSolver.matEvaluateDeltaSurfaceFlux( obj.meshUnion(1), ones(size(obj.meshUnion(1).x)), ones(size(obj.meshUnion(1).x)), qx, qy );
%             obj.Assert(ExactDeltaFlux, DeltaFlux);
%         end
        
%         function testmatEvaluateDeltaSurfaceFluxY(obj)
%             ExactDeltaFlux = obj.getExactDeltaFluxY;
%             [fm, fp] = obj.NonhydrostaticSolver.matEvaluateNonhydrostaticSurfaceValue(obj.meshUnion(1), ones(size(obj.meshUnion(1).x)));
%             PenaltyY = zeros(size(obj.meshUnion(1).eidM));
%             qY = obj.NonhydrostaticSolver.matEvaluateNonconservativeFlux( fm, fp, PenaltyY);
%             DeltaFlux = obj.NonhydrostaticSolver.matEvaluateDeltaSurfaceFluxY(obj.meshUnion(1), ones(size(obj.meshUnion(1).eidM)), qY);
%             obj.Assert(ExactDeltaFlux, DeltaFlux);
%         end
        
%         function testmatEvaluateDeltaSurfaceFluxX(obj)
%             ExactDeltaFlux = obj.getExactDeltaFluxX;
%             [fm, fp] = obj.NonhydrostaticSolver.matEvaluateNonhydrostaticSurfaceValue(obj.meshUnion(1), ones(size(obj.meshUnion(1).x)));
%             PenaltyX = zeros(size(obj.meshUnion(1).eidM));
%             qX = obj.NonhydrostaticSolver.matEvaluateNonconservativeFlux( fm, fp, PenaltyX);
%             DeltaFlux = obj.NonhydrostaticSolver.matEvaluateDeltaSurfaceFluxX(obj.meshUnion(1), ones(size(obj.meshUnion(1).eidM)), qX);
%             obj.Assert(ExactDeltaFlux, DeltaFlux);
%         end
        
%         function testmatEvaluateDivergenceDirectionX(obj)
%             ExactDivergenceX = obj.getExactDivergenceX;
%             DivergenceX = obj.NonhydrostaticSolver.matEvaluateDivergenceDirectionX(obj.meshUnion(1), ones(size(obj.meshUnion(1).x)));
%             obj.Assert(ExactDivergenceX, DivergenceX);
%         end
        
%         function testmatEvaluateDivergenceDirectionY(obj)
%             ExactDivergenceY = obj.getExactDivergenceY;
%             DivergenceY = obj.NonhydrostaticSolver.matEvaluateDivergenceDirectionY(obj.meshUnion(1), ones(size(obj.meshUnion(1).x)));
%             obj.Assert(ExactDivergenceY, DivergenceY);
%         end
        
        function testAssembleWetDryInterface(obj)
            ExactWetCellIndex = obj.getExactWetCellIndex;
            obj.NonhydrostaticSolver.getWetDryInterface(obj.meshUnion(1));
            obj.Assert(ExactWetCellIndex, obj.NonhydrostaticSolver.WetCellIndex);
            [ExactZeroFluxBoundary, ExactZeroFluxBoundaryIndex] = obj.getExactZeroFluxBoundary;
            obj.Assert(ExactZeroFluxBoundary, obj.NonhydrostaticSolver.ZeroFluxBoundary);
            obj.Assert(ExactZeroFluxBoundaryIndex, obj.NonhydrostaticSolver.ZeroFluxBoundaryIndex);
            ExactAdjacentDryCellAndFace = obj.getExactAdjacentDryCellAndFace;
            obj.Assert(ExactAdjacentDryCellAndFace, obj.NonhydrostaticSolver.AdjacentDryCellAndFace);
        end
        
        function testEidBoundaryType(obj)
            %> @brief function for testing the boundary type related non-hydrostatic information matrix
            %> EidBoundaryType is the flag matrix used to impose the boundary related Non-hydrostatic condition, at both the slip wall and non-slip wall boundaries, this value is set to 1 to indicate
            %> the zero gradient non-hydrostatic pressure condition($q^+ = q^-$) and the zero non-hydrostatic gradient condition($\frac{\partial q^+}{\partial x} =-\frac{\partial q^-}{\partial x} $), while at
            %> the clamped type boundaries, this value is set to -1 to indicate the zero non-hydrostatic pressure condition($q^+ = -q^-$) and the zero gradient non-hydrostatic gradient condition
            %> ($\frac{\partial q^+}{\partial x} =\frac{\partial q^-}{\partial x} $)
            ExactEidBoundaryType = obj.getExactEidBoundaryType;
            obj.Assert(ExactEidBoundaryType, obj.NonhydrostaticSolver.EidBoundaryType);
        end
        
        function testmatGetFaceValue(obj)
            %> @brief function for testing the Nonhydrostatic pressure at face point
            mesh = obj.meshUnion(1);
            obj.NonhydrostaticSolver.getWetDryInterface(mesh);
            %to get the NonhydroFmPoint and NonhydroFpPoint information
            [~, ~, ~, ~, ~, ~]= obj.NonhydrostaticSolver.matReconstructStiffmatrixRelatedMatrix( obj); 
            [fm, fp] = mesh.InnerEdge.matEvaluateSurfValue( obj.fphys );
            % for this case, the water depth is treated as the non-hydrostatic pressure to test the inner and outer value
            [fm, fp] = obj.NonhydrostaticSolver.GetFaceValue(fm(:,:,1), fp(:,:,1), enumNonhydroBoundaryCondition.Zero); 
            [ExactFm, ExactFp] = obj.getExactInnerOuterValue;
            obj.Assert(ExactFm, fm);
            obj.Assert(ExactFp, fp);
        end
        
        function testmatImposeNonhydroRelatedBoundaryCondition(obj)
            [ ExactZeroFp, ExactZeroGradFp ] = obj.getExactBoundaryFp;
             % for this case, the water depth is treated as the non-hydrostatic pressure to test the boundary value
            [Fm, Fp] = obj.meshUnion(1).BoundaryEdge.matEvaluateSurfValue(num2cell(obj.fphys{1}(:,:,1),[1 2]));  %
            ZeroFp = obj.NonhydrostaticSolver.GetBoundaryValue(Fm, Fp, enumNonhydroBoundaryCondition.Zero);
            ZeroGradFp = obj.NonhydrostaticSolver.GetBoundaryValue(Fm, Fp, enumNonhydroBoundaryCondition.ZeroGrad);
            obj.Assert(ExactZeroFp, ZeroFp);
            obj.Assert(ExactZeroGradFp, ZeroGradFp);
        end
        

        
        
        function testGlobalMatrixAssemble(obj)
            %> @brief function for testing the global matrix assemble methods
            %> In this study, the global matrix is assembled at the initial
            %> stage of the computation so as to reduce the computational cost.
            %> To do that, a matrix was constructed at first. Then, at each
            %> step, the matrix was constructed by multiplying each column of the matrix by
            %> the water depth at the studied point or the whole depth field
            mesh = obj.meshUnion(1);
            K = mesh.K; Np = mesh.cell.Np;
            InnerEdge = mesh.InnerEdge;
            BoundaryEdge = mesh.BoundaryEdge;
            height = obj.fphys{1}(:,:,1);

            for index = 1:K*Np
                gmat = zeros( Np, K );
                gmat(index) = 1.0/2.0;
                [tempqx, tempqy] = obj.NonhydrostaticSolver.matCalculateLDGAuxialaryTerm...
                    ( mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]), enumNonhydroBoundaryCondition.Zero);
                [Exactqx, Exactqy] = obj.NonhydrostaticSolver.matCalculateLDGAuxialaryTerm...
                    ( mesh, BoundaryEdge, InnerEdge, num2cell(height(index) .* gmat,[1 2]), enumNonhydroBoundaryCondition.Zero);
                [tempq2x, tempq2y] = obj.NonhydrostaticSolver.matCalculateLDGTerm( mesh, ...
                    BoundaryEdge, InnerEdge, num2cell( gmat,[1 2]), num2cell(tempqx,[1 2]), num2cell(tempqy,[1 2]), enumNonhydroBoundaryCondition.ZeroGrad);
                [Exactq2x, Exactq2y] = obj.NonhydrostaticSolver.matCalculateLDGTerm( mesh,...
                    BoundaryEdge, InnerEdge, num2cell(height(index) .* gmat,[1 2]), num2cell(Exactqx,[1 2]), num2cell(Exactqy,[1 2]), enumNonhydroBoundaryCondition.ZeroGrad);
                obj.Assert( height(index) .* tempqx, Exactqx );
                obj.Assert( height(index) .* tempqy, Exactqy );
                obj.Assert( height(index) .* tempq2x, Exactq2x );
                obj.Assert( height(index) .* tempq2y, Exactq2y );
            end   
            display('===========================================');
            display('LDG Assemble verified');
            display('===========================================');
            
            for index = 1:K*Np
                gmat = zeros( Np, K );
                gmat(index) = 1.0/2.0;
                [ tempqx, tempqy ] = obj.NonhydrostaticSolver.GetCharacteristicMatrix...
                    ( mesh, num2cell(gmat,[1 2]), num2cell(gmat,[1 2]), enumNonhydroBoundaryCondition.Zero);
                [ Exactqx, Exactqy ] = obj.NonhydrostaticSolver.GetCharacteristicMatrix...
                    ( mesh, num2cell(height(index) .* gmat,[1 2]), num2cell(height(index) .* gmat,[1 2]), enumNonhydroBoundaryCondition.Zero);
                
                [ tempq2x, tempq2y ] = obj.NonhydrostaticSolver.GetCharacteristicMatrix...
                    ( mesh, num2cell(tempqx,[1 2]), num2cell(tempqy,[1 2]), enumNonhydroBoundaryCondition.ZeroGrad);
                [ Exactq2x, Exactq2y ] = obj.NonhydrostaticSolver.GetCharacteristicMatrix...
                    ( mesh, num2cell(Exactqx,[1 2]), num2cell(Exactqy,[1 2]), enumNonhydroBoundaryCondition.ZeroGrad);
                          
                obj.Assert( height(index) .* tempqx, Exactqx );
                obj.Assert( height(index) .* tempqy, Exactqy );
                obj.Assert( height(index) .* tempq2x, Exactq2x );
                obj.Assert( height(index) .* tempq2y, Exactq2y );
            end   
            display('===========================================');
            display('Central Assemble verified');
            display('===========================================');      
            
            for index = 1:K*Np
                ele = ceil(index/Np);
                gmat = zeros( Np, K );
                gmat(index) = 1.0/2.0;
                [UpWindedFlag, DownWindedFlag] = AssembleWindedFlagInformation(InnerEdge, ele);
                
                [tempqx, tempqy] = obj.NonhydrostaticSolver.matCalculateFluxUpwindedTerm...
                    ( mesh, BoundaryEdge, InnerEdge, num2cell( gmat,[1 2]),UpWindedFlag, DownWindedFlag, enumNonhydroBoundaryCondition.Zero);
                [ Exactqx, Exactqy ] = obj.NonhydrostaticSolver.matCalculateFluxUpwindedTerm...
                    ( mesh, BoundaryEdge, InnerEdge, num2cell(height(index) .* gmat,[1 2]),UpWindedFlag, DownWindedFlag, enumNonhydroBoundaryCondition.Zero);
                
                [tempq2x, ~] = obj.NonhydrostaticSolver.matCalculateFluxDownwindedTerm...
                    ( mesh, BoundaryEdge, InnerEdge, num2cell(tempqx,[1 2]),UpWindedFlag, DownWindedFlag, enumNonhydroBoundaryCondition.ZeroGrad);
                [~, tempq2y] = obj.NonhydrostaticSolver.matCalculateFluxDownwindedTerm...
                    ( mesh, BoundaryEdge, InnerEdge, num2cell(tempqy,[1 2]),UpWindedFlag, DownWindedFlag, enumNonhydroBoundaryCondition.ZeroGrad);
                
                [Exactq2x, ~] = obj.NonhydrostaticSolver.matCalculateFluxDownwindedTerm...
                    ( mesh, BoundaryEdge, InnerEdge, num2cell(Exactqx,[1 2]),UpWindedFlag, DownWindedFlag, enumNonhydroBoundaryCondition.ZeroGrad);
                [~, Exactq2y] = obj.NonhydrostaticSolver.matCalculateFluxDownwindedTerm...
                    ( mesh, BoundaryEdge, InnerEdge, num2cell(Exactqy,[1 2]),UpWindedFlag, DownWindedFlag, enumNonhydroBoundaryCondition.ZeroGrad);
                          
                obj.Assert( height(index) .* tempqx, Exactqx );
                obj.Assert( height(index) .* tempqy, Exactqy );
                obj.Assert( height(index) .* tempq2x, Exactq2x );
                obj.Assert( height(index) .* tempq2y, Exactq2y );
            end   
            display('===========================================');
            display('Alternating Upwinded Assemble verified');
            display('===========================================');               
            
            for index = 1:K*Np
                gmat = zeros( Np, K );
                gmat(index) = 1.0/2.0;
                [ tempqx, tempqy ] = obj.NonhydrostaticSolver.GetCharacteristicMatrix...
                    ( mesh, num2cell(gmat,[1 2]), num2cell(gmat,[1 2]), enumNonhydroBoundaryCondition.Zero);
                [ Exactqx, Exactqy ] = obj.NonhydrostaticSolver.GetCharacteristicMatrix...
                    ( mesh, num2cell(height(index) .* gmat,[1 2]), num2cell(height(index) .* gmat,[1 2]), enumNonhydroBoundaryCondition.Zero);
                
                tempq2x = obj.NonhydrostaticSolver.matCalculatePenaltyCharacteristicMatrixX...
                    ( mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]), num2cell(tempqx,[1 2]), obj.NonhydrostaticSolver.EidBoundaryType);
                tempq2y = obj.NonhydrostaticSolver.matCalculatePenaltyCharacteristicMatrixY...
                    ( mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]), num2cell(tempqy,[1 2]), obj.NonhydrostaticSolver.EidBoundaryType);
                
                Exactq2x = obj.NonhydrostaticSolver.matCalculatePenaltyCharacteristicMatrixX...
                    ( mesh, BoundaryEdge, InnerEdge, num2cell(height(index) .* gmat,[1 2]), num2cell(Exactqx,[1 2]), obj.NonhydrostaticSolver.EidBoundaryType);
                Exactq2y = obj.NonhydrostaticSolver.matCalculatePenaltyCharacteristicMatrixY...
                    ( mesh, BoundaryEdge, InnerEdge, num2cell(height(index) .* gmat,[1 2]), num2cell(Exactqy,[1 2]), obj.NonhydrostaticSolver.EidBoundaryType);                
                          
                obj.Assert( height(index) .* tempqx, Exactqx );
                obj.Assert( height(index) .* tempqy, Exactqy );
                obj.Assert( height(index) .* tempq2x, Exactq2x );
                obj.Assert( height(index) .* tempq2y, Exactq2y );
            end   
            display('===========================================');
            display('Stable Central Assemble verified');
            display('===========================================');        
            
        end
        
        function testmatCalculateCharacteristicMatrixX(obj)
            %> @brief function for testing the way to assemble the global stiff matrix
            %> In this study, the global matrix is assembled at the initial
            %> stage of the computation so as to reduce the computational cost.
            %> To do that, a matrix was constructed at first. matCalculateCharacteristicMatrixX
            %> is the function used to calculate the global derivative with the central flux considered in the x
            %> direction
            mesh = obj.meshUnion(1);
            BoundaryEdge = mesh.BoundaryEdge;
            InnerEdge = mesh.BoundaryEdge;
            ExactMatrixX = obj.getExactMatrixX;
            MatrixX = obj.NonhydrostaticSolver.matCalculateCharacteristicMatrixX( mesh, BoundaryEdge, InnerEdge, num2cell(obj.fphys{1}(:,:,1),[1 2]), enumNonhydroBoundaryCondition.Zero);
            obj.Assert( MatrixX, ExactMatrixX );
        end
        
        function testmatCalculateCharacteristicMatrixY(obj)
            %> @brief function for testing the way to assemble the global stiff matrix
            %> In this study, the global matrix is assembled at the initial
            %> stage of the computation so as to reduce the computational cost.
            %> To do that, a matrix was constructed at first. matCalculateCharacteristicMatrixY
            %> is the function used to calculate the global derivative with the central flux considered in the y
            %> direction            
            mesh = obj.meshUnion(1);
            BoundaryEdge = mesh.BoundaryEdge;
            InnerEdge = mesh.BoundaryEdge;
            ExactMatrixY = obj.getExactMatrixY;
            MatrixY = obj.NonhydrostaticSolver.GetCharacteristicMatrixY( mesh, BoundaryEdge, InnerEdge, num2cell(obj.fphys{1}(:,:,1),[1 2]), enumNonhydroBoundaryCondition.Zero);
            obj.Assert( MatrixY, ExactMatrixY );
        end
        
        function testmatCalculateConservativeVariableRelatedMatrixX(obj)
            %> @brief function for testing the calculation of the
            %> conservative variable, [h hu hv], related global partial
            %> derivative with the central flux considered in the x direction
            %>  matCalculateConservativeVariableRelatedMatrixX is the function 
            %> used to calculate the global derivative of the conservative variable
            %> with the central flux considered in the x direction            
           mesh = obj.meshUnion(1);
           BoundaryEdge = mesh.BoundaryEdge;
           InnerEdge = mesh.BoundaryEdge;   
           ExactMatrixX = obj.getExactConservativeVariableRelatedMatrixX;
           MatrixX =  obj.NonhydrostaticSolver.GetConservativeVariableRelatedMatrixX( BoundaryEdge, InnerEdge, obj.fphys, enumNonhydroBoundaryCondition.Zero, 1);
           obj.Assert( MatrixX, ExactMatrixX );
        end
        
        function testmatCalculateConservativeVariableRelatedMatrixY(obj)
            %> @brief function for testing the calculation of the
            %> conservative variable, [h hu hv], related global partial
            %> derivative with the central flux considered in the y direction
            %>  matCalculateConservativeVariableRelatedMatrixY is the function 
            %> used to calculate the global derivative of the conservative variable
            %> with the central flux considered in the y direction                 
           mesh = obj.meshUnion(1);
           BoundaryEdge = mesh.BoundaryEdge;
           InnerEdge = mesh.BoundaryEdge;   
           ExactMatrixY = obj.getExactConservativeVariableRelatedMatrixY;
           MatrixY =  obj.NonhydrostaticSolver.GetConservativeVariableRelatedMatrixY( BoundaryEdge, InnerEdge, obj.fphys, enumNonhydroBoundaryCondition.Zero, 1);
           obj.Assert( MatrixY, ExactMatrixY );
        end 
        
        function testGlobalMatrix(obj)
            obj.NonhydrostaticSolver.dt = 1;
            StiffMatrix = obj.NonhydrostaticSolver.GetGlobalStiffMatrix(obj.NonhydrostaticSolver.PNPX, obj.NonhydrostaticSolver.PNPY, obj.NonhydrostaticSolver.SPNPX,...
    obj.NonhydrostaticSolver.SPNPY, obj.NonhydrostaticSolver.FNPBX, obj.NonhydrostaticSolver.FNPBY, obj.fphys, obj);
            ExactStiffMatrix = obj.getGlobalStiffMatrix;
            obj.Assert( full(StiffMatrix), ExactStiffMatrix );
        end
        
    end
    methods(Access = protected)
        
        function [ option ] = setOption( obj, option )
            ftime = 20;
            outputIntervalNum = 1000;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            %             option('temporalDiscreteType') = enumOutputInterval.DeltaTime;
            %             option('obcType') = NdgBCType.None;
            option('outputCaseName') = mfilename;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputNetcdfCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            %             option('CoriolisType') = SWECoriolisType.None;
            %             option('WindType') = SWEWindType.None;
            %             option('FrictionType') = SWEFrictionType.None;
            option('nonhydrostaticType') = enumNonhydrostaticType.Nonhydrostatic;
        end
        
        function Assert(obj, Exact, Numerical )
            Ind1 = size(Exact,1);Ind2 = size(Exact,2);Ind3 = size(Exact,3);
            for i = 1:Ind1
                for j = 1:Ind2
                    for k = 1:Ind3
                        assert( abs(Exact(i,j,k)- Numerical(i,j,k)) <= obj.tol );
                    end
                end
            end
        end
        
        function [ mesh ] = makeUniformMesh(obj, N, type)
            bctype = [...
                NdgEdgeType.SlipWall, ...
                NdgEdgeType.SlipWall, ...
                NdgEdgeType.SlipWall, ...
                NdgEdgeType.SlipWall];
            
            if (type == NdgCellType.Tri)
                mesh = makeUniformTriMesh(N, [-2, 2], [-1, 1], 2, 1, bctype);
            elseif(type == NdgCellType.Quad)
                mesh = makeUniformQuadMesh(N, [-2, 2], [-1, 1], 2, 1, bctype);
            else
                msgID = [mfile, ':inputCellTypeError'];
                msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
                ME = MException(msgID, msgtext);
                throw(ME);
            end
        end% func
    end
    
end

function [UpWindedFlag, DownWindedFlag] = AssembleWindedFlagInformation(InnerEdge, ele)
DownWindedFlag = zeros(size(InnerEdge.nx));
UpWindedFlag = ones(size(InnerEdge.nx));
[Row, Col] = find(InnerEdge.FToE == ele);
for i = 1:numel(Row)
    if Row(i) == 2
        UpWindedFlag(:,Col(i)) = 0;
        DownWindedFlag(:, Col(i)) = 1;
    end
end
end


