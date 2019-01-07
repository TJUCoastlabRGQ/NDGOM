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
        
        function testmatEvaluatePenaltyTerm(obj)
            [ExactPenaltyX, ExactPenaltyY] = obj.getExactPenaltyTerm;
            [fm, fp] = obj.NonhydrostaticSolver.matEvaluateNonhydrostaticSurfaceValue(obj.meshUnion(1), ones(size(obj.meshUnion(1).x)));
            PenaltyX = obj.NonhydrostaticSolver.matEvaluatePenaltyTerm(obj.meshUnion(1), fm, fp, obj.meshUnion(1).nx);
            PenaltyY = obj.NonhydrostaticSolver.matEvaluatePenaltyTerm(obj.meshUnion(1), fm, fp, obj.meshUnion(1).ny);
            obj.Assert(ExactPenaltyX, PenaltyX);
            obj.Assert(ExactPenaltyY, PenaltyY);
        end
        
        function testmatEvaluateVelocityLikeTerms(obj)
            ExactharmonicTerm = obj.getExactHarmonicTerm;
            [fm, fp] = obj.NonhydrostaticSolver.matEvaluateNonhydrostaticSurfaceValue(obj.meshUnion(1), ones(size(obj.meshUnion(1).x)));
            harmonicTerm = obj.NonhydrostaticSolver.matEvaluateVelocityLikeTerms( obj.meshUnion(1), fm, fp );
            obj.Assert(ExactharmonicTerm, harmonicTerm);
        end
        
        function testmatEvaluateNonconservativeFlux(obj)
            [ ExactQx, ExactQy ]= obj.getExactNonconservativeFlux;
            [fm, fp] = obj.NonhydrostaticSolver.matEvaluateNonhydrostaticSurfaceValue(obj.meshUnion(1), ones(size(obj.meshUnion(1).x)));
            harmonicH = obj.NonhydrostaticSolver.matEvaluateVelocityLikeTerms( obj.meshUnion(1), fm, fp );
            PenaltyX = obj.NonhydrostaticSolver.matEvaluatePenaltyTerm(obj.meshUnion(1), fm, fp, obj.meshUnion(1).nx);
            PenaltyY = obj.NonhydrostaticSolver.matEvaluatePenaltyTerm(obj.meshUnion(1), fm, fp, obj.meshUnion(1).ny);
            qx = obj.NonhydrostaticSolver.matEvaluateNonconservativeFlux( harmonicH, fm, fp, PenaltyX);
            qy = obj.NonhydrostaticSolver.matEvaluateNonconservativeFlux( harmonicH, fm, fp, PenaltyY);
            obj.Assert(ExactQx, qx);
            obj.Assert(ExactQy, qy);
        end
        
        function testmatEvaluateDeltaSurfaceFlux(obj)
            ExactDeltaFlux = obj.getExactDeltaFlux;
            [fm, fp] = obj.NonhydrostaticSolver.matEvaluateNonhydrostaticSurfaceValue(obj.meshUnion(1), ones(size(obj.meshUnion(1).x)));
            harmonicH = obj.NonhydrostaticSolver.matEvaluateVelocityLikeTerms( obj.meshUnion(1), fm, fp );
            PenaltyX = obj.NonhydrostaticSolver.matEvaluatePenaltyTerm(obj.meshUnion(1), fm, fp, obj.meshUnion(1).nx);
            PenaltyY = obj.NonhydrostaticSolver.matEvaluatePenaltyTerm(obj.meshUnion(1), fm, fp, obj.meshUnion(1).ny);
            qx = obj.NonhydrostaticSolver.matEvaluateNonconservativeFlux( harmonicH, fm, fp, PenaltyX);
            qy = obj.NonhydrostaticSolver.matEvaluateNonconservativeFlux( harmonicH, fm, fp, PenaltyY);
            DeltaFlux = obj.NonhydrostaticSolver.matEvaluateDeltaSurfaceFlux( obj.meshUnion(1), ones(size(obj.meshUnion(1).x)), ones(size(obj.meshUnion(1).x)), qx, qy );
            obj.Assert(ExactDeltaFlux, DeltaFlux);
        end
        
        function testmatEvaluateDeltaSurfaceFluxY(obj)
            ExactDeltaFlux = obj.getExactDeltaFluxY;
            [fm, fp] = obj.NonhydrostaticSolver.matEvaluateNonhydrostaticSurfaceValue(obj.meshUnion(1), ones(size(obj.meshUnion(1).x)));
            PenaltyY = zeros(size(obj.meshUnion(1).eidM));
            qY = obj.NonhydrostaticSolver.matEvaluateNonconservativeFlux( fm, fp, PenaltyY);
            DeltaFlux = obj.NonhydrostaticSolver.matEvaluateDeltaSurfaceFluxY(obj.meshUnion(1), ones(size(obj.meshUnion(1).eidM)), qY);
            obj.Assert(ExactDeltaFlux, DeltaFlux);
        end
        
        function testmatEvaluateDeltaSurfaceFluxX(obj)
            ExactDeltaFlux = obj.getExactDeltaFluxX;
            [fm, fp] = obj.NonhydrostaticSolver.matEvaluateNonhydrostaticSurfaceValue(obj.meshUnion(1), ones(size(obj.meshUnion(1).x)));
            PenaltyX = zeros(size(obj.meshUnion(1).eidM));
            qX = obj.NonhydrostaticSolver.matEvaluateNonconservativeFlux( fm, fp, PenaltyX);
            DeltaFlux = obj.NonhydrostaticSolver.matEvaluateDeltaSurfaceFluxX(obj.meshUnion(1), ones(size(obj.meshUnion(1).eidM)), qX);
            obj.Assert(ExactDeltaFlux, DeltaFlux);
        end
        
        function testmatEvaluateDivergenceDirectionX(obj)
            ExactDivergenceX = obj.getExactDivergenceX;
            DivergenceX = obj.NonhydrostaticSolver.matEvaluateDivergenceDirectionX(obj.meshUnion(1), ones(size(obj.meshUnion(1).x)));
            obj.Assert(ExactDivergenceX, DivergenceX);
        end
        
        function testmatEvaluateDivergenceDirectionY(obj)
            ExactDivergenceY = obj.getExactDivergenceY;
            DivergenceY = obj.NonhydrostaticSolver.matEvaluateDivergenceDirectionY(obj.meshUnion(1), ones(size(obj.meshUnion(1).x)));
            obj.Assert(ExactDivergenceY, DivergenceY);
        end
        
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
        
        %         function testReverseEidBoundaryType(obj)
        %             ExactReverseEidBoundaryType = obj.getExactReverseEidBoundaryType;
        %             obj.Assert(ExactReverseEidBoundaryType, obj.NonhydrostaticSolver.ReverseEidBoundaryType);
        %         end
                
        function testmatGetFaceValue(obj)
            %> @brief function for testing the Nonhydrostatic pressure at face point
            %             ExactGetFaceOuterValue = obj.getExactGetFaceOuterValue;
            mesh = obj.meshUnion(1);
            obj.NonhydrostaticSolver.getWetDryInterface(mesh);
            
            [~, ~, ~, ~, ~, ~]= obj.NonhydrostaticSolver.matReconstructStiffmatrixRelatedMatrix( obj); %to get the obj.NonhydroFmPoint and obj.NonhydroFpPoint
            [fm, fp] = mesh.InnerEdge.matEvaluateSurfValue( obj.fphys );       
            [fm, fp] = obj.NonhydrostaticSolver.GetFaceValue(fm(:,:,1), fp(:,:,1), enumNonhydroBoundaryCondition.Zero); % for this case, the water depth is treated as the non-hydrostatic pressure to test the inner and outer value
            
            [ExactFm, ExactFp] = obj.getExactInnerOuterValue;
            
            obj.Assert(ExactFm, fm);
            obj.Assert(ExactFp, fp);
        end
        
%         function testTopographyNonhydrostaticFaceOuterValue(obj)
%             %> @brief function for testing the Topography Nonhydrostatic at face point
%             ExactGetFaceOuterValue = obj.getExactGetFaceOuterValue;
%             mesh = obj.meshUnion(1);
%             obj.NonhydrostaticSolver.getWetDryInterface( mesh );
%             Inner = obj.zGrad{1}(:,:,1);
%             Outer = Inner(mesh.eidM);
%             Outer = obj.NonhydrostaticSolver.getFaceOuterValue(mesh, Outer, -Inner);
%             obj.Assert(ExactGetFaceOuterValue, Outer);
%         end
        
        function testUpdatedExactEidBoundaryType(obj)
            mesh = obj.meshUnion(1);
            ExactEidBoundaryType = obj.getUpdatedExactEidBoundaryType;
            obj.NonhydrostaticSolver.getWetDryInterface( mesh );
            EidBoundaryType = obj.NonhydrostaticSolver.getEidBoundaryType(mesh);
            obj.Assert(ExactEidBoundaryType, EidBoundaryType);
        end
        
        function testFluxTerm(obj)
            mesh = obj.meshUnion(1);
            obj.NonhydrostaticSolver.getWetDryInterface( mesh );
            ExactFluxterm = obj.GetExactFluxterm(mesh);
%             Fluxterm = ones(size(mesh.eidM));
            Fluxterm =  obj.NonhydrostaticSolver.GetFluxTerm(mesh, Fluxterm);
            obj.Assert(ExactFluxterm, Fluxterm);
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


