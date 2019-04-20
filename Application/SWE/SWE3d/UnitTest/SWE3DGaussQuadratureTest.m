classdef SWE3DGaussQuadratureTest < SWE3DAbstractTest
    %SWE3DGAUSSQUADRATURETEST 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        N = [1 2 3 4]
        Nz = [1 2 3 4]
        M = 10
        Mz = [3 4 5]
    end
    
    methods
        
        function obj = SWE3DGaussQuadratureTest()
            %             obj.Solver = StandingWaveInAClosedChannel()
        end
        
        function test_InterpolationMatrix(obj)
            for  i = 1:numel(obj.N)
                for j = 1:numel(obj.Nz)
                    for k = 1:numel(obj.M)
                        for t = 1:numel(obj.Mz)
                            str = sprintf('N = %d, Nz = %d, M = %d, Mz = %d', obj.N(i), obj.Nz(j), obj.M(k), obj.Mz(t));
                            disp(str);
                            Solver = StandingWaveInAClosedChannel(obj.N(i), obj.Nz(j), obj.M(k), obj.Mz(t));
                            obj.Assert(sum( Solver.advectionSolver.SBFVfq{1}, 2 ), ones( size(Solver.advectionSolver.SBFVfq{1},1) ,1 ));
                            obj.Assert(sum( Solver.advectionSolver.BOTFVfq{1}, 2 ), ones( size(Solver.advectionSolver.BOTFVfq{1},1) ,1 ));
                            obj.Assert(sum( Solver.advectionSolver.BBFVfq{1}, 2 ), ones( size(Solver.advectionSolver.BBFVfq{1},1) ,1 ));
                            obj.Assert(sum( Solver.advectionSolver.IEFVfq{1}, 2 ), ones( size(Solver.advectionSolver.IEFVfq{1},1) ,1 ));
                            obj.Assert(sum( Solver.advectionSolver.BEFVfq{1}, 2 ), ones( size(Solver.advectionSolver.BEFVfq{1},1) ,1 ));
                        end
                    end
                end
            end
        end
        
        function test_GaussPointDirectionVector(obj)
            for  i = 1:numel(obj.N)
                for j = 1:numel(obj.Nz)
                    for k = 1:numel(obj.M)
                        for t = 1:numel(obj.Mz)
                            str = sprintf('N = %d, Nz = %d, M = %d, Mz = %d', obj.N(i), obj.Nz(j), obj.M(k), obj.Mz(t));
                            disp(str);
                            Solver = StandingWaveInAClosedChannel(obj.N(i), obj.Nz(j), obj.M(k), obj.Mz(t));
                            %> This part is used to verify the direction vector at the Gauss point is equal to that at the interpolation point
                            obj.Assert(Solver.advectionSolver.SBnx{1}(1,:), Solver.meshUnion(1).SurfaceBoundaryEdge.nx(1,:));
                            obj.Assert(Solver.advectionSolver.SBny{1}(1,:), Solver.meshUnion(1).SurfaceBoundaryEdge.ny(1,:));
                            obj.Assert(Solver.advectionSolver.SBnz{1}(1,:), Solver.meshUnion(1).SurfaceBoundaryEdge.nz(1,:));
                            obj.Assert(Solver.advectionSolver.BOTnx{1}(1,:),Solver.meshUnion(1).BottomEdge.nx(1,:));
                            obj.Assert(Solver.advectionSolver.BOTny{1}(1,:),Solver.meshUnion(1).BottomEdge.ny(1,:));
                            obj.Assert(Solver.advectionSolver.BOTnz{1}(1,:),Solver.meshUnion(1).BottomEdge.nz(1,:));
                            obj.Assert(Solver.advectionSolver.BBnx{1}(1,:), Solver.meshUnion(1).BottomBoundaryEdge.nx(1,:));
                            obj.Assert(Solver.advectionSolver.BBny{1}(1,:), Solver.meshUnion(1).BottomBoundaryEdge.ny(1,:));
                            obj.Assert(Solver.advectionSolver.BBnz{1}(1,:), Solver.meshUnion(1).BottomBoundaryEdge.nz(1,:));
                            obj.Assert(Solver.advectionSolver.IEnx{1}(1,:), Solver.meshUnion(1).InnerEdge.nx(1,:));
                            obj.Assert(Solver.advectionSolver.IEny{1}(1,:), Solver.meshUnion(1).InnerEdge.ny(1,:));
                            obj.Assert(Solver.advectionSolver.IEnz{1}(1,:), Solver.meshUnion(1).InnerEdge.nz(1,:));
                            obj.Assert(Solver.advectionSolver.BEnx{1}(1,:), Solver.meshUnion(1).BoundaryEdge.nx(1,:));
                            obj.Assert(Solver.advectionSolver.BEny{1}(1,:), Solver.meshUnion(1).BoundaryEdge.ny(1,:));
                            obj.Assert(Solver.advectionSolver.BEnz{1}(1,:), Solver.meshUnion(1).BoundaryEdge.nz(1,:));
                            %> This part is used to verify the direction vector at the same face is the same
                            obj.Assert(Solver.advectionSolver.SBnx{1}, repmat(Solver.advectionSolver.SBnx{1}(1,:),size(Solver.advectionSolver.SBnx{1},1)));
                            obj.Assert(Solver.advectionSolver.SBny{1}, repmat(Solver.advectionSolver.SBny{1}(1,:),size(Solver.advectionSolver.SBny{1},1)));
                            obj.Assert(Solver.advectionSolver.SBnz{1}, repmat(Solver.advectionSolver.SBnz{1}(1,:),size(Solver.advectionSolver.SBnz{1},1)));
                            obj.Assert(Solver.advectionSolver.BOTnx{1}, repmat(Solver.advectionSolver.BOTnx{1}(1,:),size(Solver.advectionSolver.BOTnx{1},1)));
                            obj.Assert(Solver.advectionSolver.BOTny{1}, repmat(Solver.advectionSolver.BOTny{1}(1,:),size(Solver.advectionSolver.BOTny{1},1)));
                            obj.Assert(Solver.advectionSolver.BOTnz{1}, repmat(Solver.advectionSolver.BOTnz{1}(1,:),size(Solver.advectionSolver.BOTnz{1},1)));
                            obj.Assert(Solver.advectionSolver.BBnx{1}, repmat(Solver.advectionSolver.BBnx{1}(1,:),size(Solver.advectionSolver.BBnx{1},1)));
                            obj.Assert(Solver.advectionSolver.BBny{1}, repmat(Solver.advectionSolver.BBny{1}(1,:),size(Solver.advectionSolver.BBny{1},1)));
                            obj.Assert(Solver.advectionSolver.BBnz{1}, repmat(Solver.advectionSolver.BBnz{1}(1,:),size(Solver.advectionSolver.BBnz{1},1)));
                            obj.Assert(Solver.advectionSolver.IEnx{1}, repmat(Solver.advectionSolver.IEnx{1}(1,:),size(Solver.advectionSolver.IEnx{1},1)));
                            obj.Assert(Solver.advectionSolver.IEny{1}, repmat(Solver.advectionSolver.IEny{1}(1,:),size(Solver.advectionSolver.IEny{1},1)));
                            obj.Assert(Solver.advectionSolver.IEnz{1}, repmat(Solver.advectionSolver.IEnz{1}(1,:),size(Solver.advectionSolver.IEnz{1},1)));
                            obj.Assert(Solver.advectionSolver.BEnx{1}, repmat(Solver.advectionSolver.BEnx{1}(1,:),size(Solver.advectionSolver.BEnx{1},1)));
                            obj.Assert(Solver.advectionSolver.BEny{1}, repmat(Solver.advectionSolver.BEny{1}(1,:),size(Solver.advectionSolver.BEny{1},1)));
                            obj.Assert(Solver.advectionSolver.BEnz{1}, repmat(Solver.advectionSolver.BEnz{1}(1,:),size(Solver.advectionSolver.BEnz{1},1)));
                        end
                    end
                end
            end
        end
        
        function test_Jacobian(obj)
            for  i = 1:numel(obj.N)
                for j = 1:numel(obj.Nz)
                    for k = 1:numel(obj.M)
                        for t = 1:numel(obj.Mz)
                            str = sprintf('N = %d, Nz = %d, M = %d, Mz = %d', obj.N(i), obj.Nz(j), obj.M(k), obj.Mz(t));
                            disp(str);
                            Solver = StandingWaveInAClosedChannel(obj.N(i), obj.Nz(j), obj.M(k), obj.Mz(t));
                            switch Solver.mesh3d.cell.type
                                case enumStdCell.PrismTri
                                    % for the standard tri, sum of the  Jacobian weight is two
                                    obj.Assert(sum(Solver.advectionSolver.SBwJs{1})./Solver.meshUnion(1).SurfaceBoundaryEdge.Js(1,:), 2*ones( 1,Solver.meshUnion(1).SurfaceBoundaryEdge.Ne));
                                    obj.Assert(sum(Solver.advectionSolver.BOTwJs{1})./Solver.meshUnion(1).BottomEdge.Js(1,:), 2*ones( 1,Solver.meshUnion(1).BottomEdge.Ne));
                                    obj.Assert(sum(Solver.advectionSolver.BBwJs{1})./Solver.meshUnion(1).BottomBoundaryEdge.Js(1,:), 2*ones( 1,Solver.meshUnion(1).BottomBoundaryEdge.Ne));
                                otherwise
                                    %doing nothing, actually we need to consider the PrismQuad condition, but value Js is not a constant for this condition
                                    %if the mesh is not uniform, so we ommit it here
                            end
                            % for the standard quad, sum of the Jacobian weight is two
                            obj.Assert(sum(Solver.advectionSolver.IEwJs{1})./Solver.meshUnion(1).InnerEdge.Js(1,:), 4*ones( 1,Solver.meshUnion(1).InnerEdge.Ne));
                            obj.Assert(sum(Solver.advectionSolver.BEwJs{1})./Solver.meshUnion(1).BoundaryEdge.Js(1,:), 4*ones( 1,Solver.meshUnion(1).BoundaryEdge.Ne));
                        end
                    end
                end
            end
        end
    end
    
end

