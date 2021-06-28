function runMeshMethodsTest
%NOTE: This function is used to verify the matlab version and the C version
% are exactly the same. Generally, the matlab version in the corresponding
% call is not activated. When call this function, activate that part first,
% or the script would hault!!
N = [1 2];
Nz = [1 2];
Layer = [6 7 8 9 10];
for i = 1:numel(N)
    for j = 1:numel(Nz)
        for k = 1:numel(Layer)
            Solver = WindDrivenFlow(N(i),Nz(j),100,Layer(k));
            Field3d = rand(Solver.meshUnion.cell.Np, Solver.meshUnion.K);
            [ ~ ] = Solver.meshUnion.VerticalIntegralField( Field3d );
            Field2d = rand(Solver.meshUnion.cell.N + 1, Solver.meshUnion.mesh2d.InnerEdge.Ne);
            [ ~ ] = Solver.meshUnion.InnerEdge.matVerticalRepmatFacialValue(Field2d);
            
            Field2d = rand(Solver.meshUnion.cell.N + 1, Solver.meshUnion.mesh2d.BoundaryEdge.Ne);
            [ ~ ] = Solver.meshUnion.BoundaryEdge.matVerticalRepmatFacialValue(Field2d);
            
            Field3d = rand(Solver.meshUnion.InnerEdge.Nfp, Solver.meshUnion.InnerEdge.Ne);
            [ ~ ] = Solver.meshUnion.InnerEdge.VerticalColumnIntegralField( Field3d );
           
            Field3d = rand(Solver.meshUnion.BoundaryEdge.Nfp, Solver.meshUnion.BoundaryEdge.Ne);
            [ ~ ] = Solver.meshUnion.BoundaryEdge.VerticalColumnIntegralField( Field3d );

            clear Solver
        end
    end
end



end

