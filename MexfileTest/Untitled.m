Solver = EllipticMixedParticalDerivativeTest2d(1,20);
mesh = Solver.meshUnion;
Fmask = mesh.cell.Fmask;
for i = 1:mesh.InnerEdge.Ne
    f = mesh.InnerEdge.FToF(2*(i-1)+1);
    data = mesh.InnerEdge.FToN1(:,i) - Fmask(:,f);
    if(abs(any(data))> 0.1)
        disp("Wrong data at local side of inner edge, and the order is:")
        disp(i);
    end
    f = mesh.InnerEdge.FToF(2*(i-1)+2);
    data = mesh.InnerEdge.FToN2(:,i) - Fmask(:,f);
    if(abs(any(data))> 0.1)
        disp("Wrong data at adjacent side of inner edge, and the order is:")
        disp(i);
    end    
end