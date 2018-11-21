function [ fext ] = getInitialFunction( obj )
fext = cell( obj.Nmesh, 1 );

for m = 1:obj.Nmesh%����ѭ��
    fext{m} = zeros(obj.meshUnion(m).cell.Np, obj.meshUnion(m).K, obj.Nfield);%��СNp*K*Nfield
    a = obj.meshUnion(m).y;
    tx = -0.2*cos(3.1415927*a/1000000);
    
    fext{m}(:,:,1) = 1000;
    fext{m}(:,:,7) = tx;
end

end%func

