function [ bx, by ] = matSetBottomGradient(obj, zGrad)
bx = obj.Vq{1} * zGrad(:,:,1);
by = obj.Vq{1} * zGrad(:,:,2);
end