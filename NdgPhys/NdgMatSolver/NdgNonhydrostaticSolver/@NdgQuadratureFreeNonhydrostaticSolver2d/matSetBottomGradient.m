function [ bx, by ] = matSetBottomGradient(obj, zGrad)
bx = zGrad(:,:,1);
by = zGrad(:,:,2);
end