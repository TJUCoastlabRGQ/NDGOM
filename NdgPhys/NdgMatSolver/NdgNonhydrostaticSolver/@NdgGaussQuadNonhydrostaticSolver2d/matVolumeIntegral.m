function [ VolumeIntegral ] = matVolumeIntegral( obj, Variable )
VolumeIntegral = zeros(size(obj.Vq{1},2), size(Variable,2));
for i = 1:size(Variable,2)
    VolumeIntegral(:,i) = obj.Vq' * diag(obj.Jq{1}(:,i)) * obj.Vq * Variable;
end

end