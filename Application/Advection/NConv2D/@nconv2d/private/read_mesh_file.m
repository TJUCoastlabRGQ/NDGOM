function [ mesh ] = read_mesh_file( N, type, casename )
%READ_MESH_FILE Summary of this function goes here
%   Detailed explanation goes here

% ��ȡ��ά��������
switch type
    case ndg_lib.std_cell_type.Tri
        cell = ndg_lib.std_cell.tri(N);
        mesh = ndg_lib.mesh.tri_mesh(cell, 'file', casename);
    case ndg_lib.std_cell_type.Quad
        cell = ndg_lib.std_cell.quad(N);
        mesh = ndg_lib.mesh.quad_mesh(cell, 'file', casename);
end% func

end


