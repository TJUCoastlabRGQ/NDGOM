function drawMesh(obj)
    mesh = obj.mesh;
    plot(mesh.x, ...
        zeros(size(mesh.x)),'Marker','o');
    box on;
    grid on;
end