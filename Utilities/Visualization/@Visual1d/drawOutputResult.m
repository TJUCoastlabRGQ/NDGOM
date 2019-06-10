function drawOutputResult( obj, outputObj, step, fieldId )
    field = outputObj.readOutputResult(step);
    
    isFigureEmpty = isempty( obj.drawHandle );
    if isFigureEmpty
        drawNewFigure( obj, field(:, :, fieldId) );
        return;
    end

    if isvalid( obj.drawHandle )
        updateFigure( obj, field(:, :, fieldId) );
    else
        drawNewFigure( obj, field(:, :, fieldId) );
    end
end

function updateFigure( obj, field )
    set(obj.drawHandle, 'XData', obj.mesh.x(:),...
        'YData', field(:));
%     trisurf( obj.tri, obj.mesh.x(:), obj.mesh.y(:), field(:) )
end

function drawNewFigure(obj, field)
    obj.drawHandle = plot( obj.mesh.x(:), obj.mesh.y(:), field(:) );
    hold on;
    box on;
    grid on;
end