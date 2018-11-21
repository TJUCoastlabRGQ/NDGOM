function EvaluaetNodeCoor( obj, Nh, Nz )
    
    [ Nph, r1, s1 ] = quad_node_coor_func( Nh );
    [ Npz, t1 ] = line_node_coor_func( Nz );

    Np = Nph * Npz;
    r = repmat(r1, 1, Npz);
    s = repmat(s1, 1, Npz);
    t = repmat(t1', Nph, 1 );

    %r = r(:); s = s(:); t = t(:);
    
    obj.Np = Np; 
    obj.Nph = Nph; 
    obj.Npz = Npz;
    
    obj.r1 = r1; 
    obj.s1 = s1; 
    obj.t1 = t1;
    
    obj.r = r(:);
    obj.s = s(:);
    obj.t = t(:);
end

function [ Np, t ] = line_node_coor_func( N )
    Np = N+1;
    [ t,~ ] = zwglj(Np);
end

function [ Np, r, s ] = quad_node_coor_func(nOrder)
   np = nOrder+1;
   [x,~] = zwglj(np);

   r = x * ones(1, np);
   s = ones(np, 1) * x';

   r = r(:); 
   s = s(:);
   Np = np*np;

end
    