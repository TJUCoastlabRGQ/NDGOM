function [fm, fp] = matEvaluateNonhydrostaticSurfaceValue( obj, mesh, Variable )
%<@brief function used to get the inner and outer face value
%<@detail for points are not inner type, the value is the opposite number of the inner value
%<@param[in] mesh The mesh object 
%<@param[in] Variable The variable need to got the both the inner and outer face value
%<@param[out] fm The inner face value
%<@param[out] fp The outer face value
index = (mesh.eidtype ~= 0);
fm = Variable(mesh.eidM);
fp = Variable(mesh.eidP);
fp(index) = -fm(index);
end