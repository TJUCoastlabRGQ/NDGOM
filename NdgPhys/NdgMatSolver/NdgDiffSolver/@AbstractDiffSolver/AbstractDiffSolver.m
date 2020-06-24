classdef AbstractDiffSolver < handle
    %ABSTRACTDIFFSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Access = protected)
        %> the diffusion coefficient
        nv
        %> the penalty parameter for IP form
        tau
        %> the Prantal number for heat and sault diffusion
        Prantl
    end
    
    methods
        function   obj = AbstractDiffSolver( physClass )
            obj.nv = zeros(size(physClass.meshUnion(1).x));
        end
    end
    
    
    methods( Access = protected, Abstract)
        %> this function is used to update the viscosity
        matUpdateViscosity(obj)
        %> this function is used to update the penalty parameter adopted in
        %> IP form
        matUpdatePenaltyParameter(obj)
        
    end
    
    methods(Access = protected)
        function [ fm, fp ] = matEvaluateSurfValue(obj, edge, Kappa )
            [ fm, fp ] = mxEvaluateSurfValue( edge.FToE, edge.FToN1, edge.FToN2, Kappa );
        end
        
    end
    
end

