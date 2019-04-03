function [ advectionSolver, viscositySolver ] = initSolver2d( obj )
%INITSOLVER2D 此处显示有关此函数的摘要
%   此处显示详细说明


integralType = obj.getOption('integralType');
equType = obj.getOption('equationType');
advectionSolver = initAdvSolver2d( obj.Solver2d, integralType, equType );
viscositySolver = NdgNonVisSolver( obj.Solver2d );

end


function [ adv ] = initAdvSolver2d( physMat, integralType, equType )

if (integralType == enumDiscreteIntegral.QuadratureFree)
    if( equType == enumDiscreteEquation.Strong )
        adv = NdgQuadFreeStrongFormAdvSolver2d( physMat );
    elseif( equType == enumDiscreteEquation.Weak )
        adv = NdgQuadFreeWeakFormAdvSolver2d( physMat );
    end
elseif (integralType == enumDiscreteIntegral.GaussQuadrature)
    if( equType == enumDiscreteEquation.Strong )
        adv = NdgGaussQuadStrongFormAdvSolver2d( physMat );
    elseif( equType == enumDiscreteEquation.Weak )
        adv = NdgGaussQuadWeakFormAdvSolver2d( physMat );
    end
end

end