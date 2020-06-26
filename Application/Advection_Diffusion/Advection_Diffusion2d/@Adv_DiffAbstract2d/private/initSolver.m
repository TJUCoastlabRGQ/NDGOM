function [ Adv, HorizVis, VertVis ] = initSolver( obj )
if obj.option.isKey('AdvDiffHorizontalDiffusionType')
    type = obj.getOption('AdvDiffHorizontalDiffusionType');
    if type == enumHorizontalDiffusion.Constant
        HorizVis = NdgHorizDiffSolver( obj );
    else
        HorizVis = NdgNoneHorizDiffSolver( obj );
    end
else
    HorizVis = NdgNoneHorizDiffSolver( obj );
end

VertVis = NdgNoneVertDiffSolver( obj );


integralType = obj.getOption('integralType');
equType = obj.getOption('equationType');

if (integralType == enumDiscreteIntegral.QuadratureFree)
    if( equType == enumDiscreteEquation.Strong )
        Adv = NdgQuadFreeStrongFormAdvSolver2d( obj );
    elseif( equType == enumDiscreteEquation.Weak )
        Adv = NdgQuadFreeWeakFormAdvSolver2d( obj );
    end
elseif (integralType == enumDiscreteIntegral.GaussQuadrature)
    if( equType == enumDiscreteEquation.Strong )
        Adv = NdgGaussQuadStrongFormAdvSolver2d( obj );
    elseif( equType == enumDiscreteEquation.Weak )
        Adv = NdgGaussQuadWeakFormAdvSolver2d( obj );
    end
end
end