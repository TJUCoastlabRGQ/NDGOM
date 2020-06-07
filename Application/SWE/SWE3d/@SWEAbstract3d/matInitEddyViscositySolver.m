function  matInitEddyViscositySolver( obj )
%> @brief Function to initialize the vertical eddy viscosity sover for both the three dimensional barotropic and baroclinic shallow water model
%> @details Function to initialize the vertical eddy viscosity sover for both the three dimensional barotropic and baroclinic shallow water model
%> @param[in] mesh2d the two-dimensional mesh object
%> @param[in] mesh3d the three-dimensional mesh object
%> @param[out] EddyViscositySolver the eddy viscosity solver to be used in the three-dimensional shallow water model

    if ~obj.option.isKey('VerticalEddyViscosityType')
        msg = 'Eddy viscosity type must be set in the set option part, or the program will be terminated';
        error(msg);
    elseif enumVerticalEddyViscosity.Constant == obj.getOption('VerticalEddyViscosityType')
        obj.VerticalEddyViscositySolver = NdgVertConstantDiffSolver(obj);
%     elseif enumEddyViscosity.Parabolic == obj.getOption('EddyViscosityType')  
%         EddyViscositySolver = NdgParabolicEddyViscositySolver(obj);
    elseif enumVerticalEddyViscosity.GOTM == obj.getOption('VerticalEddyViscosityType')  
        obj.VerticalEddyViscositySolver = NdgVertGOTMDiffSolver(obj);    
    else
        msg = 'None recoganized eddy viscosity type, please confirm again';
        error(msg);
    end
    
    if obj.option.isKey('HorizontalEddyViscosityType')
        if enumHorizontalEddyViscosity.Constant == obj.getOption('HorizontalEddyViscosityType')
            obj.HorizontalEddyViscositySolver = NdgSWEHorizConstantDiffSolver(obj);
        elseif enumHorizontalEddyViscosity.Smagorinsky == obj.getOption('HorizontalEddyViscosityType')
            obj.HorizontalEddyViscositySolver = NdgSWEHorizSmagrinskyDiffSolver(obj);
        elseif enumHorizontalEddyViscosity.None == obj.getOption('HorizontalEddyViscosityType') 
            obj.HorizontalEddyViscositySolver = NdgSWEHorizNoneDiffSolver(obj);
        end
    else
        obj.HorizontalEddyViscositySolver = NdgSWEHorizNoneDiffSolver(obj);
    end

end