function EddyViscositySolver = matInitEddyViscositySolver( obj )
%> @brief Function to initialize the vertical eddy viscosity sover for both the three dimensional barotropic and baroclinic shallow water model
%> @details Function to initialize the vertical eddy viscosity sover for both the three dimensional barotropic and baroclinic shallow water model
%> @param[in] mesh2d the two-dimensional mesh object
%> @param[in] mesh3d the three-dimensional mesh object
%> @param[out] EddyViscositySolver the eddy viscosity solver to be used in the three-dimensional shallow water model

    if ~obj.option.isKey('EddyViscosityType')
        msg = 'Eddy viscosity type must be set in the set option part, or the program will be terminated';
        error(msg);
    elseif enumEddyViscosity.Constant == obj.getOption('EddyViscosityType')
        EddyViscositySolver = NdgConstantEddyViscositySolver(obj);
    elseif enumEddyViscosity.Parabolic == obj.getOption('EddyViscosityType')  
        EddyViscositySolver = NdgParabolicEddyViscositySolver(obj);
    elseif enumEddyViscosity.GOTM == obj.getOption('EddyViscosityType')  
        EddyViscositySolver = NdgGOTMEddyViscositySolver(obj);    
    else
        msg = 'None recoganized eddy viscosity type, please confirm again';
        error(msg);
    end

end