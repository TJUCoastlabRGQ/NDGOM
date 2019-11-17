classdef NdgAbstractEddyViscositySolver < handle
    %NDGABSTRACTEDDYVISCOSITYSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods( Abstract, Access = public )
        %> @brief Function to update the vertical eddy viscosity during computation
        %> @details Function to update the vertical eddy viscosity during computation
        %> @param[in] mesh3d the three-dimensional mesh object, the VTV interpolation  matrix containd in its standard cell is used to calculate its central value from the value at interpolation points
        %> @param[in] fphys2d the 2d physical field which contains the water depth field, this is used to calculate the height of each cell 
        %> @param[in] fphys the three dimensional physical field, this is used to calculate the shear production and buoyance production related terms
        %> @param[in] dt the time step of the three-dimensional shallow water model, this variable is required for both the TKE equation and the turbulent length scale related equation if GOTM adopted
        %> @param[in] time the time point at the current step, this is used to check whether this is the final point of the simulation, if so, all the resources allocated at the initial stage would be deallocated.
        matUpdateEddyViscosity( obj, mesh3d, fphys2d, fphys, dt, time );
        
    end
    
    methods( Abstract, Access = protected )
        %> @brief Function to initialize the vertical eddy viscosity
        %> @details Function to update the vertical eddy viscosity
        %> @param[in] physClass the three-dimensional shallow water solver, the viscosity will be set here
        %> @param[in] mesh2d the two-dimensional mesh object, and the cell number and interpolation point number contained in this object is used to initialize the eddy viscosity related matrix during GOTM computation
        %> @param[in] mesh3d the three-dimensional mesh object, and the cell number and interpolation point number contained in this object is used to initialize the eddy viscosity related matrix during GOTM computation
        %> @param[in] hcrit the critical water depth, if the depth is less than this value, the velocity is set to be zero.        
        matInitEddyViscosity(obj, physClass, mesh2d, mesh3d, hcrit );
        
    end
    
end

