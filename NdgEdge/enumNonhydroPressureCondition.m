%> @brief enumeration for edge types.
%
%> And here we can put some more detailed informations about the class.
% ======================================================================
%> This class is part of the NDGOM software.
%> @author li12242, Tianjin University, li12242@tju.edu.cn
% ======================================================================
classdef enumNonhydroPressureCondition < int8
    
    enumeration
        Inner           (0)
        Newmann         (1)
        Dirichlet       (2)
    end
end