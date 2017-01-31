function eta = getIntrinsicImpedance(this, frequency)
% eta = getIntrinsicImpedance(this, frequency)
%
% Input:
% 
% frequency     Nx1 vector (Hz)
    
% $Author:: kzhu                                       $
% $Rev:: 1487                                          $
% $Date:: 2011-02-11 01:56:43 -0500 (Fri, 11 Feb 2011) $
    

    EPS_O = 8.8541878176e-12;
    MU_O  = 4*pi*1e-7;
    ETA_O = sqrt(MU_O/EPS_O);
    eta   = ETA_O*sqrt(getComplexPermeability(this, frequency)./ ...
                       getComplexPermittivity(this, frequency));
    
    
    