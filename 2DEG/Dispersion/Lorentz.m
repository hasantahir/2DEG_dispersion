function epsilon = Lorentz(ep_0, ep_inf, w_t, g_0,w)
% Returns the Lorentz model based effective permittivity
% Input parameters
% ep_0 dc limit of the permittivity
% ep_inf high-frequency limit of permittivity
% w_t Phonon frequency/ Plasma Frequency
% g_0 Relaxation rate
% w angular frequency
epsilon = ep_inf + (ep_0 - ep_inf)*w_t^2./(w_t^2 - w.^2 - 1i*g_0*w);