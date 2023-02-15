function [Rarray, Carray, Rterminating, Cterminating] = makeCPE(varargin)
% [Rarray, Carray, Rterminating, Cterminating] = MakeCPEV(alpha, C_frac, fmin, fmax, k_f)
% or
% [Rarray, Carray, Rterminating, Cterminating] = MakeCPEV(alpha, Z0, f0, fmin, fmax, k_f)
% Author: Vance Farrow, Waikato University 2023
% With assistance from Marcus Wilson, Jonathan Scott, Logan Cowie


if(length(varargin) == 5)
    alpha = varargin{:, 1};
    C_frac = varargin{:, 2};
    fmin = varargin{:, 3};
    fmax = varargin{:, 4};
    k_f = varargin{:, 5};
    
    f0 = sqrt(fmin*fmax);
    Z0 = 1./(C_frac *(2*pi*f0)^alpha);
    
elseif(length(varargin) == 6)
    alpha = varargin{:, 1};
    Z0 = varargin{:, 2};
    f0 = varargin{:, 3};
    fmin = varargin{:, 4};
    fmax = varargin{:, 5};
    k_f = varargin{:, 6};

%     C_frac = 1./(Z0 *(2*pi*f0)^alpha);

else
    error("%s\n%s\n%s\n%s", "Expected usage:", ...
        "makeCPE(alpha, C_frac, fmin, fmax, k_f)", ...
        "or", "makeCPE(alpha, Z0, f0, fmin, fmax, k_f)");
end

m = 1./alpha;
k = k_f^alpha;

N_h = floor(log(fmax/f0)/log(k_f));
N_l = floor(log(f0/fmin)/log(k_f));

y_theta = (pi/(m*log(k))) * sec(pi/2 * (1-2/m));
R0 = Z0*y_theta;
C0 = 1./(2*pi*f0*R0);

Rarray = R0.*k.^(-N_h:1:N_l);
Carray = C0*(k.^(m-1)).^(-N_h:1:N_l);

Rterminating = Rarray(end)*(k-1);
Cterminating = Carray(1)/(k^(m-1)-1);

end


