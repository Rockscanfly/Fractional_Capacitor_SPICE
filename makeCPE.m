function [Rarray, Carray, Rterminating, Cterminating] = makeCPE(filename, varargin)
% [Rarray, Carray, Rterminating, Cterminating] = MakeCPEV(filename, alpha, C_frac, fmin, fmax, k_f)
% or
% [Rarray, Carray, Rterminating, Cterminating] = MakeCPEV(filename, alpha, Z0, f0, fmin, fmax, k_f)
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
        "makeCPE(filename, alpha, C_frac, fmin, fmax, k_f)", ...
        "or", "makeCPE(filename, alpha, Z0, f0, fmin, fmax, k_f)");
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


k_f_perDecade=(k_f-1)* log10(fmax/fmin); %the highest accuracy k_f value usually gives a value of around 6.1
suggested_k_f=1+(7/log10(fmax/fmin));


if(~isempty(filename))
    fileID = fopen(filename, 'w+');
    if(fileID)
        fprintf(fileID,"* CPE alpha = %3.6e\n", alpha);
		fprintf(fileID,"* C_frac = %3.6e\n", C_frac);
		f0 = sqrt(fmin*fmax);	
        fprintf(fileID,"* CPE |Z|   = %3.6e Ohms\n", Z0);
		fprintf(fileID,"*  at f0    = %3.6eHz or %3.6e rad/%3.6e\n", f0, f0.*2.*pi);
        fprintf(fileID,"* CPE fmin/fmax = %03.4e/%03.4e\n", fmin, fmax);
        fprintf(fileID,"* CPE k = %.6f\n", k);
        fprintf(fileID,"* CPE k_f = %03.6e, recommended k_f value = %03.6e\n", k_f, suggested_k_f);

        fprintf(fileID,".subckt %s cpe1 cpe2\n", filename);


        for idx = 1:length(Rarray)
            thisf = 1.0/(2.*pi.*Carray(idx).*Rarray(idx));
            fprintf(fileID,"* branch %d, centre frequency = %03.6e,\n", idx-N_l, thisf);

            fprintf(fileID, "C%d cpe1 n%d %03.17e\n", idx, idx, Carray(idx));
            fprintf(fileID, "R%d cpe2 n%d %03.17e\n", idx, idx, Rarray(idx));
        end

        fprintf(fileID, "CtermHF cpe1 cpe2 %03.17e\n", Cterminating);
        fprintf(fileID, "RtermLF cpe1 cpe2 %03.17e\n", Rterminating);

        fprintf(fileID,"\n.ends\n");

        fclose(fileID);
    end
end

end


