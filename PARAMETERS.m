classdef PARAMETERS
    
    properties
        
        % Cold & Hot parameters
        f;
        amu;
        Z;
        B_func;
        den_m3_func;
        
        % Hot only parameters
        
        T_eV_func;
        k_par;
        
    end
    
    methods
        
        function p = PARAMETERS(f,amu,Z,B_func,den_m3_func,T_eV_func,k_par)
            
            phys = constants();
            
            c = phys.('c');
            me_amu = phys.('me_amu');
            zi = complex(0,1);
            
            % default case
            
            p.amu = {me_amu, 2};
            p.Z   = {-1,1};
            
            p.B_func = @(x) x.*0 + 1.2;
            p.den_m3_func = {@(x) x.*1e19, @(x) x.*1e19};
            
            % default case - hot only parameters
            
            p.T_eV_func = {@(x) x.*0 + 0.001, @(x) x.*0 + 50};
            p.k_par = 20;
            
            % overwrite with user provided parameters
            
            if nargin>0
                if ~isempty(f)
                    p.f = f;
                end
            end
            if nargin>1
                if ~isempty(amu)
                    p.amu = amu;
                end
            end
            if nargin>2
                if ~isempty(Z)
                    p.Z = Z;
                end
            end
            if nargin>3
                if ~isempty(B_func)
                    p.B_func = B_func;
                end
            end
            if nargin>4
                if ~isempty(den_m3_func)
                    p.den_m3_func = den_m3_func;
                end
            end
            if nargin>5
                if ~isempty(T_eV_func)
                    p.T_eV_func = T_eV_func;
                end
            end
            if nargin>6
                if ~isempty(k_par)
                    p.k_par = k_par;
                end
            end
            
        end
        
        
    end
    
end