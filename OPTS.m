classdef OPTS
    
    properties
        
        use_root_finder = false;
        use_cold_eps = true;
        use_cold_det = false;
        num_init_kper = 5;        
        kper_min = -1000;
        kper_max = +1000;
        num_points = 100;
        
    end
    
    methods
        
        function opts = OPTS(varargin)
            
            if nargin == 0
            else
                input_parser = inputParser();
                input_parser.KeepUnmatched = true;
                
                addOptional(input_parser,'use_root_finder',opts.use_root_finder,@islogical);
                addOptional(input_parser,'use_cold_eps',opts.use_cold_eps,@islogical);               
                addOptional(input_parser,'use_cold_det',opts.use_cold_det,@islogical);
                addOptional(input_parser,'num_init_kper',opts.num_init_kper,@isnumeric);
                addOptional(input_parser,'num_points',opts.num_points,@isnumeric);
                addOptional(input_parser,'kper_min',opts.kper_min,@isnumeric);
                addOptional(input_parser,'kper_max',opts.kper_max,@isnumeric);

                parse(input_parser,varargin{:}{:})

                fields = fieldnames(input_parser.Results);
                
                num_fields = numel(fields);
                
                for f=1:num_fields
                    opts.(fields{f}) = input_parser.Results.(fields{f});
                end
                
            end
            
        end
        
    end
    
end
