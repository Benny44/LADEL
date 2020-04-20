classdef ladel < handle
    %Matlab interface for the LADEL C library
    %
    %ladel methods
    %
    %
    %
    %
    
    properties (SetAccess = private, Hidden = true)
        ncol;
    end
    
    methods
        
        %% Constructor - Create a new solver instance
        function this = ladel(ncol)
            ladel_mex('init', ncol);
        end
        
        function delete(~)
            ladel_mex('delete');
        end
        
        function factorize(~, M, varargin)
            M = triu(M);
            if nargin == 3
                ordering = varargin{1};
                ladel_mex('factorize', M, ordering);
            else
                ladel_mex('factorize', M);
            end    
        end
        
        function y = dense_solve(~, x)
            y = ladel_mex('solve', x);
        end

    end
    
end

