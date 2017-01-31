function names = fieldnames(this, varargin)
% this function is a way to provide a list of public member variables
% only and its format. It differs from "display" function, which
% provides a list of both public and private member variables and
% their contents.

% $Author:: kzhu                                      $
% $Rev:: 1487                                         $
% $Date: 2011-02-11 01:56:43 -0500 (Fri, 11 Feb 2011) $

if nargin == 1
    names = {};
else
    switch varargin{1}
    case '-full'
        names = {};
    case '-possible'
        names = {};
        error('Unsupported call to fieldnames');
    end
end
