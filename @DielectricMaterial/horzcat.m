function this = horzcat(varargin)
% A Comprehensive Guide to Object Oriented Programming in MATLAB
%   Chapter 13 example cStar::horzcat
%   (c) 2005 Andy Register

mismatched = varargin(~cellfun('isclass', varargin, mfilename('class')));
if ~isempty(mismatched)
    error('MATLAB:UnableToConvert', ...
        ['Conversion to ' mfilename('class') ' from ' ...
        class(mismatched{1}) ' is not possible.']);
end

this = builtin(mfilename, varargin{:});