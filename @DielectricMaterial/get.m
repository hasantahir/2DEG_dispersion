function varargout = get(this, index)
% $Author:: kzhu                                       $
% $Rev:: 1557                                          $
% $Date:: 2011-03-09 19:37:04 -0500 (Wed, 09 Mar 2011) $

if nargin == 1
    if nargout == 0
        disp(struct(this(1)));
    else
        varargout = cell(1,max([1, nargout]));
        varargout{1} = struct(this(1));
    end
    return;
end

% if index is a string, we will allow special access
called_by_name = ischar(index);

% the set switch below needs a substruct
if called_by_name
    index = substruct('.', index);
end

found = false;  % didn't find it in the public section

% special/reserved member variables, not strictly public
if ~found && called_by_name
    found = true;
    switch index(1).subs
      case 'epsilon_r'
        if isempty(this) varargout = {};
        else varargout = {this.epsilon_r}; end
      case 'sigma_e'
        if isempty(this) varargout = {};
        else varargout = {this.sigma_e}; end
      case 'mu_r'
        if isempty(this) varargout = {};
        else varargout = {this.mu_r}; end
      case 'sigma_m'
        if isempty(this) varargout = {};
        else varargout = {this.sigma_m}; end
      case 'display_func'
        if isempty(this) varargout = {};
        else varargout = {this.display_func}; end
      otherwise
        found = false;  % didn't find it in the special section
    end
end

if ~found
    error(['??? Reference to non-existent field ' index(1).subs '.']);
end

if length(varargout) > 1 & nargout <= 1
    if iscellstr(varargout) || any([cellfun('isempty', varargout)])
        varargout = {varargout};
    else
        try
            varargout = {[varargout{:}]};
        catch
            varargout = {varargout};
        end
    end
end
