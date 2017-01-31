function display(this, display_name)
% This file is class-independent. It displays the content of an
% object. The private member variable "display_func" must be present
% inside the class MATLAB does not provide a good built-in display for
% each classes, so we need this function for each class.

% $Author:: kzhu                                       $
% $Rev:: 1487                                          $
% $Date:: 2011-02-11 01:56:43 -0500 (Fri, 11 Feb 2011) $

if nargin < 2
    % assign 'ans' if inputname(1) empty
    display_name = inputname(1);
    if isempty(display_name)
        display_name = 'ans';
    end
end

% check whether display_func has a value
% if it has a value feval the value to get the display
DisplayFunc = cell(size(this));
[DisplayFunc{:}] = get(this, 'display_func');
use_standard_view = cellfun('isempty', DisplayFunc(:));
if all(use_standard_view(:))
    standard_view(this, display_name);
else
    for k = 1:length(this(:))
        if use_standard_view(k)
            standard_view(this(k), display_name);
        else
            indexed_display_name = sprintf('%s(%d)', display_name, k);
            feval(get(this(k), 'display_func'), this(k), indexed_display_name);
        end
    end    
end

% --------------------------
function standard_view(this, display_name)
if ~isempty( ...
        [strfind(display_name, '.') ...
        strfind(display_name, '(') ...
        strfind(display_name, '{')])
    display_name = 'ans';
end
% handle a scalar vs array object
% note: [] this jumps to else
if length(this) == 1  % scalar case
    % use eval to assign public structure into display_name variable
    eval([display_name ' = struct(this);']);
    % use eval to call display on the display_name structure
    eval(['display(' display_name ');']);
else  % array case
    % use eval to assign this into display_name variable
    eval([display_name ' = this;']);
    % use eval to call builtin display for size info
    eval(['builtin(''display'', ' display_name ');']);
    % still need to display variable names explicitly
    disp('    with public member variables:');
    % get list of public names with fieldname 
    names = fieldnames(this);
    % loop over the name list and display
    for name = names'
        disp(['        ' name{1}]);
    end
    % display extra line if loose
    if strcmp(get(0, 'FormatSpacing'), 'loose')
        disp(' ');
    end
end

% --------------------------
function developer_view(this, display_name)
disp('-----  Public Member Variables  -----');
full_display(struct(this), display_name);
disp('.....  Private Member Variables  .....');
full_display(this, display_name, true);
