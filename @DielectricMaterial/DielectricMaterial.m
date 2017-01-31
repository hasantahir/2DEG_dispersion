function this = DielectricMaterial(varargin)
% this = DielectricMaterial(epsilon_r, sigma_e, mu_r, sigma_m)
%
% this = DielectricMaterial(epsilon_r,sigma_e) assumes a non-magnetic
% material: mu_r =1 and sigma_m = 0
%
% The class is defined by following the convention of [Register2007].
%
%

% $Author:: kzhu                                               $
% $Rev:: 1488                                                  $  
% $Date:: 2011-02-14 17:13:57 -0500 (Mon, 14 Feb 2011)         $

  class_name = mfilename('class');
  default_this = struct([]);
  default_this(1).display_func = 'developer_view';
  switch nargin
    case 2
      default_this(1).epsilon_r = cell2mat(varargin(1));
      default_this(1).sigma_e   = cell2mat(varargin(2));
      default_this(1).mu_r      = 1;
      default_this(1).sigma_m   = 0;
      this = class(default_this, class_name);
    case 4
      default_this(1).epsilon_r = cell2mat(varargin(1));
      default_this(1).sigma_e   = cell2mat(varargin(2));
      default_this(1).mu_r      = cell2mat(varargin(3));
      default_this(1).sigma_m   = cell2mat(varargin(4));
      this = class(default_this, class_name);
    otherwise
      error('??? Invalid argument to constructor.');
  end
  
  
  
  
