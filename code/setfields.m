function structure=setfields(structure,varargin)
%SETFIELDS sets the fields of a structure. Does not overwrite existing field values
% example structure = setfields(structure,'name','david','number',1)
% Use structure = setfields([],'name','david','number',1) to create a new structure
c=1;
for i=1:floor(length(varargin)/2)
    f=varargin{c}; v=varargin{c+1};
    if isstruct(v)
        fname=fieldnames(v);
        for j=1:length(fname)
            if ~isfield(structure.(f),fname{j})
                structure.(f).(fname{j})=v.(fname{j});
            end
        end
    else
        if ~isfield(structure,varargin{c})
            structure.(f)=v;
        end
    end
    c=c+2;
end
