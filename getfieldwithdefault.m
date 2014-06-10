function val = getfieldwithdefault(struct_data, field_name, default_value)
if isfield(struct_data,field_name)
    val = struct_data.(field_name);
else
    val = default_value;
end
