function hScopes = create16QAMScopes

persistent hScope
if isempty(hScope)
    hScope = My16QAMScopes;
end
hScopes = hScope;

end

