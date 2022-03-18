function parsave(filename, v, vname)
eval([vname '= v;'])
try
    save(filename, vname, '-append')
catch
    save(filename, vname)
end
end