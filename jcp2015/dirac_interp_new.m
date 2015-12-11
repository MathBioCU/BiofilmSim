function val=dirac_interp_new(r)

%can switch between mex and non-mex file
%to switch to non-mex file, change to val=dirac_interp1(r); below


% val=(dirac_interp1(abs(r-2))+dirac_interp1(abs(r))+dirac_interp1(abs(r+2)));
if (isempty(r)==0)
    val=dirac_interp_mex(r);
else
    val=r;
end

end