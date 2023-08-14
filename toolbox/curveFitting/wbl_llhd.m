function llhd = wbl_llhd(a, b, x, r)

    p = wblcdf(x, a, b);
    llhd = zeros(size(p));

    for idx = 1:length(p)
        if r(idx) == 0
            llhd(idx) = 1 - p(idx);
        else
            llhd(idx) = p(idx);
        end
    end

    llhd = sum(log(llhd + 1e-6));

end