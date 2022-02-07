function [EMD] = envm_emd_metrics(imf,w)

for i=1:length(imf)
    EMD.pow_imf(i) = sum(abs(imf{i}));
    EMD.mu_w(i) = nanmean(w{i});
    EMD.var_w(i) = nanvar(w{i});
end

EMD.imf_ratio21 = EMD.pow_imf(2)./EMD.pow_imf(1);

end