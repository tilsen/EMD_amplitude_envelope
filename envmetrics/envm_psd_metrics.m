function [SPEC] = envm_psd_metrics(psd,freqs,par)


for j=1:size(par.powerratio_freq_bins,1)
    freq_bin1_ix = (freqs>=par.powerratio_freq_bins(j,1) & freqs<=par.powerratio_freq_bins(j,2));
    freq_bin2_ix = (freqs>par.powerratio_freq_bins(j,2) & freqs<=par.powerratio_freq_bins(j,3));
    bin1_pow = sum(psd(freq_bin1_ix,:));
    bin2_pow = sum(psd(freq_bin2_ix,:));
    SPEC.(['sbpr_' num2str(j)]) = bin1_pow./bin2_pow;
end

for j=1:size(par.centroid_freq_bins,1)
    freq_bin_ix = (freqs>=par.centroid_freq_bins(j,1) & freqs<=par.centroid_freq_bins(j,2));
    bin_centroid = sum(psd(freq_bin_ix).*freqs(freq_bin_ix)') / sum(psd(freq_bin_ix));
    SPEC.(['scntr_' num2str(j)]) = bin_centroid;
end



end