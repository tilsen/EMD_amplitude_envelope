function [metrics] = envm_metrics_batch(X,par,varargin)
%extracts measures for a set of chunks
%X: cell array of chunks
%par: structure with processing parameters (specify at least the audio
%       sample rate, i.e. par.Fs = 44100)

%some sanity checks on input
X = cellfun(@(c){c(:)'},X);
sizes = cell2mat(cellfun(@(c){size(c)},X));

if any(sizes(:,1)==0)
    fprintf('ERROR: missing audio. Aborting\n');
    return;
end
if any(sizes(:,1)>1)
    fprintf('ERROR: invalid input: audio must be vectors. Aborting\n');
    return;
end

p = inputParser;

def_bandpass = [400 nan]; %upper band edge is determined from Fs
def_lowpass = 10;
def_ds = 100; %downsampling factor
def_nfft = 2048; 
def_sm_Hz = 1;
def_powerratio_freq_bins = [1 3.5 10];
def_centroid_freq_bins = [1 10];
def_max_imf = 3;
def_imf_freq_bounds = [0 13.16]; %frequency at which magnitude response of 4th order butterworth is -10dB
def_tukeywin_param = 0.2;
def_sift_relative_tol = 0.1;
def_edge_null = 0.1;
def_imf_freq_exclusion_prctile = 99;

addRequired(p,'X',@(x)iscell(X));
addRequired(p,'par',@(x)isstruct(x) & ismember('Fs',fieldnames(x)));
addParameter(p,'bandpass',def_bandpass);
addParameter(p,'lowpass',def_lowpass);
addParameter(p,'ds',def_ds);
addParameter(p,'nfft',def_nfft);
addParameter(p,'sm_Hz',def_sm_Hz);
addParameter(p,'powerratio_freq_bins',def_powerratio_freq_bins);
addParameter(p,'centroid_freq_bins',def_centroid_freq_bins);
addParameter(p,'max_imf',def_max_imf);
addParameter(p,'tukeywin_param',def_tukeywin_param);
addParameter(p,'sift_relative_tol',def_sift_relative_tol);
addParameter(p,'edge_null',def_edge_null);
addParameter(p,'imf_freq_bounds',def_imf_freq_bounds);
addParameter(p,'imf_freq_exclusion_prctile',def_imf_freq_exclusion_prctile);
addParameter(p,'verbose',false);

parse(p,X,par,varargin{:});

res = p.Results;

ff = fieldnames(par);
for i=1:length(ff)
    res.(ff{i}) = par.(ff{i});
end

if isnan(res.bandpass(2))
    res.bandpass(2) = (res.Fs/2) - round(0.01*res.Fs);
end

res.envFs = res.Fs/res.ds;

if res.verbose
    ff = setdiff(fieldnames(res),{'par','X'});
    ffstr = pad(ff,'left');
    fprintf('extracting metrics with parameters:\n');
    for i=1:length(ff)
        if isstruct(res.(ff{i})) || iscell(res.(ff{i})), continue; end
        fprintf('\t%s: %s\n',ffstr{i},regexprep(num2str(res.(ff{i})),'\s+',' '));
    end
end

%get the envelopes, power spectra, imfs, and imf frequencies
ENV = cell(size(X));
t_ENV = ENV; smPSD = ENV; F = ENV;

%downsampled envelopes
for i=1:length(X)
    [ENV{i},t_ENV{i}] = envm_band_energy(X{i},res);
end
if res.verbose
    fprintf('extracted %d envelopes, downsampled to %1.1f\n',length(X),res.envFs); end

ENV = cellfun(@(c){c-mean(c)},ENV);
ENV = cellfun(@(c){c/max(abs(c))},ENV);
if res.verbose
    fprintf('normalized %d envelopes\n',length(X)); end

ENV = cellfun(@(c){c.*tukeywin(length(c),res.tukeywin_param)'},ENV);
if res.verbose
    fprintf('applied tukey window to %d envelopes, with param: %1.2f\n',length(X),res.tukeywin_param); end

res.audioFs = res.Fs;
res.Fs = res.envFs;

for i=1:length(ENV)
    [smPSD{i},F{i}] = envm_smoothed_psd(ENV{i},res);
end

if res.verbose
    fprintf('obtained %d smoothed power spectra, with smoothing bandwidth: %1.2f\n',length(X),res.sm_Hz); end

IMF = cellfun(@(c){emd(c,'SiftRelativeTolerance',res.sift_relative_tol,'MaxNumIMF',res.max_imf)},ENV);

n_imfs = cellfun(@(c)size(c,2),IMF);

if res.verbose
    fprintf('empirical mode decomposition summary (max imfs %d):\n\n',res.max_imf);
    T = array2table(tabulate(n_imfs),'VariableNames',{'NumIMFs','Count','Perc'});
    T = T(T.Count>0,:);    
    disp(T);
end

for i=1:length(IMF)
    [HS{i}, ~, t_IMF{i},W{i}] = hht(IMF{i},res.envFs);
end
if res.verbose
    fprintf('obtained instantaneous IMF instantaneous frequencies using Hilbert-Huang transform\n');
end

%1. replace edge values with nan
if res.edge_null>0
    for i=1:length(W)
        W{i}(1:round(res.envFs*res.edge_null)) = nan;
        W{i}(end-round(res.envFs*res.edge_null)+1:end) = nan;
    end
    if res.verbose
        fprintf('replaced IMF frequency edge values (%1.3f s) with nan\n',res.edge_null);
    end
end

%2. impose bounds on frequencies
if ~isnan(res.imf_freq_bounds(1))
    for i=1:length(W)
        W{i}(W{i}<res.imf_freq_bounds(1)) = nan;
    end
    if res.verbose
        fprintf('imposed lower bound on IMF frequencies: (%1.3f Hz)\n',res.imf_freq_bounds(1));
    end    
end
if ~isnan(res.imf_freq_bounds(2))
    for i=1:length(W)
        W{i}(W{i}>res.imf_freq_bounds(2)) = nan;
    end
    if res.verbose
        fprintf('imposed upper bound on IMF frequencies: (%1.3f Hz)\n',res.imf_freq_bounds(2));
    end    
end

%3. for each IMF, exclude values above some percentile
if ~isempty(res.imf_freq_exclusion_prctile) && ~isnan(res.imf_freq_exclusion_prctile)

    wthresh = nan(1,res.max_imf);
    for i=1:res.max_imf %loop over IMF #s
        imfW = [];
        for j=1:length(W)
            if (size(W{j},2)>=i)
                imfW = [imfW; W{j}(:,i)];
            end
        end
        wthresh(i) = prctile(imfW,res.imf_freq_exclusion_prctile);
    end

    counts = zeros(length(W),res.max_imf);
    samples = counts;
    for i=1:length(W)
        for j=1:size(W{i},2)
            w = W{i}(:,j);
            ixs = w>wthresh(j);
            counts(i,j) = sum(ixs);
            samples(i,j) = sum(~isnan(w));
            W{i}(ixs,j) = nan;
        end
    end

    rates = sum(counts)./sum(samples);

    if res.verbose
        fprintf('percentile-based exclusion summary:\n\n');
        T = array2table(rates',"RowNames",arrayfun(@(c){sprintf('imf%d',c)},1:res.max_imf), ...
            "VariableNames",{'exclusionRate'});
        disp(T);
    end  
end

%get metrics
for i=1:length(smPSD)
    SPEC(i,1) = envm_psd_metrics(smPSD{i},F{i},res);
end
SPEC = struct2table(SPEC);

%get emd metrics
for i=1:length(IMF)
    EMD(i,1) = envm_emd_metrics(num2cell(IMF{i},1),num2cell(W{i},1));
end
EMD = struct2table(EMD);

%add nan values for missing imfs
vecf = {'pow_imf' 'mu_w' 'var_w' 'sd_w'};
for i=1:height(EMD)
    for j=1:length(vecf)
        n_vals = length(EMD.(vecf{j}){i});
        if n_vals<res.max_imf
            EMD.(vecf{j}){i} = [EMD.(vecf{j}){i} nan(1,res.max_imf-n_vals)];
        end        
    end
end

for i=1:length(vecf)
    EMD.(vecf{i}) = cell2mat(EMD.(vecf{i}));
    sz = size(EMD.(vecf{i}),2);
    for j=1:sz
        EMD.([vecf{i} num2str(j)]) = EMD.(vecf{i})(:,j);
    end
    EMD.(vecf{i}) = [];
end

metrics = [SPEC EMD];

end