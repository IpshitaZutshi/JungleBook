function S1Prime = correctOneOverf(S1,f,useMean)
if nargin<3 || isempty(useMean)
    useMean = 0;
end

nt = size(S1,1);
nf = length(f);
x = [ones(nf,1) log(f)'];
S1Prime = zeros(size(S1));
if useMean
    logS1 = log(mean(S1));
    beta = regress(logS1',x);
    logS1Hat = beta(1)+beta(2)*log(f);
end
logS1 = log(S1);
for t = 1:nt
    if ~useMean
    beta = regress(logS1(t,:)',x);
    logS1Hat = beta(1)+beta(2)*log(f);
    end
    logS1Prime = logS1(t,:)-logS1Hat;
    S1Prime(t,:) = exp(logS1Prime);
end
