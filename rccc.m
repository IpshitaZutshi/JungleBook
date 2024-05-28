function [aopt phi0 rho R p]=rccc(lindat, phidat, abound, da)
iNan = isnan(lindat) | isnan(phidat);
phidat = phidat(~iNan);
lindat = lindat(~iNan);


%% phidat in radiants
%%
% starting points of maximization
Nrep=ceil(abs(abound(2)-abound(1))/da);
astart=min(abound);


Rfunc = @(a) -abs(nanmean(exp(1i*(phidat-2*pi*a*lindat))));
%Rfunc = @(a) sqrt(nanmean(sin(phidat-2*pi*a*lindat))^2+nanmean(sin(phidat-2*pi*a*lindat))^2);

aopt=nan;
R=-10^10;
for n=1:Nrep
    a0=astart+(n-1)*da;
    a1=astart+n*da;
    
    [atmp, rtmp]=fminbnd(Rfunc,a0,a1);
    
    if(-rtmp>R)
        R=-rtmp;
        aopt=atmp;
    elseif (-rtmp==R)
        aopt=cat(1,aopt,atmp);
    end
end

phis = zeros(0);
rhos = zeros(0);
ps = zeros(0);
numA = length(aopt);
if numA > 1
    if numA == Nrep
        aopt = nan;
        phi0 = nan;
        rho = nan;
        p = nan;
        R = nan;
    else
        for a = 1:numA
            %phase offset
            v=mean(exp(1i*(phidat-2*pi*aopt(a)*lindat)));
            phi0=angle(v);
            %phi0 = angle(mean(exp(1i*phidat)) / mean(exp(1i*2*pi*aopt*lindat)));
            %
            
            theta=angle(exp(2*pi*1i*abs(aopt(a))*(lindat)));
            
            thmean=angle(sum(exp(1i*theta)));
            phmean=angle(sum(exp(1i*phidat)));
            
            sthe=sin(theta-thmean);
            sphi=sin(phidat-phmean);
            
            c12=sum(sthe.*sphi);
            c11=sum(sthe.*sthe);
            c22=sum(sphi.*sphi);
            rho=c12/sqrt(c11*c22);
            
            lam22=mean(sthe.^2.*sphi.^2);
            lam20=mean(sphi.^2);
            lam02=mean(sthe.^2);
            tval=rho*sqrt(length(lindat)*lam20*lam02/lam22);
            
            p=1-erf(abs(tval)/sqrt(2));
            phis = cat(1,phis,phi0);
            rhos = cat(1,rhos,rho);
            ps = cat(1,ps,p);
        end
        if isnan(min(ps))
            aopt = nan;
            phi0 = nan;
            rho = nan;
            p = nan;
            R = nan;
        else
            
            pInd = ps==min(ps);
            aInd = abs(aopt)==min(abs(aopt(pInd)));
            
            %             p = ps(ind);
            %             if length(p) > 1
            %                 ind = find(ps==min(ps),1);
            %             end
            
            aopt = aopt(aInd);
            if min(abs(abound-aopt'))<0.01
                aopt = nan;
                phi0 = nan;
                rho = nan;
                p = nan;
                R = nan;
            else
                phi0 = phis(aInd);
                rho = rhos(aInd);
                p = ps(aInd);
            end
        end
    end
else
    if min(abs(abound-aopt))<0.01
        aopt = nan;
        phi0 = nan;
        rho = nan;
        p = nan;
        R = nan;
    else
        
        v=mean(exp(1i*(phidat-2*pi*aopt*lindat)));
        phi0=angle(v);
        %phi0 = angle(mean(exp(1i*phidat)) / mean(exp(1i*2*pi*aopt*lindat)));
        
        theta=angle(exp(2*pi*1i*abs(aopt)*(lindat)));
        
        thmean=angle(sum(exp(1i*theta)));
        phmean=angle(sum(exp(1i*phidat)));
        
        sthe=sin(theta-thmean);
        sphi=sin(phidat-phmean);
        
        c12=sum(sthe.*sphi);
        c11=sum(sthe.*sthe);
        c22=sum(sphi.*sphi);
        rho=c12/sqrt(c11*c22);
        
        lam22=mean(sthe.^2.*sphi.^2);
        lam20=mean(sphi.^2);
        lam02=mean(sthe.^2);
        tval=rho*sqrt(length(lindat)*lam20*lam02/lam22);
        
        p=1-erf(abs(tval)/sqrt(2));
    end
end

