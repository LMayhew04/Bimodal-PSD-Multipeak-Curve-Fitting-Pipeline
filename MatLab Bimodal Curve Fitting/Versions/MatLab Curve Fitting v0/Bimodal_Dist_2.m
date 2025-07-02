% Bimodal Distribution Mathematical Model

%import data for variables ParticleDiameter and SizeDistribution1
clearvars -except ParticleDiameter SizeDistribution1
close all;
x=0.001:0.01:500;

%% Normal Distribution 1
sigma1=0.186; %0.19
mu1=0.538; %0.54
pdf1 = (1/(sigma1*sqrt(2*pi))) * exp((-(x-mu1).^2)/(2*sigma1^2));
figure(1)
subplot(3,1,1)
%plot(x,pdf1,'b')
semilogx(x,pdf1,'b')
xlim([0.01 100])
ylabel('Distribution (%)')
legend('Normal Dist.')
grid on

%CDF
fun1 = @(t) exp(-(t-mu1).^2/(2*sigma1^2));
x1=0.001;
for i=1:1:size(x,2)
    int1(i) = integral(fun1,0,x1);
    cdf1(i) = (1/(sigma1*sqrt(2*pi)))*int1(i);
x1=x1+0.01;
end
figure(2)
subplot(3,1,1)
semilogx(x,cdf1)
xlim([0.01 100])
ylabel('CDF (%)')
legend('Normal Dist.')
figure(1)

%% Normal Distribution 2
% sigma2=3.0;
% mu2=8.4;
% pdf2 = 3.03 * (1/(sigma2*sqrt(2*pi))) * exp((-(x-mu2).^2)/(2*sigma2^2));
% subplot(3,1,2)
% %plot(x,pdf2,'--g')
% semilogx(x,pdf2,'--g')
% xlim([0.01 100])
% ylabel('Distribution (%)')
% %legend('Normal Dist.')
% grid on


%% Rayleigh Distribution 2
% b1=5.6;
% pdf3=(x/b1^2) .* exp (-x.^2/(2*b1^2));
% subplot(3,1,2)
% %plot(x,pdf2,'r')
% semilogx(x,pdf3,'r')
% xlim([0.01 100])
% ylabel('Distribution (%)')
% legend('Rayleigh Dist.')
% grid on


%% Lognormal Distribtuion 2
sigma4=0.585; %0.59
mu4=1.936; %1.95
pdf4 = (1./(x.*sigma4*sqrt(2*pi))) .* exp ( (-(log(x)-mu4).^2) / (2*sigma4^2) );
subplot(3,1,2)
%hold on
%plot(x,pdf4,'m')
semilogx(x,pdf4,'m')
xlim([0.01 100])
ylabel('Distribution (%)')
%legend('Normal Dist.', 'Lognormal Dist.')
legend('Lognormal Dist.')
grid on

%CDF
fun4 = @(t) (1./t) .* exp( -(log(t)-mu4).^2 / (2*sigma4^2) );
x4=0.001;
for i=1:1:size(x,2)
    int4(i) = integral(fun4,0,x4);
    cdf4(i) = (1/(sigma4*sqrt(2*pi)))*int4(i);
x4=x4+0.01;
end
figure(2)
subplot(3,1,2)
semilogx(x,cdf4)
xlim([0.01 100])
ylabel('CDF (%)')
legend('Lognormal Dist.')
figure(1)


%% Bimodal Distribution
% 2 Normal Distributions
% gain1=0.8/2;
% gain2=4.8/0.2;
% pdfsum=gain1*pdf1+gain2*pdf2;

% 1st Normal Distribution, 2nd Rayleigh Distribution
% gain1=0.8/2;
% gain3=4.8/0.2;
% pdfsum=gain1*pdf1+gain3*pdf3;

% 1st Normal Distribution, 2nd Lognormal Distribution
gain1=1.418/2; %1.4/2
gain4=7.846/0.2; %7.9/0.2
pdfsum=gain1*pdf1+gain4*pdf4;

%plot the bimodal distribution of the math model
subplot(3,1,3)
%plot(x,pdfsum,'k')
semilogx(x,pdfsum,'k')
xlim([0.01 100])
ylabel('Distribution (%)')
grid on

%CDF
cdfsum = (100/(gain1+gain4))*gain1*cdf1 + (100/(gain1+gain4))*gain4*cdf4;
figure(2)
subplot(3,1,3)
semilogx(x,cdfsum)
xlim([0.01 100])
ylabel('CDF (%)')
legend('Bimodal Dist.')
figure(1)

%plot the acutal data
%import data separately before plotting
subplot(3,1,3)
hold on
%plot(ParticleDiameter,SizeDistribution1,'--r')
semilogx(ParticleDiameter,SizeDistribution1,'--r')
legend('model','actual')
grid on



%% Compute an error quantity between math model and actual measurements

%Step through the actual measured points, find closest math model point,
%compute the vector distance.  Square the vector distances and sum, this
%single value represents the error quantity.  Vector distance is already
%squared and all positive, no need to square.

for n=1:1:size(ParticleDiameter,1)
    for m=1:1:size(x,2)
        err1(n,m)=sqrt( (x(1,m)-ParticleDiameter(n))^2 + (pdfsum(1,m)-SizeDistribution1(n))^2 );
    end
end
clear n m
%now go thru err1 row-by-row and find the minimum error value
Mins = min(err1,[],2);
%sum the error values of final error. quantity is the goodness of fit
Error=sum(Mins);

%now iterate by changing the variables in deterministic way -> HOW?
%Steepest Decent Method? Change variables in what order? None are
%independent!!!

%Initial conditions
sigma1_store(1)=sigma1;
mu1_store(1)=mu1;
sigma4_store(1)=sigma4;
mu4_store(1)=mu4;
gain1_store(1)=gain1;
gain4_store(1)=gain4;
Error_store(1)=Error;

%Divide up the 2 modes for optimization
xcut=1; %x value for "cut-line" dividing mode 1 from mode 2
m1=1;
m2=1;
for n=1:1:size(ParticleDiameter,1)
    if ParticleDiameter(n)<xcut
        ParticleDiameter_mode1(m1)=ParticleDiameter(n);
        SizeDistribution1_mode1(m1)=SizeDistribution1(n);
        m1=m1+1;
    else
        ParticleDiameter_mode2(m2)=ParticleDiameter(n);
        SizeDistribution1_mode2(m2)=SizeDistribution1(n);
        m2=m2+1;
    end
end
%the following max values for actuals and fits give a means to re-optimize
%the gains based on the difference in actual max and fitted max, re-compute
%gain after every change and apply.
[max_mode1, index_mode1]=max(SizeDistribution1_mode1);
[max_mode2, index_mode2]=max(SizeDistribution1_mode2);
[fitmax_mode1, fitindex_mode1]=max(gain1*pdf1);
[fitmax_mode2, fitindex_mode2]=max(gain4*pdf4);
%why not also compute the "exact right" x value for the fit based on the
%known x value for the actuals!!!???
PDmax_mode1=ParticleDiameter_mode1(index_mode1)+0.05;  %added fudge factor of 0.05 added to place peak at right place for actual data, actual data to sparse
PDmax_mode2=ParticleDiameter_mode2(index_mode2);
fitPDmax_mode1=x(fitindex_mode1);
fitPDmax_mode2=x(fitindex_mode2);

%possible path foward for optimization:
%Step 1: change mu's to drive difference between model max-x and
%actual max-x to minimum, must update gains, mu's don't affect each other
%Step 2: change sigma's to drive computed Error to
%minumum, must update gains. sigma's don't affect each other
%Observation: mu and sigma are not independent, sigma moves observed mu

%assume that initial guess is really good
%compute 10 next points in each direction for mu's and sigma's?
clear n m m1 m2

mustep=0.001;
nsteps=40;
for n=1:1:nsteps
    if n > (nsteps/2)
        mustep1=(mustep*-1)+((nsteps/2*mustep)/n);
    else
        mustep1=mustep;
    end
    %normal dist 1
    %sigma1=sigma1_store(1)+n*mustep1;
    sigma1=sigma1_store(1);
    %mu1=mu1_store(1)+n*mustep1;
    mu1=mu1_store(1);
    pdf1 = (1/(sigma1*sqrt(2*pi))) * exp((-(x-mu1).^2)/(2*sigma1^2));
    %lognormal dist 2
    %sigma4=sigma4_store(1)+n*mustep1;
    sigma4=sigma4_store(1);
    %mu4=mu4_store(1)+n*mustep1;
    mu4=mu4_store(1);
    pdf4 = (1./(x.*sigma4*sqrt(2*pi))) .* exp ( (-(log(x)-mu4).^2) / (2*sigma4^2) );
    %bimodal dist
    gain1=1.418/2;
    gain4=7.846/0.2;
    pdfsum=gain1*pdf1+gain4*pdf4;

    for n1=1:1:size(ParticleDiameter,1)
        for m1=1:1:size(x,2)
            err1(n1,m1)=sqrt( (x(1,m1)-ParticleDiameter(n1))^2 + (pdfsum(1,m1)-SizeDistribution1(n1))^2 );
        end
    end
    %now go thru err1 row-by-row and find the minimum error value
    Mins = min(err1,[],2);
    %sum the error values of final error. quantity is the goodness of fit
    Error=sum(Mins);


    %mu1_store(n+1)=mu1;
    %sigma1_store(n+1)=sigma1;
    %mu4_store(n+1)=mu4;
    %sigma4_store(n+1)=sigma4;
    Error_store(n+1)=Error;
end


