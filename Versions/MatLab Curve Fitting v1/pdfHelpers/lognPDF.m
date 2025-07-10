function y = lognPDF(x,mu,sigma,gain)
    if nargin<4, gain = 1; end
    y = gain * lognpdf(x,mu,sigma);
end