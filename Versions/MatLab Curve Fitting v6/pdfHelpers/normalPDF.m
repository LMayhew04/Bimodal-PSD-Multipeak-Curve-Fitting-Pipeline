function y = normalPDF(x,mu,sigma,gain)
    if nargin<4, gain = 1; end
    y = gain * (1/(sigma*sqrt(2*pi))) .* exp(-(x-mu).^2./(2*sigma.^2));
end