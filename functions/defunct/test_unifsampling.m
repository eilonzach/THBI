niter = 20000;
N = 20;
parm = 'dt';


PARMS

% profile clear
% profile on
x = ([eval(sprintf('%s.min',parm)):eval(sprintf('%s.disc',parm)):eval(sprintf('%s.max',parm))])';
y = zeros(size(x));
for ii = 1:niter
    vals = sample_prior(N,parm,1);
    y = y + hist(vals,x)';
end
% profile viewer
bar(x,y)
