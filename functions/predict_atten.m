function predict_atten( data,samprate,dT,dtstar,logA,alp,fign )
%predict_atten( data,samprate,dT,dtstar,logA,alp,fign )
%  

if nargin<7 || isempty(fign)
    fign = 3;
end

maxq = mindex(-dtstar);

dtstar_do = dtstar(maxq) - dtstar;

odata = zeros(size(data));
for is = 1:size(data,2)
    odata(:,is) = attenuate( data(:,is),samprate,dtstar_do(is),-dT(is),(10.^(-logA(is))),alp ); % attenuate each to max, undo effects of dT and gain
end

tt = ([0:size(data,1)-1]-size(data,1)/2)/samprate;

figure(fign), clf, hold on
plot(tt,data,':')
plot(tt,odata)
axis([-15,15,0,1])


end

