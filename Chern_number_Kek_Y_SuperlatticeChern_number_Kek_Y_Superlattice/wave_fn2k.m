   function y=wave_fn2k(zi,A0,At,v0,vt,k,th,hb)


T1=sqrt((A0.^2+k.^2.*v0.^2.*hb.^2).*(At.^2+k.^2.*vt.^2.*hb.^2));
T2k=sqrt(A0.^2+At.^2+k.^2.*v0.^2.*hb.^2+k.^2.*vt.^2.*hb.^2-2.*T1);

term1=-v0./vt;
term2=((T2k-(A0+At).*zi).*(A0.^2-T1+A0.*T2k.*zi+k.^2.*v0.^2.*hb.^2))...
    ./(k.^2.*T2k.*v0.*vt.*hb.^2);

y1=term1+term2;
y2=-exp(1i.*th).*(-At.^2-(k.^2.*vt.^2.*hb.^2)+T1+At.*zi.*T2k)...
    ./(k.*vt.*hb.*T2k);
y3=-exp(-1i.*th).*(T1-A0.*(A0+T2k.*zi)-k.^2.*v0.^2.*hb.^2)...
    ./(k.*v0.*hb.*T2k);
y4=1;

A=sqrt(y1.*y1'+y2.*y2'+y3.*y3'+y4.*y4');

yy=[y1;y2;y3;y4];

 y=yy./A;
