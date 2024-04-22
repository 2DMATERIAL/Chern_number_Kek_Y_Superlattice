  function y=wave_fn4k(zi,A0,At,v0,vt,k,th,hb)

T1=sqrt((A0.^2+k.^2.*v0.^2.*hb.^2).*(At.^2+k.^2.*vt.^2.*hb.^2));
T2kp=sqrt(A0.^2+At.^2+k.^2.*v0.^2.*hb.^2+k.^2.*vt.^2.*hb.^2+2.*T1);

term1=-v0./vt;
term2=((T2kp-(A0+At).*zi).*(A0.^2+T1+A0.*T2kp.*zi+k.^2.*v0.^2.*hb.^2))...
    ./(k.^2.*T2kp.*v0.*vt.*hb.^2);

y1=term1+term2;
y2=-exp(1i.*th).*(-At.^2-(k.^2.*vt.^2.*hb.^2)-T1+At.*zi.*T2kp)...
    ./(k.*vt.*hb.*T2kp);
y3=-exp(-1i.*th).*(T1+A0.*(A0+T2kp.*zi)+k.^2.*v0.^2.*hb.^2)...
    ./(k.*v0.*hb.*T2kp);
y4=1;

A=sqrt(y1.*y1'+y2.*y2'+y3.*y3'+y4.*y4');

yy=[y1;y2;y3;y4];

 y=yy./A;
