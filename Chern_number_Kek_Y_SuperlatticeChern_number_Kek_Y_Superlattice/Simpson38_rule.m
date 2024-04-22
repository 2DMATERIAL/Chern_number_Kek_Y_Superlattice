 function y=Simpson38_rule(X,Y) %Y=Y(X)

% X=linspace(0,2,4);
% Y=2.*X.^2;

n=length(X);
sum=0;
h=X(2)-X(1);
for i=1:3:n-3
    sum=sum+Y(i).*h.*3./8;
end
for i=2:3:n-2
    sum=sum+3.*Y(i).*h.*3./8;
end
for i=3:3:n-1
    sum=sum+3.*Y(i).*h.*3./8;
end
for i=4:3:n
    sum=sum+Y(i).*h.*3./8;
end
y=sum;