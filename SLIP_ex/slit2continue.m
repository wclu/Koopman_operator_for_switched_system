function Y = slit2continue(X)
x = X(:,1);
for i = 1:length(x)-1
   if (x(i+1)-x(i)) < 0
       x(i+1:end) = x(i+1:end) -x(i+1) + x(i);
   end
end

Y = X;
Y(:,1) = x;

end