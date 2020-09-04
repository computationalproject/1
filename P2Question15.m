load supernova.dat
z=supernova(:,2);
dismodul=supernova(:,3);
errmodul=supernova(:,4);
DL=10.^((dismodul/5)-5);
errDL=(log(10)/5).*errmodul.*DL;
c=3e5;
H0=70;
N = 5;
Omm = 0 : 1/N : 1;
Oml = 1 - Omm;
w0 = -2 : 1/(N/2) : 0;
Loglike2D = zeros(N + 1, N + 1);
%Do we need to add a grid for uncertainty in Omm or w0?
%If so, add another for loop including all the previous ones
for j = 1 : N + 1
for l = 1 : N + 1
Loglike(j, l) = 0;
for i = 1 : 580
zgrid = 0 : z(i)/100 : z(i);
Int = 0;
for k = 1 : 101
  #Thanks, Guillerme
Int = Int + (z(i)/100) / sqrt(Omm(j) * (1 + zgrid(k))^3 + ((Oml(j)).*(1+zgrid(k)).^(3*(1+w0(l)))));

end
DLth = (c/H0)*(1 + z(i)) * Int;
Loglike2D(j, l) = Loglike2D(j, l) - (1/2)*(DLth - DL(i)).^2/errDL(i).^2;
end
end
end
Likelihood2D=exp(Loglike2D);
Likelihood1D=sum(Likelihood2D);
Loglike1D=log(Likelihood1D);
%chisq=-2*Loglike1D
%chisq=chisq-min(chisq)
figure(2)
xlabel("w0")
ylabel("Omega M")
%colormap(w0,Omm,chisq)
contour(w0, Omm, Loglike, [max(max(Loglike)) - 5.9, max(max(Loglike)) - 3.1, max(max(Loglike)) - 1.15])
print -dpng P2Q13