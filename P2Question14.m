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
Loglike = zeros(N + 1, N + 1);
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
Loglike(j, l) = Loglike(j, l) - (1/2)*(DLth - DL(i)).^2/errDL(i).^2;
end
end
end
%chisq=-2*Loglike
%chisq=chisq-min(chisq)
figure(1)
xlabel("w0")
ylabel("Omega M")
%colormap(w0,Omm,chisq)
contour(w0, Omm, Loglike, [max(max(Loglike)) - 5.9, max(max(Loglike)) - 3.1, max(max(Loglike)) - 1.15])
H01=68;
Loglike1 = zeros(N + 1, N + 1);
for j = 1 : N + 1
for l = 1 : N + 1
Loglike1(j, l) = 0;
for i = 1 : 580
zgrid = 0 : z(i)/100 : z(i);
Int = 0;
for k = 1 : 101
Int = Int + (z(i)/100) / sqrt(Omm(j) * (1 + zgrid(k))^3 + ((Oml(j)).*(1+zgrid(k)).^(3*(1+w0(l)))));

end
DLth1 = (c/H01)*(1 + z(i)) * Int;
Loglike1(j, l) = Loglike1(j, l) - (1/2)*(DLth1 - DL(i)).^2/errDL(i).^2;
end
end
end
%chisq1=-2*Loglike1
%chisq1=chisq1-min(chisq1)
figure(1)
xlabel("w0")
ylabel("Omega M")
%colormap(w0,Omm,chisq1)
contour(w0, Omm, Loglike1, [max(max(Loglike1)) - 5.9, max(max(Loglike1)) - 3.1, max(max(Loglike1)) - 1.15])
H02=68;
Loglike2 = zeros(N + 1, N + 1);
for j = 1 : N + 1
for l = 1 : N + 1
Loglike2(j, l) = 0;
for i = 1 : 580
zgrid = 0 : z(i)/100 : z(i);
Int = 0;
for k = 1 : 101
Int = Int + (z(i)/100) / sqrt(Omm(j) * (1 + zgrid(k))^3 + ((Oml(j)).*(1+zgrid(k)).^(3*(1+w0(l)))));

end
DLth2 = (c/H02)*(1 + z(i)) * Int;
Loglike2(j, l) = Loglike2(j, l) - (1/2)*(DLth2 - DL(i)).^2/errDL(i).^2;
end
end
end
%chisq2=-2*Loglike2
%chisq2=chisq2-min(chisq2)
figure(1)
xlabel("w0")
ylabel("Omega M")
%colormap(w0,Omm,chisq2)
contour(w0, Omm, Loglike2, [max(max(Loglike2)) - 5.9, max(max(Loglike2)) - 3.1, max(max(Loglike2)) - 1.15])
print -dpng P2Q14