% convergence analysis for implicit3

% load data into variable n1
load('-ascii','n1.dat');
% load data into variable n2 ...
load('-ascii','n2.dat');
load('-ascii','n3.dat');
load('-ascii','n4.dat');
load('-ascii','n5.dat');

err = zeros(4,1);
% compute successive error
err(1) = max(abs(n1-n2));
err(2) = max(abs(n2-n3));
err(3) = max(abs(n3-n4));
err(4) = max(abs(n4-n5));

nt = 10*[1;2;4;8;16];

loglog(nt,err,'ko-')
fit = polyfit(log(nt),log(err),1);
slope = fit(1)
