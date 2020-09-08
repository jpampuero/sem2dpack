function cmap(n,cmin,cmax,fneg,fpos);
  
  switch nargin
   case 0
    n = 256;
    cmin = -1;
    cmax =  1;
    fneg = 0.5;
    fpos = 0.5;
   case 3
    fneg = 0.5;
    fpos = 0.5;
  end
  
  m = zeros(n,3);
  
  b = [0 0 1]; % blue
  w = [1 1 1]; % white
  r = [1 0 0]; % red
  g = [0 1 0]; % green
  k = [0 0 0]; % black
  e = [0 0.5 0]; % dark green
  y = [1 1 0]; % yellow
  c = [0 1 1]; % cyan
  %g = 0.9*w; % gray

  %c1 = b; % negative end
  %c2 = 0.5*b; % part-way to negative end
  %c3 = w; % middle
  %c4 = 0.5*r; % part-way to positive end
  %c5 = r; % positive end

  %c1 = k; % negative end
  %c2 = b; % part-way to negative end
  %c3 = w; % middle
  %c4 = r; % part-way to positive end
  %c5 = k; % positive end

  c1 = c; % negative end
  c2 = b; % part-way to negative end
  c3 = w; % middle
  c4 = r; % part-way to positive end
  c5 = y; % positive end

  % white-red-black one-sided scheme
  %c1 = k; % negative end
  %c2 = r; % part-way to negative end
  %c3 = w; % middle

  if cmin==0
    
    % one-sided
    
    for i=1:3
      %m(:,i) = interp1([0 cmax*fpos cmax],[c3(i) c4(i) c5(i)],linspace(0,cmax,n),'linear');
      m(:,i) = interp1([0 cmax*fpos cmax],[c3(i) c2(i) c1(i)],linspace(0,cmax,n),'linear');
    end
    
  else
  
    % two-sided
    
    for i=1:3
      m(:,i) = interp1([cmin cmin*fneg 0 cmax*fpos cmax],[c1(i) c2(i) c3(i) c4(i) c5(i)],linspace(cmin,cmax,n),'linear');
    end
    
  end

  colormap(m)