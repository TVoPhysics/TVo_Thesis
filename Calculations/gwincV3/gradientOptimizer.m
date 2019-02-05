function x0 = gradientOptimizer(scoremode,start,variate,stepsize)
%function x0 = gradientOptimizer(scoremode,start,variate,stepsize)
% gradientOptimizer 
%
% scoremode:
%       0: optimize for NS/NS
%       1: optimize for BH/BH
%       2: stocastic
% start:   vector of start values: 
%            [LASER POWER,...
%             SRC PHASE,...
%             SRM TRANS,...
%             ITM TRANS,...
%             PRM TRANS,...
%             HOMODYN PHASE]
% variate: vector of 0 or 1, specifying whether to change the corresponding value: 
%            [LASER POWER,...
%             SRC PHASE,...
%             SRM TRANS,...
%             ITM TRANS,...
%             PRM TRANS,...
%             HOMODYN PHASE]
% note that TPRM is handeled separately, because it can be simply calculated
% variate(5)=0: use the optimal TPRM calculated by gwinc
% variate(5)=1: use the TPRM from sart(5)
%
% stepsize: stepsize of (is optional):
%            [LASER POWER,...
%             SRC PHASE,...
%             SRM TRANS,...
%             ITM TRANS,...
%             HOMODYN PHASE].
%


global f_seis f_Nyq shotradmode ifo src tprmind stpsz scmode;

try, stepsize;
catch
  stepsize=[1e-2,1e-3,1e-3,1e-3,0,1e-3];
end
maxiter=20;

f_seis = 9; % Seismic Wall Cutoff frequency
f_Nyq = 8192; % Fake Nyquist frequency
tprmind=5;
shotradmode=2;
ifo=IFOModel;
src=SourceModel;
stpsz=stepsize;
scmode=scoremode;
 
% while abs(corr)>errmax
% end

  x0=reshape(start,length(start),1);
  index=find(variate==1);
  for it=1:maxiter
    [y,G,J]=gradientJacobi(x0,variate);
    fprintf('Iteration: %d,\tScore: %g',it,y);
    if norm(G)==0
      break;
    end
    TG=transpose(G);
    dx0=-(G*TG)/(G*J*TG)*TG;
    x0(index)=x0(index)+dx0;
    x0=sanityCheck(x0);
    fprintf('\tPower: %g\tDetuning: %g deg\tTSRC: %g\tTITM: %g',x0(1),x0(2)*180/pi,x0(3),x0(4));
    if length(x0)>tprmind
      fprintf('\tHomo phase: %g',x0(6)/pi*180);
    end
    fprintf('\n');
  end;
  
    shotradmode=3;
    if length(x0)>tprmind
      ifo.Optics.Quadrature.dc=x0(6);
    end;
    switch variate(tprmind)
      case 1,
        [sss,nnn] = gwinc(f_seis,f_Nyq,ifo,src,shotradmode,x0(1),x0(2),x0(3),x0(4),x0(tprmind));
      otherwise,
        [sss,nnn] = gwinc(f_seis,f_Nyq,ifo,src,shotradmode,x0(1),x0(2),x0(3),x0(4));
    end;
  
return



function varargout=gradientJacobi(x,variate);
global f_seis f_Nyq shotradmode ifo src tprmind stpsz scmode;
  vspec=variate;
  vspec(tprmind)=0;
  d=0;
  maxdim=length(x);
  y=callGwinc(x,variate);
  for ii=1:maxdim
    if(vspec(ii))
      d=d+1;
      dx=zeros(size(x));
      dx(ii)=stpsz(ii);
      yp=callGwinc(x+dx,variate);
      ym=callGwinc(x-dx,variate);
      G(d)=(yp-ym)/2/dx(ii);
      if nargout==3
        J(d,d)=(yp-2*y+ym)/(dx(ii).^2);
        for d2=1:(d-1)
          jj=min(find(cumsum(vspec)==d2));
          dx2=zeros(size(x));
   	  dx2(jj)=stpsz(jj);
          ypp=callGwinc(x+dx+dx2,variate);
          ymp=callGwinc(x-dx+dx2,variate);
          ypm=callGwinc(x+dx-dx2,variate);
          ymm=callGwinc(x-dx-dx2,variate);	
          J(d,d2)=(ypp-ymp-ypm+ymm)/4/dx(ii)/dx2(jj);
	  J(d2,d)=J(d,d2);
        end
      end
    end
  end
  switch nargout
    case 1,
      varargout{1}=y;
    case 2,
      varargout{1}=y;
      varargout{2}=G;
    case 3,
      varargout{1}=y;
      varargout{2}=G;
      varargout{3}=J;
  end
return


function y=callGwinc(x,variate)
global f_seis f_Nyq shotradmode ifo src tprmind stpsz scmode;

    if length(x)>tprmind
      ifo.Optics.Quadrature.dc=x(6);
    end;
    switch variate(tprmind)
      case 1,
        [sss,nnn] = gwinc(f_seis,f_Nyq,ifo,src,shotradmode,x(1),x(2),x(3),x(4),x(tprmind));
      otherwise,
        [sss,nnn] = gwinc(f_seis,f_Nyq,ifo,src,shotradmode,x(1),x(2),x(3),x(4));
    end;
    shotradmode=4;
   
    switch scmode
      case 1,
        score=sss.effr0bh;
      case 2,
        score=1/sss.Omega;
      otherwise,
        score=sss.effr0ns;
    end;
    y=score;
return

function x=sanityCheck(x);
global stpsz;
  if x(1)<0
     x(1)=0.1;
  end
  if x(2)<0
     x(2)=0;
  end
  if x(2)>pi/2
     x(2)=pi/2*7/9;
  end
  if x(3)<0
     x(3)=2*stpsz(3);
  end
  if x(3)>1
     x(3)=1-2*stpsz(3);
  end
  if x(4)<0
     x(4)=2*stpsz(4);
  end
  if x(4)>1
     x(4)=1-2*stpsz(4);
  end
  try,
    x(6)=mod(x(6),pi);
  end
return
