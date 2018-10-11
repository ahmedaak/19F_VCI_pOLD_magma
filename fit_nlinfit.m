

function [params,err] = fit_nlinfit(y,x,model,debug)

%% Data Fitting
% using Matlab's nlinfit
% INPUT:
% 
% y    : matrix, data (last dim for dependent variable) 
% x    : row vector, dependent variable
% est  : row vector of estimates 
% model:
%     T1_2p             : Two parameter T1 regrowth     S(t) = S0 * (1-exp(-t/T1)
%     T1_4p             : Four parameter T1 regrowth    S(t)=(a+b*exp(-t/T1))/(1+c*exp(-t/T1)); Kingsley et al. JMR 1999
%     R1_2p             : Two parameter R1 regrowth     S(t) = S0 * (1-exp(-t*R1)
%     R1_4p             : Four parameter R1 regrowth    S(t)=(a+b*exp(-t*R1))/(1+c*exp(-t*R1)); Kingsley et al. JMR 1999
%     T2_2p             : Two   parameter T2 decay      S(t) = S0 * exp(-t/T2)
%     T2_3p             : Three parameter T2 decay      S(t) = S0 * exp(-t/T2) + C
%     sin_abs           : Two parameter flip angle map  S(FA)= abs(S0*sin(c*FA))
%     ernstSignalVarFa  : Two parameter FLASH signal 
%                         neglecting T2 decay           S(FA)= S0*sin(FA)/(1-cos(FA)*exp(-TR/T1)) !TR needs to be defined manually!
%     vsiCpmgExpSe_4p   : Four parameter fit for CPMG 
%                         signal after iron oxide 
%                         injection                     S(t)= S0 * exp(-R2*t)exp(-const*t^(exponent))
%     vsiDeltaR2Exponent: 2 parameter decay for CPMG
%                         dependeny of DeltaR2          S(n)= DeltaR2*n^(exponent)
%     T1_invRec         : T1 inversion recovery         S(TI)=a+b*exp(-TI/T1)

% OUTPUT:
% 
% params: matrix of fitted values, same size as y, bad values (fit did not
%         converge) marked as NaN
% err   : relative error of fitted parameters !!is not yet implemented
%         right!!
%
% last modified by P Boehm-Sturm July 07 2015
%%

if nargin < 4, debug = 0;end
if nargin < 3, model = 'T2_3p';end    
if nargin < 2, error('More input needed'), end

dims = size(y);

if dims(end) ~= numel(x), error('dimension mismatch'), end

y = double(reshape(y,[],numel(x)));

if debug
    scrsz = get(0,'ScreenSize');    
    h = figure('Position',[1 scrsz(4)-scrsz(4)/4 scrsz(3)/4 scrsz(4)/4]);
end

% Estimates
if strcmp(model,'T1_2p')
    estT1=y(:,end).*x(1)./y(:,1);
    est = [estT1 y(:,end)];
elseif strcmp(model,'T1_4p')
    estT1=y(:,end).*x(1)./y(:,1);
    estS0=y(:,end);
    est = [estT1 estS0 -estS0 zeros(size(estT1))]; % for perfect RF pulses a=-b=S0 and c=0
elseif strcmp(model,'R1_2p')
    estT1=y(:,end).*x(1)./y(:,1);
    estR1=1./estT1;
    est = [estR1 y(:,end)];
elseif strcmp(model,'R1_4p')
    estT1=y(:,end).*x(1)./y(:,1);
    estR1=1./estT1;
    estS0=y(:,end);
    est = [estR1 estS0 -estS0 zeros(size(estT1))]; % for perfect RF pulses a=-b=S0 and c=0
elseif strcmp(model,'T2_2p')
    estT2 = (x(2)-x(1))./(log(y(:,1))./log(y(:,2)));
    est = [estT2 y(:,1) ];
elseif strcmp(model,'T2_3p')   
    estT2 = (x(2)-x(1))./(log(y(:,1))./log(y(:,2)));
    est = [estT2 y(:,1)  y(:,end)];
elseif strcmp(model,'R2_2p')
    estR2 = 1./((x(2)-x(1))./(log(y(:,1))./log(y(:,2))));
    est = [estR2 y(:,1) ];
elseif strcmp(model,'R2_3p')   
    estR2 = 1./((x(2)-x(1))./(log(y(:,1))./log(y(:,2))));
    est = [estR2 y(:,1)  y(:,end)];
elseif strcmp(model,'T2_4p')   
    %estT2 = (x(2)-x(1))./(log(y(:,1))./log(y(:,2)));
    %est = [estT2 y(:,1)./2 estT2 y(:,1)./2];
    est = [100 y(:,1) 50 y(:,1)];
elseif strcmp(model,'sin_abs')   
   est_A=max(y,[],2);
   %est_omega=pi./max(x);
   est_omega=ones(size(est_A))./10;
   est = [est_A est_omega];
elseif strcmp(model,'ernstSignalVarFa')   
   estM0=max(y,[],2);
   estT1=ones(size(estM0)).*100;           % estimated T1 is last number (watch units! usually ms)
   est = [estT1 estM0];
elseif strcmp(model,'vsiCpmgExpSe_4p')   
   estM0=y(:,1);
   estT2 = (x(2)-x(1))./(log(y(:,1))./log(y(:,2)));
   estR2=estT2.^(-1);
   estExp=ones(size(estM0))./3;
   est= [estR2 estExp estR2./1 estM0];
   %est_omega=pi./max(x);
elseif strcmp(model,'vsiDeltaR2Exponent')
    estDeltaR2=y(:,1);
    estExp=ones(size(estDeltaR2)).*(2/3);
    est = [estExp estDeltaR2];
elseif strcmp(model,'T1_invRec')
    %esta=y(:,end);
    %estb=esta.*(-2);
    %estT1=x(1)./(-y(:,1)./esta+1); %watch out with signs in inv rec from magn. images
    %est = [estT1 esta estb];
    estT1=repmat(x(floor(length(x)./2)),[size(y,1),1]);
    estInvEff=repmat(1,[size(y,1),1]);
    est=[estT1 y(:,end) estInvEff];
end


options = statset('FunValCheck','off');

params = zeros(size(y,1),size(est,2));
err    = zeros(size(y,1),size(est,2));

t=(0:x/50:1.2*max(x)).';

warning off;

for k = 1:size(y,1)

    if y(k,1) ~=0
        if strcmp(model,'T1_2p')
            [out,res,J] = nlinfit(x,y(k,:),@T1_2p,est(k,:)',options);
            ci = nlparci(out,res,J);
            cf=T1_2p(out,t);          % fitted curve
            trueci=ci;
            trueci(1,1)=ci(1,2);      % flip interval values for T1 as lower/higher T1 leads to higher/lower signal values for plot of confidence bounds in signal domain
            trueci(1,2)=ci(1,1);
            uci = T1_2p(trueci(:,1),t);   % upper ci
            lci = T1_2p(trueci(:,2),t);   % lower ci
        elseif strcmp(model,'T1_4p')
            [out,res,J] = nlinfit(x,y(k,:),@T1_4p,est(k,:)',options);
            ci = nlparci(out,res,J);
            cf=T1_4p(out,t);          % fitted curve
            trueci=ci;
            trueci(1,1)=ci(1,2);      % flip interval values for T1 as lower/higher T1 leads to higher/lower signal values for plot of confidence bounds in signal domain
            trueci(1,2)=ci(1,1);
            trueci(4,1)=ci(4,2);      % flip interval values for c as lower/higher c leads to higher/lower signal values for plot of confidence bounds in signal domain
            trueci(4,2)=ci(4,1);
            uci = T1_4p(trueci(:,1),t);   % upper ci
            lci = T1_4p(trueci(:,2),t);   % lower ci
        elseif strcmp(model,'R1_2p')
            [out,res,J] = nlinfit(x,y(k,:),@R1_2p,est(k,:)',options);
            ci = nlparci(out,res,J);
            cf=T1_2p(out,t);          % fitted curve
            trueci=ci;
            uci = T1_2p(trueci(:,1),t);   % upper ci
            lci = T1_2p(trueci(:,2),t);   % lower ci
        elseif strcmp(model,'R1_4p')
            [out,res,J] = nlinfit(x,y(k,:),@R1_4p,est(k,:)',options);
            ci = nlparci(out,res,J);
            cf=T1_4p(out,t);          % fitted curve
            trueci=ci;
            trueci(4,1)=ci(4,2);      % flip interval values for c as lower/higher c leads to higher/lower signal values for plot of confidence bounds in signal domain
            trueci(4,2)=ci(4,1);
            uci = T1_4p(trueci(:,1),t);   % upper ci
            lci = T1_4p(trueci(:,2),t);   % lower ci
        elseif strcmp(model,'T2_2p')
            [out,res,J] = nlinfit(x,y(k,:),@T2_2p,est(k,:)',options);
            ci = nlparci(out,res,J);
            trueci=ci;
            cf=T2_2p(out,t);          % fitted curve
            uci = T2_2p(ci(:,1),t);   % upper ci
            lci = T2_2p(ci(:,2),t);   % lower ci
        elseif strcmp(model,'T2_3p')
            [out,res,J] = nlinfit(x,y(k,:),@T2_3p,est(k,:)',options);
            ci = nlparci(out,res,J);
            trueci=ci;
            cf=T2_3p(out,t);          % fitted curve
            uci = T2_3p(ci(:,1),t);   % upper ci
            lci = T2_3p(ci(:,2),t);   % lower ci
        elseif strcmp(model,'R2_2p')
            [out,res,J] = nlinfit(x,y(k,:),@R2_2p,est(k,:)',options);
            ci = nlparci(out,res,J);
            trueci=ci;
            trueci(1,1)=ci(1,2);      % flip interval values for R2 as lower/higher R2 leads to higher/lower signal values for plot of confidence bounds in signal domain
            trueci(1,2)=ci(1,1);
            cf=R2_2p(out,t);          % fitted curve
            uci = R2_2p(trueci(:,1),t);   % upper ci
            lci = R2_2p(trueci(:,2),t);   % lower ci
        elseif strcmp(model,'R2_3p')
            [out,res,J] = nlinfit(x,y(k,:),@R2_3p,est(k,:)',options);
            ci = nlparci(out,res,J);
            trueci=ci;
            trueci(1,1)=ci(1,2);      % flip interval values for R2 as lower/higher R2 leads to higher/lower signal values for plot of confidence bounds in signal domain
            trueci(1,2)=ci(1,1);
            cf=R2_3p(out,t);          % fitted curve
            uci = R2_3p(trueci(:,1),t);   % upper ci
            lci = R2_3p(trueci(:,2),t);   % lower ci
            
        elseif strcmp(model,'T2_4p')
            [out,res,J] = nlinfit(x,y(k,:),@T2_4p,est(k,:)',options);
            ci = nlparci(out,res,J);
            trueci=ci;
            cf=T2_4p(out,t);          % fitted curve
            uci = T2_4p(ci(:,1),t);   % upper ci
            lci = T2_4p(ci(:,2),t);   % lower ci
        elseif strcmp(model,'sin_abs')
            [out,res,J] = nlinfit(x,y(k,:),@sin_abs,est(k,:)',options);
            ci  = nlparci(out,res,J);
            cf  =sin_abs(out,t);          % fitted curve
            trueci=ci;
            trueci(2,:)=fliplr(ci(2,:));
            uci = sin_abs(trueci(:,1),t);   % upper ci
            lci = sin_abs(trueci(:,2),t);   % lower ci
         elseif strcmp(model,'ernstSignalVarFa')
            [out,res,J] = nlinfit(x,y(k,:),@ernstSignalVarFa,est(k,:)',options);
            ci  = nlparci(out,res,J);
            cf  =ernstSignalVarFa(out,t);          % fitted curve
            trueci=ci;
            trueci(1,:)=fliplr(ci(1,:));
            uci = ernstSignalVarFa(trueci(:,1),t);   % upper ci
            lci = ernstSignalVarFa(trueci(:,2),t);   % lower ci
        elseif strcmp(model,'vsiCpmgExpSe_4p')
            [out,res,J] = nlinfit(x,y(k,:),@vsiCpmgExpSe_4p,est(k,:)',options);
            ci = nlparci(out,res,J);
            cf=vsiCpmgExpSe_4p(out,t);          % fitted curve
            trueci=zeros(4,2);                  % function is not monotenous with all parameters
            trueci(1,:)=ci(1,:);
            trueci(2,:)=fliplr(ci(2,:));
            trueci(3,:)=fliplr(ci(3,:));
            trueci(4,:)=fliplr(ci(4,:));
            uci = vsiCpmgExpSe_4p(trueci(:,1),t);   % upper ci
            lci = vsiCpmgExpSe_4p(trueci(:,2),t);   % lower ci  
        elseif strcmp(model,'vsiDeltaR2Exponent')
            [out,res,J] = nlinfit(x,y(k,:),@vsiDeltaR2Exponent,est(k,:)',options);
            ci = nlparci(out,res,J);
            trueci=ci;
            cf = vsiDeltaR2Exponent(out,t);          % fitted curve
            uci = vsiDeltaR2Exponent(ci(:,1),t);   % upper ci
            lci = vsiDeltaR2Exponent(ci(:,2),t);   % lower ci 
        elseif strcmp(model,'T1_invRec')
            [out,res,J] = nlinfit(x,y(k,:),@T1_invRec,est(k,:)',options);
            ci = nlparci(out,res,J); % cis not displayed in function plot correctly, fix later!
            trueci=ci;
            cf = T1_invRec(out,t);          % fitted curve
            uci = T1_invRec(trueci(:,1),t);   % upper ci
            lci = T1_invRec(trueci(:,2),t);   % lower ci 
        end
        
        absErr = abs(ci(:,2)-ci(:,1))/2;  %estimated error
        relErr = absErr./abs(out).*100;  % relative error [%]


        if debug
            figure(h)
            plot(x,y(k,:),'bo');hold on;
            plot(t,cf,'r-'); 
            %plot(t,uci,'k--');
            %plot(t,lci,'k-.');
            hold off;
            grid on
            axis([0  max(x) 0 max(y(k,:)) ])
            axis tight;
            xlabel('Time [ms]')
            ylabel('Signal [a.u.]')
            legend('Measurement','Fit','Upper','Lower','Location','Best')
            drawnow  
            if debug ~= 1
                pause
            end
        end 
    
        if relErr(1) < 100            
            params(k,:) = out;    
            err(k,:) = relErr;
        else
            params(k,:) = NaN;
            err(k,:)    = NaN;
        end   
        
    end  
end

warning on;

if size(params,1) > 1
    params = reshape(params,[dims(1:end-1),size(est,2)]);
    err    = reshape(err,[dims(1:end-1),size(est,2)]);
end

function y = T1_2p(x,t) 
y = x(2)*(1-exp(-t./x(1))); % (T1, M0)
  if ~all(isfinite(y(:)))
      % disp('Something bad happend in MODELFUN!!')
      y(:) = 0;
  end

function y = T1_4p(x,t) 
y = (x(2)+x(3).*exp(-t./x(1)))./(1+x(4).*exp(-t./x(1))); % (T1, M0)
  if ~all(isfinite(y(:)))
      % disp('Something bad happend in MODELFUN!!')
      y(:) = 0;
  end
function y = R1_2p(x,t) 
y = x(2)*(1-exp(-t.*x(1))); % (T1, M0)
  if ~all(isfinite(y(:)))
      % disp('Something bad happend in MODELFUN!!')
      y(:) = 0;
  end

function y = R1_4p(x,t) 
y = (x(2)+x(3).*exp(-t.*x(1)))./(1+x(4).*exp(-t.*x(1))); % (T1, M0)
  if ~all(isfinite(y(:)))
      % disp('Something bad happend in MODELFUN!!')
      y(:) = 0;
  end
  
function y = T2_2p(x,t) 
y = x(2)*exp(-t./x(1)); % (T2, M0)
  if ~all(isfinite(y(:)))
      % disp('Something bad happend in MODELFUN!!')
      y(:) = 0;
  end

function y = T2_3p(x,t) 
y = x(2)*exp(-t./x(1))+x(3);% (T2, M0, noise)
  if ~all(isfinite(y(:)))
      % disp('Something bad happend in MODELFUN!!')
      y(:) = 0;
  end

function y = R2_2p(x,t) 
y = x(2)*exp(-t.*x(1)); % (R2, M0)
  if ~all(isfinite(y(:)))
      % disp('Something bad happend in MODELFUN!!')
      y(:) = 0;
  end

function y = R2_3p(x,t) 
y = x(2)*exp(-t.*x(1))+x(3);% (R2, M0, noise)
  if ~all(isfinite(y(:)))
      % disp('Something bad happend in MODELFUN!!')
      y(:) = 0;
  end

  
function y = T2_4p(x,t) 
y = x(2)*exp(-t./x(1))+x(4)*exp(-t./x(3));% (T21,M01,T22,M02)
  if ~all(isfinite(y(:)))
      % disp('Something bad happend in MODELFUN!!')
      y(:) = 0;
end
  
function y = sin_abs(x,t)
y = abs(x(1)*sin(x(2).*t));  % (A,omega)
if ~all(isfinite(y(:)))
 % disp('Something bad happend in MODELFUN!!')
 y(:) = 0;
end  

function y=ernstSignalVarFa(x,t)                                                                 % TR needs to be defined manually
y=x(2).*(sin(t).*(1-exp(-15./x(1)))./(1-cos(t).*exp(-15./x(1)))); % (T1, M0 including T2 decay)
if ~all(isfinite(y(:)))
 % disp('Something bad happend in MODELFUN!!')
 y(:) = 0;
end


function y = vsiCpmgExpSe_4p(x,t)
y = x(4).*exp(-t.*x(1)).*exp(-(t.^(x(2))).*x(3));  % (R2 decay, exponent, CPMG contribution, M0)
if ~all(isfinite(y(:)))
 % disp('Something bad happend in MODELFUN!!')
 y(:) = 0;
end

function y = vsiDeltaR2Exponent(x,t)
y = x(2)*(t.^x(1));  % (exponent, DeltaR2)
if ~all(isfinite(y(:)))
 % disp('Something bad happend in MODELFUN!!')
 y(:) = 0;
end

function y = T1_invRec(x,t)
y = abs(x(2).*(1-2.*x(3).*exp(-t./x(1))));  % (a+b exp(-TI/T1))
if ~all(isfinite(y(:)))
 % disp('Something bad happend in MODELFUN!!')
 y(:) = 0;
end
