function out=tLST(wing,vel,mangle,nu,aero_data)
%-----ó‹µ‚ğ“ü—Í‚µ‚Ä—g—ÍŒW”‚ğ•Ô‚·
%-----wing‚Í[ƒXƒpƒ“•ûŒüˆÊ’uC—ƒŒ·’·]‚Ì‡


%-----ŠeƒŒƒCƒmƒ‹ƒY”‚É‚Ä—ë—g—ÍŠp‚Æ—g—ÍŒXÎ‚ğ‹‚ß‚é
%foil=1;
%for i=1:12
	%poly1(i,:)=polyfit(aero_data{foil,i}(1:25,1),aero_data{foil,i}(1:25,2),1);
	%Cl0(i,1)=roots(poly1(i,:));%deg
	%relist(i,1)=50000*i;
%end
	%Clalpha(:,1)=poly1(:,1);%1/deg
	%-----Å‰‚ÌŒ}‚¦Šp,ƒŒƒCƒmƒ‹ƒY”
	aoa(1:rows(wing),1)=mangle;
	re(1:rows(wing),1)=vel.*wing(:,2)./nu;
	
disp("calculating gamma")

for c=1:2
	
	c
	%Cl(:,1)=spline(relist,Clalpha,re).*(aoa-spline(relist,Cl0,re));
	%ˆÀS‚ÆM—Š‚ÌRe_spline
	for i=1:rows(wing)
		%printf('%.3d ',i)
		Cl(i,1)=Re_spline(re(i,1),1,1,aoa(i,1),2,aero_data);
	end
		gamma(:,1)=1./2.*vel.*wing(:,2).*Cl;
		
		dgamma_dy(:,1)=gradient(gamma,wing(:,1));
		
		for i=1:rows(wing)
			for n=1:rows(wing)
				if n==i
					dwh(n,1)=dgamma_dy(n,1)/2;
				else
					dwh(n,1)=dgamma_dy(n,1)/abs(wing(n,1)-wing(i,1))/pi/(-4);
				end
				
				%‘ÎÌ—ƒ‚Ì‚½‚ß‚Ì—U“±‘¬“xŒ©Ï‚à‚è
				%dwh(n+rows(wing),1)=dgamma_dy(n,1)/abs(wing(n,1)+wing(i,1))/pi/(-4);
				
			end
			dw(i,1)=trapz(wing(:,1),dwh);
		end
		
		aoa(:,1)=mangle-atan(dw./vel).*180./pi;
		re(:,1)=sqrt(vel.^2+dw.^2).*wing(:,2)./nu;
				
end		

for i=1:rows(wing)
	%printf('%.3d ',i)
	Cl(i,1)=Re_spline(re(i,1),1,1,aoa(i,1),2,aero_data);
end
			
out=[Cl,re,aoa,dw];
			
endfunction