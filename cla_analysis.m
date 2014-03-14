%å—ƒ‚Ì—g—ÍŒXÎ‚ğ‹‚ß‚é
%—ƒŒ^’è‹`Œã‚Ì«”\‰ğÍ

clear cla_lst_data

defad_vel=9.5;	%‰ğÍ‹@‘¬
ang_adj=0;	%‘Œ¸æ‚è•t‚¯Šp

cla_lst_gamma=Hex_gamma;

if(c!=1)
	disp("regenerting cla_lst_gamma");
	cla_lst_gamma(:,1)=cla_lst_data(:,1).*lst_fdata(:,2).*cla_lst_comp_vel(:,1)./2;
endif

%‚’¼•ûŒü—U“±‘¬“x
disp("recalculating Vn_i.");
for i=1:(wld_c+we_div)
	for j=1:(wld_c+we_div)
		cla_lst_Vn_i_d(j,i)=Q_ij(i,j)*cla_lst_gamma(j,1)./2;
	endfor
	cla_lst_Vn_i(i,1)=sum(cla_lst_Vn_i_d(:,i));
	%printf("%d/100\n",i);
endfor

%—U“±Šp“x
disp("calculating i_ang.");
for i=1:(wld_c+we_div)
	cla_lst_i_ang(i,1)=atan(cla_lst_Vn_i(i,1)/defad_vel)*180/pi;
endfor

%‡¬‘¬“xŒvZ
disp("calculating comp_vel");
for i=1:(wld_c+we_div)
	cla_lst_comp_vel(i,1)=sqrt(defad_vel^2+cla_lst_Vn_i(i,1)^2);
endfor

%—U“±R—Í
disp("calculating di_d.");
for i=1:(wld_c+we_div)
	cla_lst_di_d(i,1)=cla_lst_gamma(i,1)*cla_lst_Vn_i(i,1)*line_e_d(i,1);
endfor

%‹ß—•”
for i=1:5
	cla_lst_di_d(i*dw_divn-1:i*dw_divn,1)=cla_lst_di_d(i*dw_divn-2,1);
	cla_lst_di_d(i*dw_divn+1:i*dw_divn+2,1)=cla_lst_di_d(i*dw_divn+3,1);			
endfor


cla_Di_lst=2*a_density*sum(cla_lst_di_d(:,1));

disp("generating cla_lst_data");


errc=1;
for n=1:we_div+wld_c
	printf("%d,",n);
	cla_lst_re(n,1)=lst_fdata(n,2)*cla_lst_comp_vel(n,1)/nu;
	[bufx foil err]=xfoil([x_foil lst_ydata(:,n)],f_mtang(n,1)-cla_lst_i_ang(n,1)+ang_adj,cla_lst_re(n,1),0);
	if(err==0)
		cla_lst_data(n,:)=[bufx.CL bufx.CD bufx.Cm];
	else
		errlist(errc,1)=n;
		errc++;
	endif
end

if errc!=1
	for n=1:errc-1
		ncc=errlist(n,1)
		cla_lst_data(ncc,:)=[interp1(nang_cp(ncc-1:2:ncc+1,1),cla_lst_data(ncc-1:2:ncc+1,1),nang_cp(ncc,1),"spline") interp1(nang_cp(ncc-1:2:ncc+1,1),cla_lst_data(ncc-1:2:ncc+1,2),nang_cp(ncc,1),"spline") interp1(nang_cp(ncc-1:2:ncc+1,1),cla_lst_data(ncc-1:2:ncc+1,3),nang_cp(ncc,1),"spline")];
	endfor
endif

disp("");
disp("generating cla_L_lst_d");

for i=1:we_div
	cla_L_lst_d(i,1)=1./2.*a_density.*lst_fdata(i,2).*(cla_lst_comp_vel(i,1)^2).*(cla_lst_data(i,1).*cos(cla_lst_i_ang(i,1)./180.*pi)-cla_lst_data(i,2).*sin(cla_lst_i_ang(i,1)./180.*pi)).*cos(line_e_ang(i,1)).*line_e_d(i,1);

endfor

for i=we_div:we_div+wld_c
	cla_L_lst_d(i,1)=1./2.*a_density.*lst_fdata(i,2).*(cla_lst_comp_vel(i,1)^2).*(cla_lst_data(i,1).*cos(cla_lst_i_ang(i,1)./180.*pi)-cla_lst_data(i,2).*sin(cla_lst_i_ang(i,1)./180.*pi)).*cos(line_e_ang(i,1)).*line_e_d(i,1);

endfor

%‹ß—•”
for i=1:5
	cla_L_lst_d(i*dw_divn-1:i*dw_divn,1)=cla_L_lst_d(i*dw_divn-2,1);
	cla_L_lst_d(i*dw_divn+1:i*dw_divn+2,1)=cla_L_lst_d(i*dw_divn+3,1);			
endfor

cla_L_lst=sum(cla_L_lst_d)*2/9.8;

for i=1:we_div
	cla_Dp_lst_d(i,1)=1./2.*a_density.*lst_fdata(i,2).*(cla_lst_comp_vel(i,1)^2).*(cla_lst_data(i,2).*cos(cla_lst_i_ang(i,1)./180.*pi)).*line_e_d(i,1);
	
endfor

for i=we_div:we_div+wld_c
	cla_Dp_lst_d(i,1)=1./2.*a_density.*lst_fdata(i,2).*(cla_lst_comp_vel(i,1)^2).*(cla_lst_data(i,2).*cos(cla_lst_i_ang(i,1)./180.*pi)).*line_e_d(i,1);
	
endfor

%‹ß—•”
for i=1:5
	cla_Dp_lst_d(i*dw_divn-1:i*dw_divn,1)=cla_Dp_lst_d(i*dw_divn-2,1);
	cla_Dp_lst_d(i*dw_divn+1:i*dw_divn+2,1)=cla_Dp_lst_d(i*dw_divn+3,1);			
endfor

cla_Dp_lst=sum(cla_Dp_lst_d)*2;