%混合率ゼロの翼型のCl*Cを得る関数（翼端用）
function clc=zerofoil_def_wl(cl_c,aoa,vel,adata)
nu=1.604*10^(-5);
c_list=0.1:0.1:0.5;
ang_list=-3:1:7;
re_list=c_list.*vel/nu;

disp("generating zerofoil");
cl_c

errct=1;
for rec=1:length(c_list)
	printf("%d,",rec);
	% for angc=1:11
	% 	angc
	% 	% angc_cl(angc,1)=Re_spline(re_list(rec),6,1,ang_list(angc),2,adata);
	% 	angc_cl(angc,1)=xfoil(adata,ang_list(angc),re_list(rec),0).CL;
		
	% endfor
	
	% a_cl(1,rec)=spline(ang_list,angc_cl(:,1),aoa);
	% a_cl(1,rec)=xfoil(adata,aoa,re_list(rec),0).CL;
	%あるRE数におけるあるAOAのCL
	[bufx foil err]=xfoil(adata,aoa,re_list(rec),0);
	if(err==0)
		a_cl(rec,1)=bufx.CL;
	else
		errlistt(errct,1)=rec
		errct++;
	endif
	
endfor

if errct!=1
	for n=1:errct-1
		ncct=errlistt(n,1)
		a_cl(ncct,1)=interp1(re_list(1,ncct-1:2:ncct+1),a_cl(ncct-1:2:ncct+1,1),re_list(1,ncct));
	endfor
endif

disp("")

%CLをもとにC＊CLを計算する
for i=1:length(c_list)
	c_cl_o(1,i)=a_cl(i,1)*c_list(1,i);
endfor

%入力したc*clにあうCをsplineでもとめる
% length(c_cl_o)
% length(c_list)

clc=interp1(c_cl_o,c_list,cl_c,"spline","extrap")

endfunction
