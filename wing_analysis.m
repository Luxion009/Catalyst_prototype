%翼型定義後の性能解析

clear lst_data errlist

do

	%定義したデータをコントロールポイントに保管する

	lst_fdata(:,1)=nang_cp(:,1);
	lst_fdata(:,2)=interp1(nang_line(:,1),f_data(:,2),nang_cp(:,1),"extrap");
	lst_fdata(:,3)=interp1(nang_line(:,1),f_data(:,3),nang_cp(:,1),"extrap");

	for i=1:50
		lst_ydata(:,i)=y_foil(:,1);
	endfor

	for i=2:w_div_c-1
		for j=1:dw_divn
			lst_ydata(:,(i-1)*dw_divn+j)=y_foil(:,i-1)*lst_fdata(i,3)+y_foil(:,i)*(1-lst_fdata(i,3));
		endfor
	endfor

	for j=1:dw_divn
		lst_ydata(:,4*dw_divn+j)=y_foil(:,4);
	endfor

	for j=we_div+1:we_div+wld_c
		lst_ydata(:,j)=y_foil(:,4);
	endfor

	lst_gamma=Hex_gamma;
	
	%m=input("calc count=");

	c=0;
	do
	c++;
	
		if(c!=1)
			disp("regenerting lst_gamma");
			lst_gamma(:,1)=lst_data(:,1).*lst_fdata(:,2).*lst_comp_vel(:,1)./2;
		endif
		
		%垂直方向誘導速度
		disp("recalculating Vn_i.");
		for i=1:(wld_c+we_div)
			for j=1:(wld_c+we_div)
				lst_Vn_i_d(j,i)=Q_ij(i,j)*lst_gamma(j,1)./2;
			endfor
			lst_Vn_i(i,1)=sum(lst_Vn_i_d(:,i));
			%printf("%d/100\n",i);
		endfor
		
		%誘導角度
		disp("calculating i_ang.");
		for i=1:(wld_c+we_div)
			lst_i_ang(i,1)=atan(lst_Vn_i(i,1)/d_vel)*180/pi;
		endfor

		%合成速度計算
		disp("calculating comp_vel");
		for i=1:(wld_c+we_div)
			lst_comp_vel(i,1)=sqrt(d_vel^2+lst_Vn_i(i,1)^2);
		endfor

		%誘導抗力
		disp("calculating di_d.");
		for i=1:(wld_c+we_div)
			lst_di_d(i,1)=lst_gamma(i,1)*lst_Vn_i(i,1)*line_e_d(i,1);
		endfor

		%近似部
		for i=1:5
			lst_di_d(i*dw_divn-1:i*dw_divn,1)=lst_di_d(i*dw_divn-2,1);
			lst_di_d(i*dw_divn+1:i*dw_divn+2,1)=lst_di_d(i*dw_divn+3,1);			
		endfor


		Di_lst=2*a_density*sum(lst_di_d(:,1));

		disp("generating lst_data");
		

		errc=1;
		for n=1:we_div+wld_c
			lst_re(n,1)=lst_fdata(n,2)*lst_comp_vel(n,1)/nu;
			[bufx foil err]=xfoil([x_foil lst_ydata(:,n)],f_mtang(n,1)-lst_i_ang(n,1),lst_re(n,1),0);
			if(err==0)
				printf("%d,",n);
				lst_data(n,:)=[bufx.CL bufx.CD bufx.Cm];
			else
				printf("e,");
				errlist(errc,1)=n;
				errc++;
			endif
		end

		if errc!=1
			for n=1:errc-1
				ncc=errlist(n,1)
				lst_data(ncc,:)=[interp1(nang_cp(ncc-1:2:ncc+1,1),lst_data(ncc-1:2:ncc+1,1),nang_cp(ncc,1)) interp1(nang_cp(ncc-1:2:ncc+1,1),lst_data(ncc-1:2:ncc+1,2),nang_cp(ncc,1)) interp1(nang_cp(ncc-1:2:ncc+1,1),lst_data(ncc-1:2:ncc+1,3),nang_cp(ncc,1))];
			endfor
		endif
		buff=input("")

		disp("");
		disp("generating L_lst_d");
		
		for i=1:we_div
			L_lst_d(i,1)=1./2.*a_density.*lst_fdata(i,2).*(lst_comp_vel(i,1)^2).*(lst_data(i,1).*cos(lst_i_ang(i,1)./180.*pi)-lst_data(i,2).*sin(lst_i_ang(i,1)./180.*pi)).*cos(line_e_ang(i,1)).*line_e_d(i,1);

		endfor

		for i=we_div:we_div+wld_c
			L_lst_d(i,1)=1./2.*a_density.*lst_fdata(i,2).*(lst_comp_vel(i,1)^2).*(lst_data(i,1).*cos(lst_i_ang(i,1)./180.*pi)-lst_data(i,2).*sin(lst_i_ang(i,1)./180.*pi)).*cos(line_e_ang(i,1)).*line_e_d(i,1);

		endfor

		%近似部
		for i=1:5
			L_lst_d(i*dw_divn-1:i*dw_divn,1)=L_lst_d(i*dw_divn-2,1);
			L_lst_d(i*dw_divn+1:i*dw_divn+2,1)=L_lst_d(i*dw_divn+3,1);			
		endfor

		L_lst=sum(L_lst_d)*2/9.8

		for i=1:we_div
			Dp_lst_d(i,1)=1./2.*a_density.*lst_fdata(i,2).*(lst_comp_vel(i,1)^2).*(lst_data(i,2).*cos(lst_i_ang(i,1)./180.*pi)).*line_e_d(i,1);
			
		endfor

		for i=we_div:we_div+wld_c
			Dp_lst_d(i,1)=1./2.*a_density.*lst_fdata(i,2).*(lst_comp_vel(i,1)^2).*(lst_data(i,2).*cos(lst_i_ang(i,1)./180.*pi)).*line_e_d(i,1);
			
		endfor

		%近似部
		for i=1:5
			Dp_lst_d(i*dw_divn-1:i*dw_divn,1)=Dp_lst_d(i*dw_divn-2,1);
			Dp_lst_d(i*dw_divn+1:i*dw_divn+2,1)=Dp_lst_d(i*dw_divn+3,1);			
		endfor

		Dp_lst=sum(Dp_lst_d)*2

		printf("Di_lst");
		disp(Di_lst);
		D_lst=Dp_lst+Di_lst
		W_LST=(D_lst+6)*d_vel/0.85/0.9
		c
		
		figure(16)
			plot(f_data(:,1),ed_Hex_gamma(:,1),"-;design gamma;",lst_fdata(:,1),lst_gamma(:,1),"-;lst gamma;");
			xlabel("y[m]");
			grid on;
			cd("design_data")
			cd(pj_name)
			print("lst_gamma.png",'-dpng','-r100')
			cd(now_work)

		figure(17)
			plot(f_data(:,1),ed_Hex_Vn_i(:,1),"-;ed_Hex_Vn_i;",lst_fdata(:,1),lst_Vn_i(:,1),"-;lst_Vn_i;");
			xlabel("y[m]");
			grid on;
			cd("design_data")
			cd(pj_name)
			print("lst_Vn_i.png",'-dpng','-r100')
			cd(now_work)

		sign=input("ok:1 more:0 ---");
		
	until(sign==1)

	sign=input("ok:1 return:0 ---");
	
until(sign==1)