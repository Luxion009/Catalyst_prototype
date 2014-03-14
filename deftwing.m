%尾翼定義
%１セクションにつきリブ4本以下
disp("generate h_tail");

do

	%諸元入力
	lhd=4.32	%主翼取付け部から水平尾翼取付け部まで4.52以下でなければならない

	Vh=input("Vh=");
	Sh=Vh*mw_s*cbar/lhd;
	Fh=Sh*lhd^2/mw_s/(cbar^2)
		% lim_hre=input("水平尾翼限界レイノルズ数=");
		% rs_h=input("リブ間(mm)=");
		% brn_h=input("BOX部リブ数=");
		% trn_h=input("テーパー部リブ数=");
		% lim_hchord=lim_hre*nu/d_vel;%限界翼弦長naca0012だと27万くらいかな
		% root_hchord=(Sh-(lim_hchord*rs_h/1000*trn_h))/(rs_h/1000*trn_h+rs_h/1000*(brn_h-1));

		% %テーパー比決定
		% ut=(root_hchord-lim_hchord)/((trn_h)*rs_h)*1000;
		% ot=-ut;
		% drv=rs_h/10000*5;
		% for i=1:round((brn_h+trn_h+trn_h-1)*rs_h/drv/1000)+1
		% 	h_wing(i,1)=drv*(i-1);
		% 	if(h_wing(i,1)<(trn_h)*rs_h/1000)
		% 		h_wing(i,2)=ut*h_wing(i,1)+lim_hchord;
		% 	elseif(h_wing(i,1)>=(trn_h-1)*rs_h/1000&&h_wing(i,1)<=(brn_h+trn_h-1)*rs_h/1000)
		% 		h_wing(i,2)=root_hchord;
		% 	else
		% 		h_wing(i,2)=ot*(h_wing(i,1)-(trn_h+brn_h-1)*rs_h/1000)+root_hchord;
		% 	end
		% end

		% Sh=trapz(h_wing(:,1),h_wing(:,2))

	%矩形水平尾翼を作る
	%尾翼スパンから拘束
	rs_h=input("リブ間(m)=");
	rn_h=input("セクション数(奇数)=");
	span_h=rs_h*rn_h
	chord_h=Sh/span_h

	%水平尾翼を50分割
	for i=1:51
		h_wing(i,1)=span_h/50*(i-1);
		h_wing(i,2)=chord_h;
	endfor

	%全機揚力係数算出
	CL=L_lst*9.8*2/a_density/d_vel/mw_s;


	%尾翼揚力傾斜算出
	% dClh=tLST(h_wing,d_vel,5,nu,ahdata)(:,1);
	% alphawt=mean(dClh./5.*180./pi);

	%alphawt=6.818/(1+6.818/(pi*(span_h/chord_h)));
	alphawt=4.136755

	%尾翼効率
	etat=(d_vel^2+(2*mean(lst_Vn_i(1:we_div/5)))^2)/(d_vel^2+mean(lst_Vn_i(1:we_div/5))^2);

	%全機揚力傾斜
	CLalpha=alphaw+etat*Sh/mw_s*alphawt*(1-dep_dal);

	%全機空力中心
	xnp_cbar=0.25+etat*Vh*alphawt/CLalpha*(1-dep_dal);

	% %中央翼前縁からCbar位置前縁までのｘ方向長さ
	% cbar_ledge=f_data(1,2)*sparpos-cbar*spline(f_data(:,1),sparpos(:,1),cbar_pos);

	% %cbar位置における重心位置
	% cbar_centG=spline(f_data(:,1),sparpos(:,1),cbar_pos);

	%動圧
	qbar=0.5*a_density*d_vel^2;

	%飛行機効率
	efficiency=(L_lst^2/(1/2*a_density*d_vel^2*mw_s)^2)/(Di_lst/(1/2*a_density*d_vel^2*mw_s))/pi/aspect_ratio;

	%無次元化用動圧
	Kp=0.5*a_density*d_vel^2*mw_s;

	Cxu=-2/a_density/d_vel/mw_s*(D_lst+6)/d_vel-2*(Dp_lst/(1/2*a_density*d_vel^2*mw_s)+L_lst/(1/2*a_density*d_vel^2*mw_s)*tan(0));
	Cxa=CL*(1-2*CLalpha/pi/efficiency/aspect_ratio);
	Czu=0;
	Cza=-CLalpha;
	Cmu=0;
	Cma=CLalpha*(cog-xnp_cbar)
	Cmq=-2*Vh*(lhd/cbar)*alphawt;
	Cmadot=-2*Vh*(lhd/cbar)*alphawt*dep_dal;
	Czde=-Sh/mw_s*alphawt;
	Cmde=-Vh*alphawt;

	y_dcoefh=a_density*d_vel^2*mw_s*cbar/2/Ine(2,2);

	theta0=0;
	alpha0=0;

	Xa=Kp/L_lst*(Cxa+2*CL*tan(theta0)*tan(alpha0))
	Zu=Kp/d_vel/L_lst*(Czu-2*CL)
	Za=Kp/L_lst*(Cza-2*CL*tan(alpha0))
	Zde=Kp/L_lst*Czde
	Ma=y_dcoefh*cbar*Cma
	Mq=y_dcoefh/d_vel*cbar/2*Cmq
	Madot=y_dcoefh/d_vel*cbar/2*Cmadot
	Mde=y_dcoefh*cbar*Cmde

	omegasp=sqrt(-Ma+(Za/d_vel)*Mq);
	zetasp=(-Za/d_vel-Mq-Madot)/2/omegasp;
	base=(Ma/Za-Mq/d_vel);
	omega2na=L_lst*cbar/Ine(2,2)*(xnp_cbar-cog);

	base
	Sh
	Vh
	kfr=Vh*lhd/cbar
	omegasp
	zetasp
	omega2na

	sign=input("1:ok 0:return ---");

until(sign==1)