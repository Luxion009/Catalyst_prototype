do
	ss2=input("‚PFC³@‚QF•”•ªÏ‘w‘‰Á@‚RFˆêüÏ‘w‘‰Á@‚SFÏ‘wŒ¸­ ---");

	if ss2==1
		rw=input("C³—ƒ”Ô†=");
		rs=input("‘w”Ô†=");
		rl=input("Ä’è‹`”ÍˆÍ 1:1-1 2:1-2 3:1-3 4:1-4 ---");
		ply_data{rw}=ply_leng(ply_data{rw},rs,mand_data(rw,1)/4*rl*1000);

	elseif ss2==2
		rw=input("‘‰Á—ƒ”Ô†=");
		ra=input("Ï‘wŠp“x=");
		rl=input("’è‹`”ÍˆÍ 1:1-1 2:1-2 3:1-3 4:1-4 ---");

		ply_data{rw}=ply_ex(ply_data{rw},1,ra,C_T,mand_data(rw,1)./4*rl*1000)

	elseif ss2==3
		rw=input("‘‰Á—ƒ”Ô†=");
		rs=input("‘w”Ô†=");
		rl=input("’è‹`”ÍˆÍ 1:1-1 2:1-2 3:1-3 4:1-4 ---");

		ply_data{rw}=ply_ins(ply_data{rw},rs,C_T,mand_data(rw,1)./4*rl*1000);

	elseif ss2==4
		rw=input("Œ¸­—ƒ”Ô†=");
		rs=input("‘w”Ô†=");
		ply_data{rw}(rs,:)=[];
		for j=2:rows(ply_data{rw})
			ply_data{rw}(j,1)=ply_data{rw}(j-1,1)+C_T(lookup(C_T(:,1),ply_data{rw}(j,2)),2)*2;
		end

	endif

	%Ï‘w\¬‚ğ•\¦
	for i=1:w_div_c
		ply_disp(ply_data{i},i);
	endfor

	bend_disp
	x_bend_disp

	figure(19);
	plot(S_mat(:,1),S_mat_swang(:,1),"-;swang;",S_mat(:,1),ID_swang(:,1),"-;id;");
	xlabel('y[m]');
	ylabel('z[m]');
	grid on;
	axis equal;
	% print("output/defspar.png","-dpng","-r100")


	figure(43)
	plot(S_mat(:,1),S_mat_Dswang(:,1),"-;swang;");
	xlabel('y[m]');
	ylabel('z[m]');
	grid on;

	figure(20);
	plot(S_mat(:,1),S_max_Z(:,1),"-");
	xlabel('y[m]');
	ylabel('s_max_z[Mpa]');
	grid on;
	% print("output/defspar.png","-dpng","-r100")

	figure(21)
	plot(S_mat(:,1),D_mat(:,1),"-");
	xlabel("y[m]");
	ylabel("line density");
	grid on;

	%Œ…‚ª’Ê‚é‚©ÅIŠm”F

	for i=1:w_div_c
		ps_gap(i,1)=(real_thn_xcp(i,1)*1000-ply_data{i}(rows(ply_data{i}),1))/2-p_thn;
		if(ps_gap(i,1)>=lim_gap)
			printf("[%d—ƒ]@Œ…-ƒvƒ‰ƒ“ƒNŠÔ—]—T : %f \n",i,ps_gap(i,1));
		else
			printf("[%d—ƒ]@ŒÀŠE—]—TˆÈ‰º@—ƒ’è‹`•s‰Â”\ : %f \n",i,ps_gap(i,1));
		end
	end

	W_Spar=W_Spar
	W_wing=W_wing

	ssign=input("next:1 end:0 ---")
		
	until(ssign==0)