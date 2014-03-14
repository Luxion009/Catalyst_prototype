%使用するマンドレルを定義

clear mand_data w_leng dw_ep w_gap w_thn dw_len;

% ml=fopen('mandrellist.txt');
% mandrel_list=(fscanf(ml,'%f',[2,Inf]))';#２行無限列で読み込んだ後転置

do
% 	%翼の分割数を定義
% 	% w_div_c=input("半翼の分割数=");

% 	%翼の等分割時の位置、各パラメータ変化点を毎回表示しマンドレルを定義する
% 	for i=1:w_div_c
		
% 		do
% 			do
% 				printf("%d翼定義\n",i);
				
% 				printf("無変化矩形部位置 0-%f [m]\n",f_data(nc_p,1));
% 				printf("第一テーパー部位置 %f-%f [m]\n",f_data(nc_p,1),f_data(t1,1));
% 				printf("第二テーパー部位置 %f-%f [m]\n",f_data(t1,1),f_data(t2,1));
% 				printf("途中矩形部位置 %f-%f [m]\n",f_data(t2,1),f_data(m1,1));
% 				printf("第三テーパー部位置 %f-%f [m]\n",f_data(m1,1),f_data(t3,1));
% 				printf("第四テーパー部位置 %f-%f [m]\n",f_data(t3,1),f_data(t4,1));
				
% 				w_leng(i,1)=input("部分翼端位置(m)=");
% 				dw_ep(i,1)=lookup(f_data(:,1),w_leng(i,1))
				
% 				disp("マンドレルからプランク上面までの余裕(mm)=");
% 				w_gap(i,1)=input("w_gap=");
% 				if(i==1)
% 					w_thn(i,1)=min(foil_thn(1:dw_ep(i,1),1))*1000-w_gap(i,1)*2;%必要厚さ
% 					dw_len(i,1)=f_data(dw_ep(i,1),1);%部分翼長さ
% 				endif
% 				if(i!=1)				
% 					w_thn(i,1)=min(foil_thn(dw_ep(i-1,1)+1:dw_ep(i,1),1))*1000-w_gap(i,1)*2;
% 					dw_len(i,1)=f_data(dw_ep(i,1),1)-f_data(dw_ep(i-1,1),1);
% 				endif
				
% 				%接合部長さを決定
% 				mand_data(i,3)=input("接合部長さ(mm)=");
% 				if(i==1)
% 					junc=0;
% 				endif
				
% 				%必要長さの決定
% 				% mand_data_r(i,1)=dw_len(i,1)*1000+300+junc;
% 				% mand_data_r(i,2)=w_thn(i,1);
% 				% mand_data_r(i,3)=junc;
				
% 				mand_data(i,1)=dw_len(i,1);%部分翼長さ
% 				mand_data(i,2)=w_thn(i,1)+w_gap(i,1)*2;%最小翼厚
				
% 				printf("最小翼厚=%f(mm)\n",w_thn(i,1)+w_gap(i,1)*2);
% 				printf("想定限界桁径=%f(mm)\n",w_thn(i,1));
% 				printf("必要長さ=%f(mm)\n",mand_data(i,1)*1000+mand_data(i,3)+300);
% 				printf("接合部長さ=%f(mm)\n",mand_data(i,3));
				
% 				sign=input("ok:1 return:0 ---");
				
% 			until(sign==1)
			
			
% 			md(i,1)=lookup(mandrel_list(:,1),(mand_data(i,2)-w_gap(i,1)*2))
			
% 			%最大径以下で必要長さを満たすマンドレルを探索
% 			do
% 				if(mandrel_list(md(i,1),2)<mand_data(i,1)*1000+mand_data(i,3)+300)
% 					md(i,1)=md(i,1)-1;
% 				endif
% 			until(mandrel_list(md(i,1),2)>=mand_data(i,1)*1000+mand_data(i,3)+300||md(i,1)==1)
			
% 			if(md(i,1)>=1)
% 				mand_data(i,4)=mandrel_list(md(i,1),1);
% 				mand_data(i,5)=mandrel_list(md(i,1),2);
% 				rem=1;
% 			endif

% 			if(md(i,1)<=1)
% 				disp("該当マンドレル無し");
% 				disp("------------------------------------------------");
% 				rem=0;
% 			endif
			
% 			if(rem==1)
% 				disp("マンドレル確定");
% 				printf("径=%d(mm)		長さ=%d(mm)\n",mand_data(i,4),mand_data(i,5));
% 				rem=input("ok:1 return:0 ---");
% 				disp("------------------------------------------------")
% 			endif	
			
% 			if(rem!=1)
% 				md(i,:)=[];
% 				mand_data(i,:)=[];
% 			endif
		
% 		until(rem==1)
		
% 	endfor
	
	for i=1:w_div_c
		printf("%d翼\n",i);
		mand_data(i,4)=input("マンドレル径(mm)=");
		mand_data(i,5)=input("マンドレル長さ(mm)=");
		mand_data(i,1)=divw_span(i,1);
		mand_data(i,3)=input("接合部長さ(mm)=");
		mand_data(i,2)=real_thn_xcp(i,1);
	endfor

	disp("確定マンドレルデータ");
	printf("------------------------------------------------\n");
	for i=1:w_div_c
		printf("%d翼\n",i);
		% printf("必要径=%f(mm)	必要長さ=%f(mm)\n",mand_data(i,2),mand_data(i,1));
		printf("径=%f(mm)		長さ=%f(mm)\n",mand_data(i,4),mand_data(i,5));
	endfor
	printf("------------------------------------------------\n");

	sign=input("ok:1 return:0 ---");
			
until(sign==1)