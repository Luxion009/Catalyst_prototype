%�g�p����}���h�������`

clear mand_data w_leng dw_ep w_gap w_thn dw_len;

% ml=fopen('mandrellist.txt');
% mandrel_list=(fscanf(ml,'%f',[2,Inf]))';#�Q�s������œǂݍ��񂾌�]�u

do
% 	%���̕��������`
% 	% w_div_c=input("�����̕�����=");

% 	%���̓��������̈ʒu�A�e�p�����[�^�ω��_�𖈉�\�����}���h�������`����
% 	for i=1:w_div_c
		
% 		do
% 			do
% 				printf("%d����`\n",i);
				
% 				printf("���ω���`���ʒu 0-%f [m]\n",f_data(nc_p,1));
% 				printf("���e�[�p�[���ʒu %f-%f [m]\n",f_data(nc_p,1),f_data(t1,1));
% 				printf("���e�[�p�[���ʒu %f-%f [m]\n",f_data(t1,1),f_data(t2,1));
% 				printf("�r����`���ʒu %f-%f [m]\n",f_data(t2,1),f_data(m1,1));
% 				printf("��O�e�[�p�[���ʒu %f-%f [m]\n",f_data(m1,1),f_data(t3,1));
% 				printf("��l�e�[�p�[���ʒu %f-%f [m]\n",f_data(t3,1),f_data(t4,1));
				
% 				w_leng(i,1)=input("�������[�ʒu(m)=");
% 				dw_ep(i,1)=lookup(f_data(:,1),w_leng(i,1))
				
% 				disp("�}���h��������v�����N��ʂ܂ł̗]�T(mm)=");
% 				w_gap(i,1)=input("w_gap=");
% 				if(i==1)
% 					w_thn(i,1)=min(foil_thn(1:dw_ep(i,1),1))*1000-w_gap(i,1)*2;%�K�v����
% 					dw_len(i,1)=f_data(dw_ep(i,1),1);%����������
% 				endif
% 				if(i!=1)				
% 					w_thn(i,1)=min(foil_thn(dw_ep(i-1,1)+1:dw_ep(i,1),1))*1000-w_gap(i,1)*2;
% 					dw_len(i,1)=f_data(dw_ep(i,1),1)-f_data(dw_ep(i-1,1),1);
% 				endif
				
% 				%�ڍ�������������
% 				mand_data(i,3)=input("�ڍ�������(mm)=");
% 				if(i==1)
% 					junc=0;
% 				endif
				
% 				%�K�v�����̌���
% 				% mand_data_r(i,1)=dw_len(i,1)*1000+300+junc;
% 				% mand_data_r(i,2)=w_thn(i,1);
% 				% mand_data_r(i,3)=junc;
				
% 				mand_data(i,1)=dw_len(i,1);%����������
% 				mand_data(i,2)=w_thn(i,1)+w_gap(i,1)*2;%�ŏ�����
				
% 				printf("�ŏ�����=%f(mm)\n",w_thn(i,1)+w_gap(i,1)*2);
% 				printf("�z����E���a=%f(mm)\n",w_thn(i,1));
% 				printf("�K�v����=%f(mm)\n",mand_data(i,1)*1000+mand_data(i,3)+300);
% 				printf("�ڍ�������=%f(mm)\n",mand_data(i,3));
				
% 				sign=input("ok:1 return:0 ---");
				
% 			until(sign==1)
			
			
% 			md(i,1)=lookup(mandrel_list(:,1),(mand_data(i,2)-w_gap(i,1)*2))
			
% 			%�ő�a�ȉ��ŕK�v�����𖞂����}���h������T��
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
% 				disp("�Y���}���h��������");
% 				disp("------------------------------------------------");
% 				rem=0;
% 			endif
			
% 			if(rem==1)
% 				disp("�}���h�����m��");
% 				printf("�a=%d(mm)		����=%d(mm)\n",mand_data(i,4),mand_data(i,5));
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
		printf("%d��\n",i);
		mand_data(i,4)=input("�}���h�����a(mm)=");
		mand_data(i,5)=input("�}���h��������(mm)=");
		mand_data(i,1)=divw_span(i,1);
		mand_data(i,3)=input("�ڍ�������(mm)=");
		mand_data(i,2)=real_thn_xcp(i,1);
	endfor

	disp("�m��}���h�����f�[�^");
	printf("------------------------------------------------\n");
	for i=1:w_div_c
		printf("%d��\n",i);
		% printf("�K�v�a=%f(mm)	�K�v����=%f(mm)\n",mand_data(i,2),mand_data(i,1));
		printf("�a=%f(mm)		����=%f(mm)\n",mand_data(i,4),mand_data(i,5));
	endfor
	printf("------------------------------------------------\n");

	sign=input("ok:1 return:0 ---");
			
until(sign==1)