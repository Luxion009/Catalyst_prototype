%�׏d�����p�̂Ȃɂ�

%���f�����ʒu��`

w_div=dlmread('weight_div5.csv');

w_div=w_div/1000;
%���ʒu
w_pos(1,1)=w_div(1,1)/2;

for i=2:rows(w_div)
	w_pos(i,1)=(w_div(i,1)-w_div(i-1,1))/2+w_div(i-1,1);
end

%-----���f�ɔ�������cl�ƌ��d�ʂ��v�Z����B
clear yst
yst=linspace(0,w_div(1,1),20)';

Lst=spline(S_mat(:,1),S_mat_dL(:,1),yst);
wst=interp1(S_mat(:,1),S_mat(:,6),yst,'linear');

wst(1,1)=wst(2,1);

wst=wst./1000;
%1G
weight(1,1)=trapz(yst(:,1),(Lst(:,1)./9.8-wst(:,1)));%N
%1.25G
weight(1,2)=trapz(yst(:,1),(1.25.*Lst(:,1)./9.8-1.25.*wst(:,1)));

for i=2:rows(w_div)
		clear yst
		yst=linspace(w_div(i-1,1),w_div(i,1),20)';

		Lst=interp1(S_mat(:,1),S_mat_dL(:,1),yst);
		wst=interp1(S_mat(:,1),S_mat(:,6),yst,'linear');

		wst(1,1)=wst(2,1);
		wst=wst./1000;
	%-----�ϕ�
		%1G
		weight(i,1)=trapz(yst(:,1),(Lst(:,1)./9.8-wst(:,1)));%N
		%1.25G
		weight(i,2)=trapz(yst(:,1),(1.25.*Lst(:,1)./9.8-1.25.*wst(:,1)));

endfor

%20�Ԗڂ̂�����ɂQ�ŊO��WL�̗g�͂�����
% org_w(20,1:2)=weight(20,1:2);

% weight(20,1)=weight(20,1)+trapz(S_mat(401:500,1),S_mat_dL(401:500,1))./9.8+0.31666/9.8;
% weight(20,2)=weight(20,2)+1.25.*trapz(S_mat(401:500,1),S_mat_dL(401:500,1))./9.8+1.25.*0.31666/9.8;

weight(20,1)=weight(20,1)+0.31666/9.8;
weight(20,2)=weight(20,2)+1.25.*0.31666/9.8;


%-----�d���U�蕪����
for i=1:rows(weight)
	if mod(i,2)==0
		weight13(i,1)=sum(weight(i-1:i,2)-weight(i-1:i,1));
	else
		weight13(i,1)=0;
	end
end
[w_pos,weight(:,1),weight13]



%���������ꂽ����
for gst=1:2

	disp("generating Q");
	for i=1:div_c*dp_num-1
		S_mat_dQ(i,gst)=(S_mat(i,6).*9.8./1000+S_mat(i+1,6).*9.8./1000)./2*(S_mat(i+1,1)-S_mat(i,1));
	end
	% S_mat_dQ(div_c*dp_num,gst)=interp1(S_mat(1:div_c*dp_num-1,gst),S_mat_dQ(:,gst),S_mat(div_c*dp_num,1),'linear','extrap');
	S_mat_dQ(div_c*dp_num,gst)=(S_mat_dQ(div_c*dp_num-1,gst)-S_mat_dQ(div_c*dp_num-2,gst))/(S_mat(div_c*dp_num-1,gst)-S_mat(div_c*dp_num-2,gst))*(S_mat(div_c*dp_num,gst)-S_mat(div_c*dp_num-1,gst))+S_mat_dQ(div_c*dp_num-1,gst);

	%�����蕪���Ǐ�����f�͂ɉ��Z����
	for j=1:rows(weight)
		for i=2:rows(S_mat(:,1))
			if (w_pos(j,1)<=S_mat(i,1))
				if (w_pos(j,1)>=S_mat(i-1,1))
					if (gst==1)
						S_mat_dQ(i,gst)=S_mat_dQ(i,gst)+weight(j,1).*9.8;
					elseif(gst==2)
						S_mat_dQ(i,gst)=S_mat_dQ(i,gst)+weight13(j,1).*9.8+weight(j,1).*9.8;
					end
						
				end
			end
		end
	end


	S_mat_Q(div_c*dp_num+1,gst)=0;%���݂��Ȃ��_�����v�Z�̂��߂ɒ�`

	for i=div_c*dp_num:-1:1
		S_mat_Q(i,gst)=S_mat_Q(i+1,gst)+S_mat_dQ(i,gst);
	end
	%-----����f��(kgf)
	S_mat_Q=S_mat_Q./9.8;

	disp("generating M");
	for i=1:div_c*dp_num-1
		S_mat_dM(i,gst)=(S_mat_Q(i,gst)+S_mat_Q(i+1,gst))/2*(S_mat(i+1,1)-S_mat(i,1))*1000;
	end
	%-----���[�����g(kgf.m)
	S_mat_M(div_c*dp_num,gst)=0;
	for i=div_c*dp_num-1:-1:1
		S_mat_M(i,gst)=S_mat_M(i+1,gst)+S_mat_dM(i,gst);
	end

	disp("generating Zz");
	for n=1:div_c
		Zz(dp_num*(n-1)+1:dp_num*n,1)=S_I(dp_num*(n-1)+1:dp_num*n,2)./ply_data{n}(rows(ply_data{n}),1).*1000;%�f�ʌW��
		S_max_Z(1:dp_num*n,gst)=S_mat_M(1:dp_num*n,gst)./Zz(1:dp_num*n,1).*9.8./1000;%�ő�Ȃ�����
	end
			
	for i=1:div_c*dp_num-1
		S_mat_dtheta(i,gst)=(S_mat_M(i,gst)/S_H(i,2)+S_mat_M(i+1,gst)/S_H(i+1,2))/2*(S_mat(i+1,1)-S_mat(i,1));%�A�������̋ȗ��̕��ςɂ��̊Ԃ̋�����������
		%�ڐ��̊p�x�ω�
	end
	S_mat_theta(1,gst)=0;%�ڐ��Ɛ������̊p�x
	for i=2:div_c*dp_num-1
		S_mat_theta(i,gst)=S_mat_theta(i-1,gst)+S_mat_dtheta(i,gst);
	end
	% S_mat_theta(div_c*dp_num,gst)=interp1(S_mat(1:div_c*dp_num-1,gst),S_mat_theta(:,gst),S_mat(div_c*dp_num,1),'linear','extrap');
	S_mat_theta(div_c*dp_num,gst)=(S_mat_theta(div_c*dp_num-1,gst)-S_mat_theta(div_c*dp_num-2,gst))/(S_mat(div_c*dp_num-1,gst)-S_mat(div_c*dp_num-2,gst))*(S_mat(div_c*dp_num,gst)-S_mat(div_c*dp_num-1,gst))+S_mat_theta(div_c*dp_num-1,gst);



	for i=1:div_c*dp_num-1
		S_mat_dswang(i,gst)=(S_mat_theta(i,gst)+S_mat_theta(i+1,gst))./2.*(S_mat(i+1,1)-S_mat(i,1));%sin x=x�Ƌߎ��ł���̂�Z�����ω����Z�o�ł���
	end
	S_mat_swang(1,gst)=S_mat_dswang(1,gst);
	for i=2:div_c*dp_num-1
		S_mat_swang(i,gst)=S_mat_swang(i-1,gst)+S_mat_dswang(i,gst);%Z�����ʒu
	end
	% S_mat_swang(div_c*dp_num,gst)=interp1(S_mat(1:div_c*dp_num-1,gst),S_mat_swang(:,gst),S_mat(div_c*dp_num,1),'linear','extrap');
	S_mat_swang(div_c*dp_num,gst)=S_mat_swang(div_c*dp_num-1,gst);
	S_mat_swang(div_c*dp_num,gst)=(S_mat_swang(div_c*dp_num-1,gst)-S_mat_swang(div_c*dp_num-2,gst))/(S_mat(div_c*dp_num-1,gst)-S_mat(div_c*dp_num-2,gst))*(S_mat(div_c*dp_num,gst)-S_mat(div_c*dp_num-1,gst))+S_mat_swang(div_c*dp_num-1,gst);


	S_mat_swang(:,gst)=S_mat_swang(:,gst).*10^(-3);

endfor

disp("output")

% figure(124)
% 	% subplot(2,1,1)
% 		plot(S_mat(:,1),S_mat_swang(:,1),S_mat(:,1),S_mat_swang(:,2))
% 		grid on 
% 		axis('equal')
% 		title('deflectiton curve')

% 		ylabel('v mm')
% 		xlabel('y mm')
	% subplot(2,1,2)
	% 	plot(S_mat(:,1),S_mat_dQ(:,1),'marker','*')
	% 	hold on
	% 	plot(S_mat(:,1),S_mat_dQ(:,2),'marker','x','linestyle','none')
	% 	grid on
	% 	title('SFD')
	% 	ylabel('q kgf')
	% 	xlabel('y mm')

		% for i=1:div_c
		% 	sdiv(i,1)=S_mat(100*i,1);
		% end
		% plot(sdiv,1,'marker','+','linestyle','none')
		% legend('1G','1.3G','positions of spar division')
		% hold off
		

%�t�@�C���o��

dlmwrite('weight_position.csv',[w_pos,weight(:,1),weight13,w_pos-0,w_pos-1.61,w_pos-4.83,w_pos-8.05,w_pos-11.04], ',');
