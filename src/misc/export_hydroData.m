function export_hydroData(DIR,hydroParameters,nBody,DoF)

if nargin<4
	DoF=[1 2 3]; nDoF=length(DoF);
end

W=hydroParameters.Omega; nw=length(W);
A=hydroParameters.A;
Ainf=hydroParameters.A_inf;
B=hydroParameters.B;
Fe=hydroParameters.Fe;
Fe_Haskind=hydroParameters.Fe_Haskind;

waterDepth=150;

if 1
%% EXPORT des efforts d'excitation
t='%.3f';
line=['#-w--'];
line2=['w'];
for i=1:length(DoF)
	t=[t,[' %.5e %.5e']];
	line=[line,[['----|Fe(',num2str(DoF(i)),')|'],['---arg(Fe(',num2str(DoF(i)),'))']]];
	line2=[line2,[[' Fe',num2str(DoF(i))],[' phi',num2str(DoF(i))]]];
end
t=[t,'\n']; line=[line,'\n']; line2=[line2,'\n'];

for k=1:nBody
	vFe=[]; vFe_Haskind=[];
	for i=1:length(DoF)
		vFe=[vFe,[abs(Fe(:,DoF(i)+(k-1)*3)),angle(Fe(:,DoF(i)+(k-1)*3))*(180/pi)]];
		vFe_Haskind=[vFe_Haskind,[abs(Fe_Haskind(:,DoF(i)+(k-1)*3)),angle(Fe_Haskind(:,DoF(i)+(k-1)*3))*(180/pi)]];
	end

	fid=fopen([DIR,filesep,'Fe',num2str(k),'.dat'],'w');
	% Écriture entête
	fprintf(fid,'# ---- Environment ---------------------------------------------------------- \n');
	fprintf(fid,'# 1025.0	! WATER DENSITY RHO 			! KG/M**3 	! Fluid specific volume \n');
	fprintf(fid,'# 9.81	! GRAVITY ACCELERATION G			! M/S**2	! Gravity \n');
	fprintf(fid,'# %.1f	! WATER DEPTH			! M		! Water depth\n',waterDepth);
	fprintf(fid,'# 0.	0.	! XEFF YEFF		! M		! Wave measurement point\n');
	fprintf(fid,'# ---- Data description ----------------------------------------------------- \n');
	fprintf(fid,'# OMEGA,DoF(1)->AMP/PHASE,...\n');
	fprintf(fid,'# UNIT\n');
	fprintf(fid,'#	* OMEGA: [rad/s]\n');
	fprintf(fid,'#	* AMP: [N/m]\n');
	fprintf(fid,'#	* PHASE: [degre]\n');
	fprintf(fid,line);
	fprintf(fid,line2);
	% écriture dans le fichier
	for w=1:nw
		fprintf(fid,t,[W(w),vFe(w,:)]);
	end
	fclose('all');

	fid=fopen([DIR,filesep,'Fe',num2str(k),'_Haskind.dat'],'w');
	% Écriture entête
	fprintf(fid,'# ---- Environment ---------------------------------------------------------- \n');
	fprintf(fid,'# 1025.0	! WATER DENSITY RHO 			! KG/M**3 	! Fluid specific volume \n');
	fprintf(fid,'# 9.81	! GRAVITY ACCELERATION G			! M/S**2	! Gravity \n');
	fprintf(fid,'# %.1f	! WATER DEPTH			! M		! Water depth\n',waterDepth);
	fprintf(fid,'# 0.	0.	! XEFF YEFF		! M		! Wave measurement point\n');
	fprintf(fid,'# ---- Data description ----------------------------------------------------- \n');
	fprintf(fid,'# OMEGA,DoF(1)->AMP/PHASE,...\n');
	fprintf(fid,'# UNIT\n');
	fprintf(fid,'#	* OMEGA: [rad/s]\n');
	fprintf(fid,'#	* AMP: [N/m]\n');
	fprintf(fid,'#	* PHASE: [degre]\n');
	fprintf(fid,line);
	fprintf(fid,line2);
	% écriture dans le fichier
	for w=1:nw
		fprintf(fid,t,[W(w),vFe_Haskind(w,:)]);
	end
	fclose('all');
	
end
end

if 1
%% EXPORT des masses ajoutées et des amortissements de radiation
for i=1:nBody
	for j=1:nBody
		fidA=fopen([DIR,filesep,'A',num2str(i),num2str(j),'.dat'],'w');
		fprintf(fidA,'w A11 A13 A15 A31 A33 A35 A51 A53 A55\n');
		fidB=fopen([DIR,filesep,'B',num2str(i),num2str(j),'.dat'],'w');
		fprintf(fidB,'w B11 B13 B15 B31 B33 B35 B51 B53 B55\n');
		for w=1:nw
			fprintf(fidA,'%.3f',W(w));
			fprintf(fidB,'%.3f',W(w));
			for k=1:nDoF
				for l=1:nDoF
					fprintf(fidA,' %.3e',A(3*(i-1)+DoF(k),3*(j-1)+DoF(l),w));
					fprintf(fidB,' %.3e',B(3*(i-1)+DoF(k),3*(j-1)+DoF(l),w));
				end
			end
			fprintf(fidA,'\n'); fprintf(fidB,'\n');
		end
		% Add one more line for infinite added mass
		fprintf(fidA,'\n');
		fprintf(fidA,'%.3f',Inf);
		for k=1:nDoF
			for l=1:nDoF
				fprintf(fidA,' %.3e',Ainf(3*(i-1)+DoF(k),3*(j-1)+DoF(l)));
			end
		end
		fclose('all');
	end
end

end
%% EXPORT des réponses impulsionelles IRF
if 1
IRF=hydroParameters.IRF.irf;
t=hydroParameters.IRF.t; nt=length(t);

nIRF=size(IRF,1);
vIRF=[t];
line='t';
format='%.2f';
for i=1:nBody
for j=1:nBody
	for p=1:nDoF
	for q=1:nDoF
		format=[format,' %.4e'];
		line=[line,[' IRF',num2str(3*(i-1)+DoF(p)),num2str(3*(j-1)+DoF(q))]];
		irf=reshape(IRF(3*(i-1)+DoF(p),3*(j-1)+DoF(q),:),nt,1);
		vIRF=[vIRF,irf];
	end
	end
% 		irf=reshape(IRF(i,j,:),nt,1);
% 		vIRF=[vIRF,irf];
end
end
format=[format,'\n']; line=[line,'\n'];

fid=fopen([DIR,filesep,'IRF.dat'],'w');
fprintf(fid,line);
for i=1:nt
	fprintf(fid,format,vIRF(i,:));
end
fclose('all');
end	
		
% 		
% 		fprintf(fidA,'w A11 A13 A15 A31 A33 A35 A51 A53 A55\n');
% 		fidB=fopen([DIR,filesep,'B',num2str(i),num2str(j),'.dat'],'w');
% 		fprintf(fidB,'w B11 B13 B15 B31 B33 B35 B51 B53 B55\n');
% 		for w=1:nw
% 			fprintf(fidA,'%.3f',W(w));
% 			fprintf(fidB,'%.3f',W(w));
% 			for k=1:nDoF
% 				for l=1:nDoF
% 					fprintf(fidA,' %.3e',A(3*(i-1)+DoF(k),3*(j-1)+DoF(l),w));
% 					fprintf(fidB,' %.3e',B(3*(i-1)+DoF(k),3*(j-1)+DoF(l),w));
% 				end
% 			end
% 			fprintf(fidA,'\n'); fprintf(fidB,'\n');
% 		end
% 		% Add one more line for infinite added mass
% 		fprintf(fidA,'\n');
% 		fprintf(fidA,'%.3f',Inf);
% 		for k=1:nDoF
% 			for l=1:nDoF
% 				fprintf(fidA,' %.3e',Ainf(3*(i-1)+DoF(k),3*(j-1)+DoF(l)));
% 			end
% 		end
% 
% 		fclose('all');

end

