function [irf,t]= IRFfcn(W,B,tmax,dt)

t=(0:dt:tmax)'; nt=length(t);
nw=length(W); dw=W(2)-W(1);


[m,n,p]=size(B);
if p==1
	B=reshape(B,1,1,nw);
end


irf=zeros(size(B,1),size(B,1));
for i=1:size(B,1)
	for j=1:size(B,1)
		for k=1:nt
		Bw=reshape(B(i,j,:),nw,1);
		irf(i,j,k)=(2/pi)*sum(Bw.*cos(W*t(k)))*dw;
		end
	end
end

end

