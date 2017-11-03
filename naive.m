function ann = naive(A, B)
%%
patch_len = 7;
% A = double(uint8(A*255));
% B = double(B);
A = (double(A));
B = (double(B));

%%
[ma, na, ~] = size(A);
[mb, nb, ~] = size(B);

% A = padarray(A, [patch_len-1, patch_len-1], 'post');
% B = padarray(B, [patch_len-1, patch_len-1], 'post');

countA = 1;
for u = 1:ma
    for v = 1:na
        minx = min(u, ma-patch_len); miny = min(v,na-patch_len);
        patchA{u,v} = A(minx:minx+patch_len-1,miny:miny+patch_len-1,:);
    end
end
countB = 1;        
for u = 1:mb
    for v = 1:nb
        minx = min(u,mb-patch_len); miny = min(v,nb-patch_len);
        patchB{u, v} = B(minx:minx+patch_len-1,miny:miny+patch_len-1,:);
    end
end
%%
xcoordinate = zeros(ma,na);
ycoordinate = zeros(ma,na);
L2 = inf*ones(ma,na);
%{
% for ii = 1:numel(patchA)
%     for u = 1:size(patchB,1)
%         for v = 1:size(patchB,2)
%             cost = sum(sum((patchA{ii} - patchB{u,v}).^2));
%             if cost < L2(ii)
%                 xcoordinate(ii) = u;
%                 ycoordinate(ii) = v;
%                 L2(ii) = cost;
%             end
%         end
%     end
% end
%}
for ii = 1:numel(patchA)
    for jj = 1:numel(patchB)
        cost = sum(sum(sum((patchA{ii} - patchB{jj}).^2)));
        if cost < L2(ii)
            L2(ii) = cost;
            if mod(ii,ma) == 0
                xcoordinate(ii) = ma;
                ycoordinate(ii) = ii/ma;
            else
                xcoordinate(ii) = mod(ii,ma);
                ycoordinate(ii) = ceil(ii/ma);
            end
        end
    end
end
%%
ann(:,:,2) = xcoordinate;
ann(:,:,1) = ycoordinate;
ann(:,:,3) = L2;
%%
ann(:,end-(patch_len-2):end,1:2) = 0;
ann(:,end-(patch_len-2):end,3) = 2147483647;
ann(end-(patch_len-2):end,:,1:2) = 0;
ann(end-(patch_len-2):end,:,3) = 2147483647;


ann = int32(ann);
end








