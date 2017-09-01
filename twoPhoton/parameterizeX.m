function [ X ] = parameterizeX( licks )
%PARAMETERIZEX Summary of this function goes here
%   Detailed explanation goes here
    
    onsets = zeros(size(licks));
    offsets = zeros(size(licks));
    midbout = zeros(size(licks));
    licking = licks > 0;
    for j=1:size(licks,2)
        for k=1:size(licks,1)
            if j<22
                if j==1
                    if licking(k,j)
                        onsets(k,j) = 1;
                    end
                else
                    if licking(k,j) && sum(licking(k,1:j-1))==0
                        onsets(k,j) = 1;
                    end
                end 
            elseif j>size(licks,2)-21

                if licking(k,j) && sum(licking(k,j-21:j-1))==0
                    onsets(k,j) = 1;
                end

                if licking(k,j) && sum(licking(k,j+1:end))==0
                    offsets(k,j) = 1;
                end
            else
                if licking(k,j) && sum(licking(k,j-21:j-1))==0
                    onsets(k,j) = 1;
                end

                if licking(k,j) && sum(licking(k,j+1:j+21))==0
                    offsets(k,j) = 1;
                end
            end
        end
    end

    for k=1:size(licks,1)
        midbout(k,find(onsets(k,:)):find(offsets(k,:))) = 1;
    end

    nframes = size(licks,2);
    ntrials = size(licks,1);

    % Set feature matrix
    X = zeros(numel(licks),1+43+18+18+40);
    %               task coefficients-on-off-bout
    for k=1:ntrials
        tstart = ((k-1)*nframes);
        X(tstart+18:tstart+60,1:43) = eye(numel(tstart+18:tstart+60));
        onind = find(onsets(k,:));
        offind = find(offsets(k,:));
        ontoend = nframes-onind;
        for l=1:numel(onind)
            X(tstart+onind(l)-min(onind(l)-1,3):tstart+onind(l)+min(ontoend(l),14),46-min(onind(l)-1,3):43+min(ontoend(l)+3,17)) = eye(1+min(ontoend(l)+3,17)-3+min(onind(l)-1,3));
        end
        offtoend = nframes-offind;
        for l=1:numel(offind)
            X(tstart+offind(l)-min(offind(l)-1,3):tstart+offind(l)+min(offtoend(l),14),61:61+min(offtoend(l)+3,17)) = eye(1+min(offtoend(l)+3,17));
        end
        if (onind(1)+15 < offind(1)-3)
            boutind = onind(1)+15;
            bouttooff = (offind(1)-3)-boutind;
            X(tstart+boutind:tstart+boutind+min(bouttooff,39),79:79+min(bouttooff,39)) = eye(1+min(bouttooff,39));
        end
        if numel(offind)>2
            boutind = onind(2)+15;
            bouttooff = (offind(2)-3)-boutind;
            X(tstart+boutind:tstart+boutind+min(bouttooff,39),79:79+min(bouttooff,39)) = eye(1+min(bouttooff,39));
        end

    end
    X(:,end) = 1;
    assert(size(X,1)==numel(licks));

end

