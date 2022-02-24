function stableStdMean = findstable(cells)
ref = cells{1};
setOfNearest = zeros(length(ref), length(cells));

setOfNearest(:, 1) = ref;

for i = 2:length(cells)
    for j = 1:length(ref)
        setOfNearest(j, i) = findNearest(ref(j), cells{i});
    end
end

m = mean(setOfNearest, 2);
s = std(setOfNearest, 0, 2);

[s, I] = sort(s, 'ascend');
stableStdMean = [s m(I)];
end

function nearest = findNearest(ref, set)
if isempty(set)
    nearest = 90;
else
    nearest = set(1);
    for i = 2:length(set)
        if abs(nearest - ref) > abs(set(i) - ref)
            nearest = set(i);
        end
    end
end
end

