STC = ((accSTAStims-stackSTA)' * (accSTAStims-stackSTA)) ./(totalAPs-1);

accSTC = zeros(30,30);
for i=1:totalAPs
    accSTC = accSTC + (accSTAStims(i,:)-STA)'*(accSTAStims(i,:)-STA);
end

STC2 = accSTC./(totalAPs-1);

