echo "Count Start"

# rm Strange_Count/Strange_Count.root -f
mkdir Strange_Count
for i in {1..5};do bsub -q l "basf2 Strange-Count.py -i C00-Generation/D_s+/ccbarDs+EventGeneration_$i.root -o Strange_Count/Strange-Count_$i.root";done
sleep 5m # Waits 5 minutes.
bjobs
hadd -f Strange_Count.root Strange_Count/*.root
rm -r Strange_Count -f

echo "Count Complete"