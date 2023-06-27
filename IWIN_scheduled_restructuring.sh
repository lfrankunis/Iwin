cd /mnt/c/Github/Iwin/
python ./IWIN_scheduled_restructuring.py
for s in Narveneset Bohemanneset Daudmannsodden Gasoyane
do
	for t in 1min 10min
	do
		rsync -avh -t -O --delete /mnt/c/Data/sorted_by_location/lighthouse_AWS_${s}/${t} unismet@157.249.75.166:/data/unismet/iwin/Data_transfer_MET/lighthouse_AWS_${s}/
	done
done
for s in MSBerg MSPolargirl MSBillefjord MSBard
do
	for t in 20sec 1min 10min
	do
		rsync -avh -t -O --delete /mnt/c/Data/sorted_by_location/mobile_AWS_${s}/${t} unismet@157.249.75.166:/data/unismet/iwin/Data_transfer_MET/mobile_AWS_${s}/
	done
done
