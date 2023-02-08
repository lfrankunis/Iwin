source /home/unismet/miniconda3/etc/profile.d/conda.sh
conda activate iwin
cd /mnt/c/Github/Iwin/
python ./IWIN_scheduled_restructuring.py
python ./IWIN_copy_files_for_transfer_MET.py 5
conda deactivate
rsync -avh --delete /mnt/c/Data_transfer_MET unismet@157.249.75.166:/data/unismet/iwin