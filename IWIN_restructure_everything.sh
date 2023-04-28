ENDDAY=$(date -d "-2 day" "+ %Y%m%d")
#python ./IWIN_manual_restructuring.py -s 1884 -t 20220617 "$ENDDAY"
#python ./IWIN_manual_restructuring.py -s 1885 -t 20210819 "$ENDDAY"
#python ./IWIN_manual_restructuring.py -s 1886 -t 20220708 "$ENDDAY"
#python ./IWIN_manual_restructuring.py -s 1887 -t 20220903 "$ENDDAY"

python ./IWIN_manual_restructuring.py -s 1883 -t 20210505 20211021
python ./IWIN_manual_restructuring.py -s 1872 -t 20210528 20211020

python ./IWIN_manual_restructuring.py -s 1883 -t 20220530 20221027
python ./IWIN_manual_restructuring.py -s 1872 -t 20220531 20221017
python ./IWIN_manual_restructuring.py -s 1924 -t 20220321 20221027

python ./IWIN_manual_restructuring.py -s 1883 -t 20230322 "$ENDDAY"
