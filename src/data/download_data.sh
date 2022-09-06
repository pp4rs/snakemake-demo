set -e

cd data/raw
wget --no-check-certificate -O QOB.rar https://economics.mit.edu/files/2853
unrar x QOB.rar
mv QOB QOB.txt
rm QOB.rar
