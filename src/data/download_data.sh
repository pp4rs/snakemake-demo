set -e

URL=$1
FILEPATH=$2

DIR="$(dirname "${FILEPATH}")"
FILE="$(basename "${FILEPATH}")"

cd "$DIR"
wget --no-check-certificate -O "QOB.rar" "$URL"
unrar x QOB.rar
mv QOB "$FILE"
rm QOB.rar
