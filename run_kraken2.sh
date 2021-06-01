KRAKEN="kraken2/kraken2"
KRAKEN_DB="kraken2/kraken_db"
NUM_THREADS=8
DB=""
FILE=""

while getopts t:f:d: flag
do
    case "${flag}" in
        t) NUM_THREADS=${OPTARG};;
        d) DB=${OPTARG};;
    esac
done

DATA_DIR=datasets/$DB

$KRAKEN --db $KRAKEN_DB --report $DATA_DIR/$DB.report --output _ --threads $NUM_THREADS $DATA_DIR/*.fastq
