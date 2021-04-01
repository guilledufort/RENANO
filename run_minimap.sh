MINIMAP="minimap2/minimap2"
MINI_PARAMS="-x map-ont -t 30 --secondary=no --cs"

DATA_DIR=datasets/$1
REF=$DATA_DIR/$1_genome.fna
RES_DIR=res_minimap

function mini {
    $MINIMAP $MINI_PARAMS $REF $1 > $DATA_DIR/$2.paf
}

foo () {
    filename="${1##*/}"
    filename="${filename%.*}"
    echo $filename
    echo $DATA_DIR/$filename.paf
    if [ ! -e "$DATA_DIR/$filename.paf" ];
    then
        mini $1 $filename
    else
        echo ".paf already exists..."
    fi

}
###############################################################################

#rm -f $TEMP_DIR/output.* $TEMP_DIR/temp.fastq

for fullfile in `ls -lS $DATA_DIR/*.fastq | awk '{ print $9; }'`
do
   foo $fullfile
done