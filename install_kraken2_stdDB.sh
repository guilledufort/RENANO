cd kraken2
curl ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v2_8GB_201904.tgz -o kraken_db.tgz
tar zxvf kraken_db.tgz
rm kraken_db.tgz
mv minikraken2_v2_8GB_201904_UPDATE kraken_db