# script for downloading and compressing Seismic datase
curl http://snap.stanford.edu/seismic/data.csv | gzip > data.csv.gz
curl http://snap.stanford.edu/seismic/index.csv | gzip > index.csv.gz
