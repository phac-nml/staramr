# AMR Detection

## Preparing Databases

for i in *.fsa; do makeblastdb -in $i -parse_seqids -dbtype nucl; done
