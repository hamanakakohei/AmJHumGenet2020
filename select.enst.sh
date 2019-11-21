CODING="snpEff_genes.txt"
CANONICAL="snpeff.canonical.txt"
GENCODE="gencodev19.comprehensive.40bp.seq.txt"

awk '$1 !~ /^#/ && $4=="protein_coding"{print $3}' $CODING|sort -u > enst.coding.txt
awk '!a[$3]++{print $3}' $CANONICAL > enst.canonical.txt
paste <(awk -F"[_ ]" '$1 ~ /^>/{print $3}' $GENCODE) <(awk -F"[=: ]" '$1 ~ /^>/{print $3}' $GENCODE)|awk '$2 !~ /^_/{print $1}'|cut -d"." -f1|sort -u > enst.gencode.noalthaplo.txt
join <(sort enst.coding.txt) <(sort enst.canonical.txt)|join - <(sort enst.gencode.noalthaplo.txt) > enst.coding.canonical.gencode.noalthaplo.txt
