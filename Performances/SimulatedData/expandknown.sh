var=$1 #number of repetitions to add

shuf KnownRepetitions.bed | head -1 > random.region.bed #choose random region

awk -v val=$var 'OFS=FS="\t"''{split($4,a,"x"); print $1,$2,$3, "tandem repeat expansion", a[2]":"val}' random.region.bed > h1.bed #true expansion
awk 'OFS=FS="\t"''{split($4,a,"x"); print $1,$2,$3, "tandem repeat expansion", a[2]":"0}' random.region.bed > h2.bed #fake contraction, no expansion
