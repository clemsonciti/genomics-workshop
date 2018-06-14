module load gnu-parallel
ls dracula.txt.* | parallel -j15 "wordfreq {} > count.{}"
