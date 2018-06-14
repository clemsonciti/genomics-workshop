for file in dracula.txt.*
do
    wordfreq $file > count.$file &
done
wait
