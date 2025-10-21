#测试std数据
g++ -std=c++20 source.cpp -o source

t=10
maxround=20
while (($t <= $maxround))
do
    in="$t.in"
    out="$t.myout"
    ans="$t.out"
    ./source < "$in" > "$out"    
    if diff -Z "$out" "$ans"; then
        echo "AC on test $t"
    else 
        echo "WA on test $t"
        #break
    fi
    ((t++))
done