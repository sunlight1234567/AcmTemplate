#!/bin/bash
# 编译运行cpp

# AddressSanitizer (ASan): 检测内存越界、悬挂指针等问题。 启用方式：-fsanitize=address
# ThreadSanitizer (TSan): 检测多线程数据竞争。 启用方式：-fsanitize=thread
# UndefinedBehaviorSanitizer (UBSan): 检测未定义行为（如整数溢出）。 启用方式：-fsanitize=undefined
# LeakSanitizer (LSan): 检测内存泄漏。 启用方式：-fsanitize=leak
# -fno-omit-frame-pointer : 检测到内存错误时打印函数调用栈，这个参数一直都带上
# g++ -std=c++23 a.cpp -o a
g++ -fsanitize=address -g -std=c++17 a.cpp -o a

if [ $? -ne 0 ]; then
    echo "编译失败"
    exit 1
fi
# 运行并打印时间信息
# /usr/bin/time -f "real: %e\nuser: %U\nsys: %S" ./a < a.in > a.out
/usr/bin/time -f "time: %e s\nmem: %M KB" ./a < a.in > a.out
# ./a < a.in > a.out
# echo "Input:"
# cat a.in
# echo
# echo "Output:"
# cat a.out