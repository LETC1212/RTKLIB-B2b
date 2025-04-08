#!/bin/bash

# 确保脚本在出现任何错误时终止执行
set -e

# 打印正在删除的文件类型的信息
echo "Deleting all .trace, .log, .tag, and .pos files in the current directory..."

# 删除当前目录下的指定文件
rm -f ./*.trace
rm -f ./*.log
rm -f ./*.tag
rm -f ./*.pos
rm -f ./*.sp3
rm -f ./*.B2bssr
./clear_nav.sh

echo "Deletion complete."
