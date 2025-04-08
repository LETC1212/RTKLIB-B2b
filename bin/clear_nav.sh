#!/bin/bash

# 定义文件路径
FILE="rtkrcv.nav"

# 检查文件是否存在
if [ -f "$FILE" ]; then
    # 清空文件内容
    > "$FILE"
    echo "文件内容已清空: $FILE"
else
    echo "错误：文件 $FILE 不存在。"
fi
