#!/bin/bash
# make_datalist.sh - 生成数据集清单文件

# 定义输入文件匹配模式
input_dirs_regex=".*triOniaVtxValid.*202[2-4].*\/[0-9]\{4\}"

# 查找输入文件并生成清单
find /eos/user/c/chiw/JpsiJpsiUps/rootNtuple/ParkingDoubleMuonLowMass[0-7] \
     -type d -regextype 'posix-basic' -regex "$input_dirs_regex" > datadir.list

echo "Generated dataset list with $(wc -l < datadir.list) entries"