#!/bin/bash
INPUT_DIR=$1

# 生成处理脚本（无需处理输出路径）
sed -e "s|JOB_DATA|${INPUT_DIR}/mymultilep*.root|g" runReadTreeX.C.tpl \
  | sed -e 's|// #define RUN_JOB|#define RUN_JOB|g' > runReadTreeX.C

# 直接运行（输出文件自动由HTCondor传输到EOS）
root -l -b -q runReadTreeX.C