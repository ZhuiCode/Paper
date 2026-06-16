#!/bin/bash
# 下载文档中的图片
urls=(
  "https://cdn.nlark.com/yuque/0/2026/png/34609537/1780361074924-95cc8d02-952e-4c4e-ace2-8bbe598b20d7.png"
  "https://cdn.nlark.com/yuque/0/2026/png/34609537/1780381453603-d8d9d708-146b-432a-8445-1c321fc3d147.png"
  "https://cdn.nlark.com/yuque/0/2026/png/34609537/1780883752567-b8b39e64-05b6-4092-b706-f5951ff9059d.png"
  "https://cdn.nlark.com/yuque/0/2026/png/34609537/1780907915802-77c67cec-d72b-40c1-91a4-0b52b3ec2caf.png"
  "https://cdn.nlark.com/yuque/0/2026/png/34609537/1780883941268-f05977cf-56a6-4b6e-b937-cc9daa27c037.png"
)

names=("doc_figure_1.png" "doc_figure_2.png" "doc_figure_3.png" "doc_figure_4.png" "doc_figure_5.png")

for i in "${!urls[@]}"; do
  echo "下载 ${names[$i]}..."
  if command -v curl &>/dev/null; then
    curl -L -o "${names[$i]}" "${urls[$i]}" 2>/dev/null && echo "✓ 已下载 ${names[$i]}" || echo "✗ 下载失败 ${names[$i]}"
  elif command -v wget &>/dev/null; then
    wget -q -O "${names[$i]}" "${urls[$i]}" && echo "✓ 已下载 ${names[$i]}" || echo "✗ 下载失败 ${names[$i]}"
  fi
done
