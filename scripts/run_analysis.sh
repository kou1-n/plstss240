#!/bin/bash
# PLSTss解析実行スクリプト - ターミナル出力をファイルに保存

if [ $# -lt 1 ]; then
    echo "Usage: ./run_analysis.sh input_name [output_dir]"
    echo "Example: ./run_analysis.sh 1elem_f"
    echo "         ./run_analysis.sh 1elem_f output/my_logs"
    exit 1
fi

INPUT_NAME=$1
OUTPUT_DIR=${2:-"output/analysis_logs"}
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_FILE="${OUTPUT_DIR}/${INPUT_NAME}_${TIMESTAMP}.txt"

# 出力ディレクトリ作成
mkdir -p "$OUTPUT_DIR"

echo "=== Running PLSTss Analysis ==="
echo "Input: $INPUT_NAME"
echo "Output will be saved to: $OUTPUT_FILE"
echo "========================================"

# teeコマンドを使用してコンソールとファイルに同時出力
echo "$INPUT_NAME" | ~/bin/plstss2393_nagasaku 2>&1 | tee "$OUTPUT_FILE"

echo ""
echo "========================================"
echo "Analysis complete. Output saved to: $OUTPUT_FILE"
echo ""
echo "Generated result files:"
for file in STS DIS NOR ENE TMP RES; do
    if [ -f "output/${file}_${INPUT_NAME}.txt" ] || [ -f "output/${file}_${INPUT_NAME}.cml" ]; then
        echo "  ✓ output/${file}_${INPUT_NAME}.*"
    fi
done