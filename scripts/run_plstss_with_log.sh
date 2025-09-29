#!/bin/bash
# PLSTss実行スクリプト - 出力をターミナルとログファイルに同時保存

if [ $# -lt 1 ]; then
    echo "Usage: ./run_plstss_with_log.sh input_name"
    echo "Example: ./run_plstss_with_log.sh 1elem_f"
    exit 1
fi

INPUT_NAME=$1
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="output/LOG_${INPUT_NAME}_${TIMESTAMP}.txt"

echo "Running PLSTss analysis for: $INPUT_NAME"
echo "Log will be saved to: $LOG_FILE"
echo "========================================="

# teeコマンドでターミナルとファイルに同時出力
echo "$INPUT_NAME" | ~/bin/plstss2393_nagasaku 2>&1 | tee "$LOG_FILE"

echo ""
echo "========================================="
echo "Analysis complete!"
echo "Log saved to: $LOG_FILE"