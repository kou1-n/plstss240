#!/bin/bash
# PLSTss一括解析実行スクリプト
# 使用方法: ./scripts/batch_analysis.sh [input_files_pattern]
# 例: ./scripts/batch_analysis.sh "1elem_*"

echo "=========================================="
echo "PLSTss Batch Analysis Tool"
echo "=========================================="

# デフォルトパターン（引数なしの場合、すべての.cmlファイル）
PATTERN=${1:-"*"}

# 実行可能ファイルのパス
PLSTSS_EXE="$HOME/bin/plstss2393_nagasaku"

# 解析するファイルのリストを取得
FILES=$(ls input_files/${PATTERN}.cml 2>/dev/null)

if [ -z "$FILES" ]; then
    echo "Error: No input files found matching pattern: ${PATTERN}.cml"
    echo "Available input files:"
    ls input_files/*.cml 2>/dev/null | sed 's|input_files/||g' | sed 's|.cml||g'
    exit 1
fi

# ファイル数をカウント
FILE_COUNT=$(echo "$FILES" | wc -l)
echo "Found $FILE_COUNT input file(s) to analyze:"
echo "$FILES" | sed 's|input_files/||g' | sed 's|.cml||g'
echo ""

# 解析開始時刻
START_TIME=$(date +%s)
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# ログディレクトリ作成
LOG_DIR="output/logs/batch_${TIMESTAMP}"
mkdir -p "$LOG_DIR"

echo "Log directory: $LOG_DIR"
echo "=========================================="

# 各ファイルに対して解析を実行
SUCCESS_COUNT=0
FAIL_COUNT=0

for FILE in $FILES; do
    # ファイル名を抽出（パスと拡張子を除く）
    BASENAME=$(basename "$FILE" .cml)

    echo ""
    echo "[$((SUCCESS_COUNT + FAIL_COUNT + 1))/$FILE_COUNT] Analyzing: $BASENAME"
    echo "----------------------------------------"

    # 入力ファイルを作業ディレクトリにコピー
    cp "$FILE" "${BASENAME}.cml"

    # 解析実行（ログファイルに保存）
    LOG_FILE="$LOG_DIR/${BASENAME}.log"
    echo "$BASENAME" | $PLSTSS_EXE > "$LOG_FILE" 2>&1

    # 実行結果をチェック
    if [ $? -eq 0 ] && grep -q "Program \"PLSTss\" has done" "$LOG_FILE" 2>/dev/null; then
        echo "✓ Success: Analysis completed"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))

        # 結果ファイルを移動
        if [ -d "output/results/latest" ]; then
            mv output/STS_${BASENAME}.txt output/results/latest/ 2>/dev/null
            mv output/DIS_${BASENAME}.txt output/results/latest/ 2>/dev/null
            mv output/NOR_${BASENAME}.txt output/results/latest/ 2>/dev/null
            mv output/ENE_${BASENAME}.txt output/results/latest/ 2>/dev/null
            mv output/TMP_${BASENAME}.txt output/results/latest/ 2>/dev/null
            mv output/RES_${BASENAME}.cml output/results/latest/ 2>/dev/null
        fi
    else
        echo "✗ Failed: Check log file for errors"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi

    # 作業ディレクトリの入力ファイルを削除
    rm -f "${BASENAME}.cml"
done

# 解析完了
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "=========================================="
echo "Batch Analysis Complete!"
echo "=========================================="
echo "Total files: $FILE_COUNT"
echo "Successful: $SUCCESS_COUNT"
echo "Failed: $FAIL_COUNT"
echo "Time elapsed: $ELAPSED seconds"
echo ""
echo "Logs saved to: $LOG_DIR"
echo "Results saved to: output/results/latest/"
echo "=========================================="