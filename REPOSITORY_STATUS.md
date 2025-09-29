# Repository Status Summary

## 現在の状態 (2025-09-30)

### Gitリポジトリ
- **現在のブランチ**: main
- **最新コミット**: 5a6639c "Fix Drucker-Prager consistent tangent and file handling issues"
- **2つ前のコミット**: 906e459 "LOVE CLAUDE CODE"

### stress_dp_rm.f の状態

**現在のmainブランチ (5a6639c) - 理論的に正しい定式化:**
```fortran
c Line 163: 正しいハードニング項
sqrt3 = dsqrt(3.d0)
Dg = -vmu -eta_dp*vkp*etabar_dp -xi_dp*dhdtmp/sqrt3

c Line 171: 正しいα更新
alptmp = alpeg + deltag/sqrt3

c Line 262: 正しいα最終更新
alpeg = alpeg + deltag/sqrt3

c Line 308: 正しい一貫接線係数
A = 1.d0/(vmu + vkp*etabar_dp*eta_dp + xi_dp*dhard/sqrt3)
```

**906e459 (LOVE CLAUDE CODE) - 理論的に誤りだが動作する:**
```fortran
c 誤ったハードニング項
Dg = -vmu -eta_dp*vkp*etabar_dp -xi_dp*xi_dp*dhdtmp

c 誤ったα更新
alptmp = alpeg + xi_dp*deltag

c 誤った一貫接線係数
A = 1.d0/(vmu + vkp*etabar_dp*eta_dp + xi_dp*xi_dp*dhard)
```

### 違いと選択

| 項目 | 906e459 (LOVE CLAUDE CODE) | 5a6639c (現在のmain) |
|------|---------------------------|---------------------|
| 理論的正しさ | ❌ 誤り (H=ξ²×∂K/∂α) | ✅ 正しい (H=ξ×∂K/∂α/√3) |
| 非関連流れ則の収束 | ✅ 収束する (偶然) | ❌ 収束しない |
| 収束の種類 | 線形収束 (16-17反復) | - |
| 根本原因 | 誤った定式化が偶然に非対称性を減少 | 正しい非対称接線が対称格納と不整合 |

### 推奨事項

1. **理論的正確性を重視する場合**:
   - 現在のmainブランチ (5a6639c) を維持
   - ただし非関連流れ則 (φ≠ψ) は使用不可
   - 関連流れ則 (φ=ψ) のみ使用

2. **非関連流れ則を使いたい場合**:
   - 選択肢A: マトリックス格納を非対称に拡張（大規模改修必要）
   - 選択肢B: 906e459の「誤った」定式化を使用（理論的には誤りだが動作）

3. **現在の選択**:
   - mainブランチは理論的に正しい実装を保持
   - リモートサーバーも同期済み
   - 非関連流れ則が必要な場合は別途検討が必要

### ログファイル保存済み
- `output/logs/LOVE_CLAUDE_phi10psi5.log` - 元のコード、収束
- `output/logs/LOVE_CLAUDE_phi20psi15.log` - 元のコード、収束
- `output/logs/corrected_phi10psi5.log` - 修正コード、収束失敗
- `output/logs/corrected_phi20psi15.log` - 修正コード、収束失敗
- `output/COMPARISON_ORIGINAL_VS_CORRECTED.md` - 詳細な比較