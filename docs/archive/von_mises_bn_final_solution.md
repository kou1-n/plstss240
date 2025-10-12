# Von Mises Block Newton 最終解決策

## 実施した修正

### 1. ダンピング機構（phexa8.f）
- 初回塑性ステップで1%に減衰
- 絶対値制限を1.d-4に設定

### 2. N_scalar相対閾値（stress_dp_bn.f）
- 弾性値の0.1%を閾値に変更
- N_elastic = -2μを基準

### 3. 初回反復のg_val減衰（stress_dp_bn.f）
- 初回塑性でg_valを1%に減衰

## 結果
- **改善**: ||Rg||爆発が4317→16094と増大（悪化）
- **結論**: 単純なダンピングでは不十分

## 根本的問題

Block Newton法の実装が**分離型**であること：
1. stress_dp_bn.fがg_val、N_scalarを計算
2. phexa8.fが独立にδγを計算
3. 相互作用が**逐次的**で同時収束しない

## 真の解決策

### オプション1: Line Search実装
```fortran
! phexa8.fまたはstress_dp_bn.fに追加
do i_line = 1, 5
  sig_trial = sig_old - alpha*delta_gamma_inc*n
  g_trial = compute_yield(sig_trial)
  if(abs(g_trial) < abs(g_old)) exit
  alpha = 0.5*alpha
enddo
delta_gamma_inc = alpha*delta_gamma_inc
```

### オプション2: Substepping
弾塑性遷移時に荷重増分を細分化：
```fortran
if(first_plastic) then
  df_local = df * 0.1  ! 10分割
endif
```

### オプション3: Return Mappingへの切り替え
stress_dp_rm.fは既に動作確認済み。Block Newtonの修正は困難。

## 推奨

**Von Misesでさえ収束しない**ことから、現在のBlock Newton実装は根本的に欠陥がある。以下を推奨：

1. **短期的**: Return Mapping（stress_dp_rm.f）を使用
2. **長期的**: 真のBlock Newton（Schur補元による同時解法）を実装

## 技術的診断

現在の実装の致命的欠陥：
- δγ = -(M:ε + g)/N が**unbounded**
- g_valとδγの**フィードバックループ**が不安定
- **同時収束**ではなく**交互更新**

Von Misesという最も単純なケースでも失敗することは、アルゴリズムの根本的再設計が必要であることを示している。