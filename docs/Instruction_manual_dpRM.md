# Drucker–Prager 戻り写像コード修正 指示書（Claude Code 用）

**対象ファイル**: `stress_dp_rm.f`（Drucker–Prager：戻り写像）  
**前提**: 弾性パートの `stress.f` は正しい（`vkp`=体積弾性率K, `vmu`=せん断弾性率μ の定義と利用は正）。  
**目的**: 体積成分と偏差成分の厳密な分離、純静水圧トライアル時の特異点回避、塑性仕事の式の是正、および一貫接線の安全化を行い、数値安定性と物理整合性を高める。

---

## 0. 作業準備
- [ ] 作業ブランチを作成: `git checkout -b fix/dp-return-mapping`
- [ ] バックアップ: `cp stress_dp_rm.f stress_dp_rm.before_fix.f`
- [ ] ファイル内の主要シンボルの実体を確認：`vkp, vmu, DELTA(3x3), str(3x3), plstrg(3x3), etabar_dp, eta_dp, xi_dp, yld, hrdtmp, dhdtmp, deltag, stno, oun(3x3), ctens(3,3,3,3), gg, Dg, A, theta, thetab`  
  （命名が多少異なる場合は**意味**で対応するものに読み替えること）

---

## 1. 体積・偏差の“定義”を明確化（必須）
**修正方針**: 「弾性トライアル量＝総ひずみ − 過去塑性ひずみ」で定義し、体積（I）と偏差（dev）を分離する。

### 1.1 平均（体積）ひずみ
以下の2量を導入・利用する：
```fortran
emean  = (str(1,1)+str(2,2)+str(3,3)) / 3.d0
epmean = (plstrg(1,1)+plstrg(2,2)+plstrg(3,3)) / 3.d0
```

### 1.2 試行“体積弾性ひずみ” `etrs`
> **置換**: 旧 `etrs = tr(str)` → 新 `etrs = tr(str - plstrg)`
```fortran
etrs = (str(1,1)+str(2,2)+str(3,3)) - (plstrg(1,1)+plstrg(2,2)+plstrg(3,3))
```

### 1.3 試行“偏差応力” `stry`
> **意図**: `s_tr = 2μ * ( (ε - ε^p_n)_dev )`。塑性ひずみの**体積成分**は dev から除外。
```fortran
! 偏差部分: X_dev = X - mean(X)*DELTA
stry = 2.d0*vmu * ( (str - emean*DELTA) - (plstrg - epmean*DELTA) )
```

---

## 2. 偏差ノルムと単位テンソル `oun`
```fortran
stno = 0.d0
Do j=1,3; Do i=1,3
  stno = stno + stry(i,j)*stry(i,j)
End Do; End Do
stno = dsqrt(stno)

tiny = 1.0d-14
If (stno .gt. tiny) Then
  oun(:,:) = stry(:,:) / stno
Else
  oun(:,:) = 0.d0  ! 純静水圧近傍：偏差方向は未定義→0 を採用
End If
```

---

## 3. 一致条件（Newton 更新）分岐の是正
**狙い**: 純静水圧（`stno≈0`）では偏差寄与を消し、体積＋硬化のみで Δγ を解く。

### 3.1 既存コードの `gg` / `Dg` を置換
- **偏差あり (`stno>tiny`)**:
```fortran
! gg = 0 を Newton で解いて Δγ を更新（Δγ ≥ 0 を強制）
! hrdtmp = H(α_n), dhdtmp = H'(α_n) の意味

gg = dsqrt(1.d0/2.d0)*stno - vmu*deltag  &
   + eta_dp*(vkp*etrs - vkp*etabar_dp*deltag)  &
   - xi_dp*(yld + hrdtmp)

Dg = -vmu - eta_dp*vkp*etabar_dp - xi_dp*xi_dp*dhdtmp
```
- **純静水圧近傍 (`stno≤tiny`)**:
```fortran
! 偏差寄与（-vmu*deltag）は入れない

gg =        eta_dp*(vkp*etrs - vkp*etabar_dp*deltag)  &
   - xi_dp*(yld + hrdtmp)

Dg = -      eta_dp*vkp*etabar_dp - xi_dp*xi_dp*dhdtmp
```
- **Newton 更新の安全化**:
```fortran
! 例: 単調減少ステップ、Δγの非負制約
step = gg / Dg
if (Dg .eq. 0.d0) step = 0.d0
new_dg = max(0.d0, deltag - step)
deltag = new_dg
```

> **注**: 既存の `A, theta, thetab` を用いた記法を残す場合も、上記の**分岐思想**（`stno` に応じて偏差寄与を抑制）を必ず反映させる。

---

## 4. 応力更新式の整合化
**目的**: 弾性と完全整合する形で球・偏差を更新。
```fortran
sig(:,:) =  stry(:,:) - dsqrt(2.d0)*vmu*deltag*oun(:,:)  &
          + (vkp*etrs - vkp*etabar_dp*deltag)*DELTA(:,:)
```
- `trace(sig) = 3 * (vkp*etrs - vkp*etabar_dp*deltag)` となること（偏差項はトレース0）。

---

## 5. 塑性ひずみ更新（流れ則）
**目的**: 体積・偏差の寄与係数を厳密化（係数抜けを解消）。
```fortran
plstrg(:,:) = plstrg(:,:)  &
  + deltag * ( dsqrt(1.d0/2.d0)*oun(:,:) + (etabar_dp/3.d0)*DELTA(:,:) )
```

---

## 6. 塑性仕事密度 `p_dns` の式（係数是正）
**旧**: `p_dns += deltag * oun(i,j)*sig(i,j)` （→係数不足＆体積寄与欠落）

**新**:
```fortran
tmp = 0.d0
Do j=1,3; Do i=1,3
  tmp = tmp + sig(i,j) * ( dsqrt(1.d0/2.d0)*oun(i,j) + (etabar_dp/3.d0)*DELTA(i,j) )
End Do; End Do
p_dns = p_dns + deltag * tmp
```

---

## 7. 一貫接線 `ctens` の安全化
- **弾性側**（降伏外）: 既存の `K I⊗I + 2μ (I⁴ − 1/3 I⊗I)` を維持。
- **塑性側**: 既存式を尊重しつつ、**`stno≤tiny` 分岐**で偏差寄与の塑性項を**ゼロ化**（体積側のみ）。
  - 例: `oun` を含む項は `stno≤tiny` では寄与させない。
  - `theta, thetab` を利用する計算系なら、`stno≤tiny` のとき `theta=1.d0`、`thetab=-vmu*A` 等の安全値にフォールバック。
- **NaN/Inf フェイルセーフ**：更新直後に `isnan`/`isinf` 相当のチェックを入れ、異常時は弾性接線へフォールバック（ログを出す）。

---

## 8. 命名・可読性（任意だが推奨）
- `etrs` → `evol_tr`（trial volumetric strain）
- `emean, epmean` を導入して dev 分離を常時**見える化**。
- `plstrg_dev = plstrg - epmean*DELTA` を明示変数に（誤用防止）。

---

## 9. 受け入れ基準（Acceptance Criteria）
1) **静水圧試験**（`ε = κ`倍 I のみ）
- [ ] `stno→0` で特異が発生しない（`oun=0`）。
- [ ] Δγ は体積寄与のみで決定される。
- [ ] `sig = K*(etrs - etabar*Δγ) I` を満たす。

2) **純せん断**（体積0、偏差のみ）
- [ ] 体積塑性は発生しない（`trace(plstrg)` の増分≈0）。
- [ ] 降伏後、偏差応力の大きさは戻り写像で一貫（`f≈0`）。

3) **単軸圧縮/引張**
- [ ] 4ステップ目の応力状態が降伏面上に収束（|f| < 許容誤差）。
- [ ] 弾性ステップでは旧版と一致（回帰）。

4) **数値安定性**
- [ ] Δγ≥0 が常に満たされる（Newton 更新で強制）。
- [ ] `|oun|=1`（`stno>tiny`）が維持される。
- [ ] NaN/Inf 検出時は弾性接線へフェイルセーフ。

---

## 10. テスト実行（例）
- 既存の単体ドライバ／要素テストがあればそれを流用。なければ最小限のドライバを用意：
  - **Case A**: 静水圧圧縮（`str = diag(κ,κ,κ)` を時間発展）
  - **Case B**: 純せん断（`str = [[0,γ/2,0],[γ/2,0,0],[0,0,0]]`）
  - **Case C**: 単軸（`str = diag(ε,0,0)`）
- 収束判定: `|gg| < 1e-10`、`|Δγ_k - Δγ_{k-1}| < 1e-12` など。
- ログ出力: 主要量（`stno, etrs, deltag, trace(sig), ||s||, f`）を各ステップで出力。

---

## 11. コンパイル/警告レベル
- Fortran: 可能なら `implicit none` を有効化。
- 例: `ifx/ifort` なら `-warn all -stand f18 -fpe0`（環境に合わせて）。

---

## 12. 変更概要（コミットメッセージ雛形）
```
DP: fix return-mapping (dev/vol split, hydrostatic branch, plastic work)
- etrs := tr(eps - eps_p)  (trial volumetric strain)
- stry := 2μ*((eps-eps_p)_dev)  (dev-only)
- hydrostatic branch: stno→0 ⇒ oun=0, gg/Dg without dev terms
- stress update: s_dev & p parts consistent with elastic law
- plastic strain/work: correct coefficients (√1/2, etabar/3)
- tangent: safe fallback & branch-specific updates
- add guards: Δγ≥0, |oun|=1, NaN/Inf fallback
```

---

## 13. 注意点
- 既存の `A, theta, thetab` 記法を残す場合は、**本指示の分岐ロジック**（偏差寄与の on/off）を忠実に反映させる。
- 体積塑性の取り扱いはこの指示に統一（`etrs = tr(ε-ε^p)`、`Δε^p_vol = (etabar/3)Δγ I`）。
- `stress.f` と定数・記号の整合（`vkp, vmu, DELTA`）を常に確認。

---

## 14. 完了チェック
- [ ] すべての置換/追加コードにコメントを付与（`! DP FIX:` プレフィクス）
- [ ] 3種類のテストケースで合格（§9）
- [ ] 主要関数に回帰テストを追加（弾性域完全一致）
- [ ] ブランチにコミット＆PR 作成

---

### 付録：参考実装スニペット（貼り付け用）
> **目的**: 主要箇所の置換テンプレート（行番号依存にしない）。
```fortran
! --- DP FIX: volumetric/deviatoric split ----------------------------------
emean  = (str(1,1)+str(2,2)+str(3,3)) / 3.d0
epmean = (plstrg(1,1)+plstrg(2,2)+plstrg(3,3)) / 3.d0

etrs = (str(1,1)+str(2,2)+str(3,3)) - (plstrg(1,1)+plstrg(2,2)+plstrg(3,3))

stry(:,:) = 2.d0*vmu * ( (str(:,:) - emean*DELTA(:,:)) - (plstrg(:,:) - epmean*DELTA(:,:)) )

stno = 0.d0
Do j=1,3; Do i=1,3
  stno = stno + stry(i,j)*stry(i,j)
End Do; End Do
stno = dsqrt(stno)

tiny = 1.0d-14
If (stno .gt. tiny) Then
  oun(:,:) = stry(:,:) / stno
Else
  oun(:,:) = 0.d0
End If

! --- DP FIX: consistency (Newton) skeleton --------------------------------
If (stno .gt. tiny) Then
  gg = dsqrt(1.d0/2.d0)*stno - vmu*deltag  &
     + eta_dp*(vkp*etrs - vkp*etabar_dp*deltag)  &
     - xi_dp*(yld + hrdtmp)
  Dg = -vmu - eta_dp*vkp*etabar_dp - xi_dp*xi_dp*dhdtmp
Else
  gg =        eta_dp*(vkp*etrs - vkp*etabar_dp*deltag)  &
     - xi_dp*(yld + hrdtmp)
  Dg = -      eta_dp*vkp*etabar_dp - xi_dp*xi_dp*dhdtmp
End If

! Newton update example
step   = gg / Dg
If (Dg .eq. 0.d0) step = 0.d0
new_dg = max(0.d0, deltag - step)
deltag = new_dg

! --- DP FIX: stress & plastic updates -------------------------------------
sig(:,:) =  stry(:,:) - dsqrt(2.d0)*vmu*deltag*oun(:,:)  &
          + (vkp*etrs - vkp*etabar_dp*deltag)*DELTA(:,:)

plstrg(:,:) = plstrg(:,:) + deltag * ( dsqrt(1.d0/2.d0)*oun(:,:) + (etabar_dp/3.d0)*DELTA(:,:) )

tmp = 0.d0
Do j=1,3; Do i=1,3
  tmp = tmp + sig(i,j) * ( dsqrt(1.d0/2.d0)*oun(i,j) + (etabar_dp/3.d0)*DELTA(i,j) )
End Do; End Do
p_dns = p_dns + deltag * tmp
```

> **以上**。上記に従って `stress_dp_rm.f` を修正し、§9 の受け入れ基準を満たすこと。

