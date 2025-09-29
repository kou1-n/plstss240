# `stress_dp_rm` 修正設計ノート（Drucker–Pragerの体積塑性を整合化）

**目的**  
既存 `stress_dp_rm.f` は J2 前提の名残として、降伏判定・残差・応力更新の体積項に **全体積ひずみ** \( \varepsilon_v = \mathrm{tr}\,\boldsymbol\varepsilon \) をそのまま用いています。  
Drucker–Prager（DP）では **塑性体積ひずみ** \( \varepsilon_v^p = \mathrm{tr}\,\boldsymbol\varepsilon^p \) が非ゼロとなるため、**平均応力は弾性体積ひずみ**  
\[
p = K\,\varepsilon_v^{e} = K\bigl(\varepsilon_v - \varepsilon_v^{p}\bigr)
\]
から求める必要があります。

以下では、**既存コードを最大限活かしつつ**、`p_tr = K(\varepsilon_v - \varepsilon_v^{p(n)})` を導入して、試行降伏関数／Newton 残差／応力更新に一貫して反映させる最小修正を示します。

---

## 変更点サマリ（箇条書き）

1. **試行平均応力**の導入  
   ```fortran
   epv0 = plstrg(1,1) + plstrg(2,2) + plstrg(3,3)
   p_tr = vkp*(etrs - epv0)     ! K * (ε_v - ε_v^p)
   ```
   - `etrs` は既存コード（全体積ひずみ）  
   - `epv0 = tr(ε^p_n)`（履歴から）

2. **試行降伏関数 `ftreg`** の体積項置換  
   **旧**：`eta_dp*vkp*etrs` → **新**：`eta_dp*p_tr`

3. **Newton 残差 `gg`** の体積項置換  
   **旧**：`eta_dp*(vkp*etrs - vkp*etabar_dp*deltag)`  
   **新**：`eta_dp*(p_tr - vkp*etabar_dp*deltag)`

4. **応力更新（塑性時）** の体積項置換  
   **旧**：`+ (vkp*etrs - vkp*etabar_dp*deltag) * DELTA`  
   **新**：`+ (p_tr - vkp*etabar_dp*deltag) * DELTA`

5. **応力更新（弾性時）** の体積項置換  
   **旧**：`sig = stry + vkp*etrs*DELTA`  
   **新**：`sig = stry + p_tr*DELTA`

6. **塑性仕事密度 `p_dns`**（任意だが推奨）  
   **旧**：`p_dns += deltag * (oun : sig)`  
   **新**：\[ \sigma:\mathrm{d}\varepsilon^p = \Delta\lambda(\sqrt{1/2}\,\|s\| + \bar\eta\,p) \] に合わせ、
   ```fortran
   smean = (sig(1,1)+sig(2,2)+sig(3,3))/3.d0
   stno_u = dsqrt( s(1,1)**2 + ... + s(3,3)**2 )      ! 更新後 ||s||
   p_dns  = deltag*( dsqrt(0.5d0)*stno_u + etabar_dp*smean )
   ```

7. **一貫接線 `ctens`** は現行の式のままで整合的  
   `Dg = -μ - ηK\barη - ξ^2 H'` を使っており、今回の `p_tr` 置換と整合。直接 `p_tr` を参照していないため**変更不要**です。

---

## 具体的な変更（検索アンカーごとのパッチ）

### 1) 体積項の前処理（`etrs` 直後）

- **検索アンカー**（すでに存在）  
  ```fortran
  etrs = str(1,1) + str(2,2) + str(3,3)
  ```

- **すぐ下に追記**  
  ```fortran
  epv0 = plstrg(1,1) + plstrg(2,2) + plstrg(3,3)    ! tr(ε^p_n)
  p_tr = vkp*(etrs - epv0)                          ! K(ε_v - ε_v^p)
  ```

> 参考：既に偏差応力の試行値は
> \(
\boldsymbol{s}^{tr} = 2\mu \bigl[(\boldsymbol\varepsilon - \mathrm{mean}\,\boldsymbol\varepsilon\mathbf{I})
 - (\boldsymbol\varepsilon^p - \mathrm{mean}\,\boldsymbol\varepsilon^p\mathbf{I})\bigr]
\)
として**塑性体積成分を除外**しているので、その設計方針に**平均応力側も揃える**のが本修正の意図です。

---

### 2) 試行降伏関数 `ftreg`

- **旧**（例）
  ```fortran
  ftreg = dsqrt(1.d0/2.d0)*stno + eta_dp*vkp*etrs 
  &      - xi_dp*(yld + hard)
  ```

- **新**
  ```fortran
  ftreg = dsqrt(1.d0/2.d0)*stno + eta_dp*p_tr
  &      - xi_dp*(yld + hard)
  ```

**理由**：  
\(
f = \sqrt{J_2} + \eta\,p - \xi\bigl(\sigma_y + H(\alpha)\bigr),\;
p = K(\varepsilon_v - \varepsilon_v^p)
\)  
DP では \( \varepsilon_v^p\neq 0 \)。J2 のように \(p=K\varepsilon_v\) とできない。

---

### 3) Newton 残差 `gg`

- **旧**
  ```fortran
  gg = dsqrt(1.d0/2.d0)*stno 
  &  - vmu*deltag
  &  + eta_dp*(vkp*etrs - vkp*etabar_dp*deltag)
  &  - xi_dp*(yld + hrdtmp)
  ```

- **新**
  ```fortran
  gg = dsqrt(1.d0/2.d0)*stno 
  &  - vmu*deltag
  &  + eta_dp*(p_tr - vkp*etabar_dp*deltag)
  &  - xi_dp*(yld + hrdtmp)
  ```

**理由（式）**：  
\(
p_{n+1} = K\bigl(\varepsilon_v - \varepsilon_v^{p(n)} - \bar\eta\,\Delta\lambda\bigr)
       = p_{tr} - K\bar\eta\,\Delta\lambda
\)  
\(
g(\lambda)=\sqrt{\tfrac12}\|s^{tr}\| - \mu\lambda + \eta\,(p_{tr}-K\bar\eta\,\lambda)
 - \xi\bigl(\sigma_y+H(\alpha_n+\xi\lambda)\bigr)
\)

---

### 4) 応力更新（塑性時）

- **旧**
  ```fortran
  sig(:,:) = stry(:,:) - dsqrt(2.d0)*vmu*deltag*oun(:,:)
  &        + (vkp*etrs - vkp*etabar_dp*deltag)*DELTA(:,:)
  ```

- **新**
  ```fortran
  sig(:,:) = stry(:,:) - dsqrt(2.d0)*vmu*deltag*oun(:,:)
  &        + (p_tr - vkp*etabar_dp*deltag)*DELTA(:,:)
  ```

**理由（式）**：  
\(
\boldsymbol{s}_{n+1}=\boldsymbol{s}^{tr}-\sqrt{2}\mu\,\Delta\lambda\,\mathbf{n},\quad
p_{n+1}=p_{tr}-K\bar\eta\,\Delta\lambda,\quad
\boldsymbol{\sigma}_{n+1}=\boldsymbol{s}_{n+1}+p_{n+1}\mathbf{I}
\)

---

### 5) 応力更新（弾性時）

- **旧**
  ```fortran
  sig(:,:) = stry(:,:) + vkp*etrs*DELTA(:,:)
  ```

- **新**
  ```fortran
  sig(:,:) = stry(:,:) + p_tr*DELTA(:,:)
  ```

**理由**：弾性でも \(p=K(\varepsilon_v-\varepsilon_v^{p(n)})\)。\( \varepsilon_v^p \) を無視しない。

---

### 6) 塑性仕事密度 `p_dns`（推奨置換）

- **旧**
  ```fortran
  p_dns = 0.d0
  do jj=1,3; do ii=1,3
    p_dns = p_dns + deltag*oun(ii,jj)*sig(ii,jj)
  enddo; enddo
  ```

- **新**
  ```fortran
  smean  = (sig(1,1)+sig(2,2)+sig(3,3))/3.d0      ! = p_{n+1}
  stno_u = 0.d0
  do jj=1,3; do ii=1,3
    stno_u = stno_u + s(ii,jj)*s(ii,jj)           ! s は更新後の偏差応力
  enddo; enddo
  stno_u = dsqrt(stno_u)
  p_dns  = deltag*( dsqrt(0.5d0)*stno_u + etabar_dp*smean )
  ```

**理由（式）**：  
非関連ポテンシャル \( g=\sqrt{J_2}+\bar\eta\,p \) より  
\(
\mathrm{d}\boldsymbol\varepsilon^p=\Delta\lambda\left(\sqrt{\tfrac12}\mathbf{n}
+\tfrac{\bar\eta}{3}\mathbf{I}\right)\Rightarrow
\sigma:\mathrm{d}\varepsilon^p=\Delta\lambda\bigl(\sqrt{\tfrac12}\|s\|+\bar\eta\,p\bigr)
\)

---

## 数式による一貫性の確認

- **降伏関数**（実装形）  
  \[
  f = \sqrt{J_2} + \eta\,p - \xi\bigl(\sigma_y + H(\alpha)\bigr),\quad
  p = K(\varepsilon_v-\varepsilon_v^p)
  \]
  ここで \( \varepsilon_v^p\neq 0 \) なのが J2 との本質的相違。

- **流動則**（非関連）  
  \[
  g = \sqrt{J_2}+\bar\eta\,p,\quad
  \mathrm{d}\varepsilon^p=\Delta\lambda\left(\sqrt{\tfrac12}\mathbf{n}
  +\tfrac{\bar\eta}{3}\mathbf{I}\right)
  \]

- **体積応力更新**  
  \[
  p_{n+1} = K\bigl(\varepsilon_v - \varepsilon_v^{p(n)} - \bar\eta\,\Delta\lambda\bigr)
           = p_{tr} - K\bar\eta\,\Delta\lambda
  \]

- **Newton 残差**  
  \[
  g(\lambda)=\sqrt{\tfrac12}\|s^{tr}\| - \mu\lambda
             + \eta\,(p_{tr}-K\bar\eta\,\lambda)
             - \xi\bigl(\sigma_y+H(\alpha_n+\xi\lambda)\bigr)
  \]
  \[
  g'(\lambda)= -\mu - \eta K\bar\eta - \xi^2 H'(\alpha_n+\xi\lambda)
  \]
  ※ 既存実装の `Dg` と一致（**接線は変更不要**）。

---

## 追加のデバッグ出力（任意）

```fortran
write(*,'(A,1PE12.5)') ' tr(plstrg) (epv0) = ', epv0
write(*,'(A,1PE12.5)') ' p_tr = K(εv - εvp) = ', p_tr
```

- 弾性ステップ：`sig = stry + p_tr*DELTA` で `tr(plstrg)` の影響が入ることを確認。
- 塑性ステップ：`p_now = p_tr - K*etabar_dp*deltag` が負圧縮側で増減する様子を確認。

---

## 簡易検証シナリオ

1. **J2 退化**：`phi=psi=0` → `eta=etabar=0` で元の J2 に一致。  
2. **単軸圧縮（φ,ψ>0）**：`tr(plstrg)` が正に蓄積し、`p_tr` がそれを差し引くため、
   弾性時と比べて平均応力が小さく更新されることを確認。  
3. **apex 近傍**：`stno→0` での収束を確認（既存の保護分岐があればそのまま利用）。

---

## 参考（今回**変更しない**箇所）

- `Dg`、`ctens`（一貫接線）ブロックは、今回の \(p_{tr}\) 置換で**式構造が不変**のため変更不要。  
- 偏差方向ベクトル `oun` の更新ロジック（試行偏差応力ベース）は現状のままで問題なし。

---

## 付録：挿入・置換スニペットまとめ

```fortran
! --- after computing etrs ---
epv0 = plstrg(1,1) + plstrg(2,2) + plstrg(3,3)
p_tr = vkp*(etrs - epv0)

! --- ftreg ---
ftreg = dsqrt(1.d0/2.d0)*stno + eta_dp*p_tr
&      - xi_dp*(yld + hard)

! --- residual gg ---
gg = dsqrt(1.d0/2.d0)*stno - vmu*deltag
&  + eta_dp*(p_tr - vkp*etabar_dp*deltag)
&  - xi_dp*(yld + hrdtmp)

! --- plastic stress update ---
sig(:,:) = stry(:,:) - dsqrt(2.d0)*vmu*deltag*oun(:,:)
&        + (p_tr - vkp*etabar_dp*deltag)*DELTA(:,:)

! --- elastic stress update ---
sig(:,:) = stry(:,:) + p_tr*DELTA(:,:)

! --- plastic work density (optional but recommended) ---
smean  = (sig(1,1)+sig(2,2)+sig(3,3))/3.d0
stno_u = 0.d0
do j=1,3; do i=1,3
  stno_u = stno_u + s(i,j)*s(i,j)
enddo; enddo
stno_u = dsqrt(stno_u)
p_dns  = deltag*( dsqrt(0.5d0)*stno_u + etabar_dp*smean )
```
