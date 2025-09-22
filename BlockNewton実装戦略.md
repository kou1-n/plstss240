# Block Newton法による弾塑性解析の実装戦略
## 山本先生の論文（2021）に基づくstress_dp_bn.fコーディング方針

---

## 1. 論文の核心概念

### 1.1 Block Newton法の特徴
山本ら(2021)の論文「Simultaneously iterative procedure based on block Newton method for elastoplastic problems」は、弾塑性問題を**平衡方程式と降伏条件の連成問題**として定式化し、Block Newton法により同時に解く革新的手法である。

#### 従来法（Return Mapping）との主要な違い：
- **従来法**: 各積分点で局所反復計算により内部変数（consistency parameter Δγ）を決定
- **Block Newton法**: 内部変数を代数的に消去し、大域的Newton反復のみで解を得る
- **利点**: 平衡方程式と降伏条件の誤差が同時に減少

### 1.2 数学的定式化
連成問題の基本構造（論文式45-46）：
```
∇·{C : ε(δu) + Lδγ} = f_{n+1} - ∇·σ     (平衡方程式)
M : ε(δu) + Nδγ = -g                      (降伏条件)
```
ここで：
- **C**: 接線係数テンソル（幾何学的軟化含む）
- **L, M**: 流動則に関する係数
- **N**: 硬化則に関するスカラー係数

---

## 2. BOX 1アルゴリズム（応力更新手順）

### 2.1 手順概要
BOX 1は各材料点における応力更新の詳細手順を示す。

```fortran
! STEP 1: Trial elastic stress計算
s^tr = 2μ{dev[ε] - ε^p_n}
ξ^tr = s^tr - β_n
n = ξ^tr / |ξ^tr|

! STEP 2: Consistency parameter Δγのチェック
IF (Δγ > 0) THEN
  ! 塑性状態
  ! 2.1 状態変数更新
  ε^p = ε^p_n + Δγ·n
  α = α_n + √(2/3)·Δγ
  β = β_n + (2/3)H'Δγ·n

  ! 2.2 応力更新
  σ = κtr[ε]I + s^tr - 2μΔγ·n

  ! 2.3 Stress corrector計算
  g = |ξ| - √(2/3)[σ_Y + K(α)]
  σ_g = -(N^{-1}·g)·L
ELSE
  ! 弾性状態
  σ = C^e : ε
  g ≡ 0  (論文式84)
  σ_g ≡ 0
ENDIF

! STEP 3: Tangent moduli計算
C^ep = C - N^{-1}L⊗M
```

### 2.2 実装上の重要点
1. **初回反復時（itr=1）**: BOX 2の初期化手順に従う
2. **幾何学的軟化**: θ = 1 - (2μΔγ)/|ξ^tr|を考慮
3. **特異性回避**: N_scalarが小さすぎる場合の処理

---

## 3. BOX 2アルゴリズム（有限要素解析の変数更新）

### 3.1 手順概要
BOX 2は有限要素解析全体の反復手順を示す。

```fortran
! STEP 1: 初期化
Δγ^(0) = 0
g^(0) = 0
ε^p = ε^p_n
α = α_n
β = β_n
n = n_n

! STEP 2: 初回用tangent moduli計算
C^ep = C - N^{-1}L⊗M
（ここでΔγ=0なのでC=C^e）

! STEP 3: 反復ループ
DO UNTIL (収束)
  k = k + 1

  ! 3.1: 大域平衡方程式を解く
  KU = R

  ! 3.2: δγを計算（要素レベルで）
  δγ = -N^{-1}(M:ε(δu) + g)

  ! 3.3: 変数更新
  Δu^(k+1) = Δu^(k) + δu
  Δγ^(k+1) = Δγ^(k) + δγ
  IF (Δγ < 0) Δγ = 0  ! 非負制約

  ! 3.4: BOX 1により状態更新
  ! 3.5: 内力計算
ENDDO
```

---

## 4. stress_dp_bn.fの実装戦略

### 4.1 現在の実装状況
- **良好な点**：
  - BOX 1, BOX 2の基本構造は論文に準拠
  - Drucker-Prager則への適用も適切
  - デバッグ情報の出力機能

- **改善が必要な点**：
  - N_scalar計算の数値安定性
  - 収束判定基準の調整
  - 履歴変数（histi）の管理

### 4.2 Drucker-Prager則への適用
論文のvon Mises則をDrucker-Prager則に拡張：

```fortran
! D-P parameters
eta_dp = 6sin(φ)/(√3(3-sin(φ)))
xi_dp = 6cos(φ)/(√3(3-sin(φ)))
etabar_dp = 6sin(ψ)/(√3(3-sin(ψ)))

! Yield function
g = √(1/2)|s| + η·p - ξ(σ_Y + K(α))

! L, M, N matrices for D-P
L = -2μn - κη/3·I
M = 2μn + 2μη̄/3·I
N = -(2μ + κη̄η + ξ²K')
```

### 4.3 数値安定性の確保
```fortran
! N_scalarの特異性チェック
IF (abs(N_scalar) < 1.d-12) THEN
  N_scalar = -2.d0*vmu  ! 最小値で置換
ENDIF

! Δγの非負制約
IF (deltag < 0.d0) deltag = 0.d0
```

---

## 5. phexa8.fとの連携

### 5.1 呼び出し構造
```fortran
! phexa8.f内での呼び出し（MATYPE=5の場合）
nel_current = nel
ig_current = ig
CALL stress_dp_bn(itrmax, idepg,
                  prope, sig, str, ehist,
                  ctol, vons, e_dns, p_dns,
                  ctens, g_val,
                  ierror, itr, histi)
```

### 5.2 データフロー
1. **入力**: str（ひずみ）, ehist（履歴変数）, itr（反復回数）
2. **処理**: Block Newton法による応力更新
3. **出力**: sig（応力）, ctens（接線剛性）, g_val（降伏関数値）

### 5.3 履歴変数の管理
- **ehist(20)**: 塑性ひずみ、等価塑性ひずみを保存
- **histi(50)**: Δγ、前回のひずみ、感度行列xa等を保存

---

## 6. 実装時の注意点

### 6.1 収束性の改善
1. **適応的緩和係数**: 収束が悪い場合、Δγの更新に緩和係数を導入
2. **Line search**: 降伏関数値が増加する場合の対策
3. **初期値の工夫**: 前ステップの解を有効活用

### 6.2 デバッグ戦略
```fortran
! 重要変数の監視
IF (nel_current.eq.1 .and. ig_current.eq.1) THEN
  WRITE(*,*) 'N_scalar=', N_scalar
  WRITE(*,*) 'g_val=', g_val
  WRITE(*,*) 'deltag=', deltag
ENDIF
```

### 6.3 パフォーマンス最適化
1. **行列演算の効率化**: BLAS/LAPACKの活用
2. **メモリアクセス**: 配列の連続性を考慮
3. **並列化**: OpenMPによるガウス点並列

---

## 7. 推奨される実装改良

### 7.1 短期的改善
1. **N_scalar計算の安定化**
   - 各項（2μ, κη̄η, ξ²K'）の寄与を個別に評価
   - 条件数のチェック機能追加

2. **収束判定の明確化**
   - 降伏関数の残差: |g| < tol_g
   - 平衡方程式の残差: ||R|| < tol_R

### 7.2 中長期的拡張
1. **Adaptive step control**: Δγの自動調整機能
2. **Multi-surface plasticity**: 複数降伏面への対応
3. **Large deformation**: 大変形への拡張

---

## 8. 【重要】現状と論文の相違点および具体的修正案

### 8.1 現在の実装の根本的問題

#### ❌ **最大の問題: δγの計算方法が不完全**
現在の実装では、`deltagi`（histi(1)から読み込まれた値）をそのまま使用しているが、**どこでどのように計算されているか不明**である。

**論文の正しい流れ（BOX 2 Step 3.2）**:
```fortran
δγ = -N^{-1}(M:ε(δu) + g)
```

**現在の実装（問題箇所）**:
```fortran
! line 269: 単にhistiから読み込んだ値を使用
deltag = deltagi
```

### 8.2 具体的な修正案

#### 📝 **修正案1: phexa8.fでの要素レベル縮退実装**

```fortran
! phexa8.f内に追加すべきコード（line 200付近）
! === Block Newton: 要素レベルでδγを計算 ===
if(MATYPE.eq.5 .and. itr.gt.1) then
  ! 各ガウス点でδγを計算
  do ig=1,ngaus
    ! M:ε(δu)を計算
    M_eps = 0.d0
    do jj=1,3
      do ii=1,3
        ! δuから計算されたδεを使用
        M_eps = M_eps + M_mat(ii,jj,ig) * delta_strain(ii,jj,ig)
      enddo
    enddo

    ! δγ = -N^{-1}(M:ε(δu) + g)
    if(abs(N_scalar(ig)).gt.1.d-12) then
      delta_gamma(ig) = -(M_eps + g_vals(ig)) / N_scalar(ig)
    else
      delta_gamma(ig) = 0.d0
    endif

    ! histiに格納してstress_dp_bnに渡す
    histi0(1,ig,nel) = delta_gamma(ig)
  enddo
endif
```

#### 📝 **修正案2: stress_dp_bn.fでのδγ処理改善**

```fortran
! stress_dp_bn.f内の修正（line 54-77付近）
if(itr.gt.1) then
  ! === 重要: δγを受け取る（累積値ではなく増分） ===
  delta_gamma_inc = histi(1)  ! phexa8.fから渡されたδγ増分

  ! 前回の累積値を読み込み
  deltag_prev = histi(2)

  ! 累積値を更新
  deltagi = deltag_prev + delta_gamma_inc

  ! その他の履歴変数
  alpegi  = histi(3)
  ftri    = histi(4)
  ! ... 以下同様 ...
else
  ! 初回はBOX 2初期化に従う
  deltagi = 0.d0
  deltag_prev = 0.d0
  delta_gamma_inc = 0.d0
  alpegi  = alpeg
  ftri    = 0.d0
endif
```

#### 📝 **修正案3: Stress Correctorの実装**

```fortran
! stress_dp_bn.f内に追加（line 310付近）
! === Stress Corrector計算（論文式53） ===
! σ_g = -(N^{-1}·g)·L
if(idepg.eq.1 .and. dabs(N_scalar).gt.1.d-12) then
  factor = -g_val / N_scalar
  do jj=1,3
    do ii=1,3
      psig(ii,jj) = factor * L_mat(ii,jj)
    enddo
  enddo
else
  psig(:,:) = 0.d0
endif
```

#### 📝 **修正案4: 残差力ベクトルへの反映**

```fortran
! phexa8.f内での内力計算修正（finte計算部分）
do ig=1,ngaus
  ! ... 通常の内力計算 ...

  ! === Stress Correctorの寄与を追加 ===
  if(MATYPE.eq.5) then
    ! psigはstress_dp_bnから返される
    do no=1,node
      do ii=1,3
        finte_correction = 0.d0
        do jj=1,3
          finte_correction = finte_correction
     &                     + psig(ii,jj,ig)*dndx(jj,no)*detwxyz
        enddo
        finte(3*(no-1)+ii) = finte(3*(no-1)+ii) + finte_correction
      enddo
    enddo
  endif
enddo
```

### 8.3 修正の優先順位

1. **🔴 最優先**: δγの計算メカニズム実装（修正案1,2）
2. **🟡 重要**: Stress Corrector実装（修正案3,4）
3. **🟢 推奨**: N_scalar計算の安定化強化

### 8.4 検証方法

```fortran
! デバッグコードの追加
if(nel_current.eq.1 .and. ig_current.eq.1) then
  write(*,'(A)') '=== Block Newton Debug ==='
  write(*,'(A,E12.4)') 'delta_gamma_inc=', delta_gamma_inc
  write(*,'(A,E12.4)') 'deltag_prev    =', deltag_prev
  write(*,'(A,E12.4)') 'deltagi        =', deltagi
  write(*,'(A,E12.4)') 'g_val          =', g_val
  write(*,'(A,E12.4)') 'N_scalar       =', N_scalar
  write(*,'(A,E12.4)') '||psig||       =', dsqrt(sum(psig**2))
endif
```

### 8.5 期待される効果

これらの修正により：
- ✅ 真のBlock Newton法が実現（局所反復の完全排除）
- ✅ 平衡方程式と降伏条件の同時収束
- ✅ 収束速度の向上（論文で示された2次収束）
- ✅ 大規模問題への適用性向上

---

## 9. まとめ

山本らのBlock Newton法は画期的手法だが、**現在の実装には根本的な修正が必要**である。特にδγの計算メカニズムとStress Correctorの実装が不完全であり、これらを修正することで初めて論文の理論通りの性能が発揮される。

### 実装チェックリスト（修正後）
- [x] BOX 1アルゴリズムの基本実装
- [x] BOX 2アルゴリズムの基本実装
- [x] Drucker-Prager則への拡張
- [x] phexa8.fとの基本連携
- [ ] **δγ計算メカニズムの実装（修正案1,2）**
- [ ] **Stress Correctorの実装（修正案3,4）**
- [ ] **要素レベル縮退の実装**
- [ ] **残差力ベクトルへのStress Corrector反映**
- [ ] 数値安定性の完全な保証
- [ ] 大変形問題への拡張
- [ ] 並列化の最適化

### 現在の実装の評価
**理論的完成度**: 60%（基本構造は正しいが、核心部分が未実装）
**実用性**: 30%（δγ計算メカニズムが不完全なため、真のBlock Newton法として機能していない）

---

## 参考文献
Yamamoto, T., Yamada, T., & Matsui, K. (2021). Simultaneously iterative procedure based on block Newton method for elastoplastic problems. International Journal for Numerical Methods in Engineering, 122(9), 2145-2178.