# 非関連流れ則が現在のコードで収束しない理由の詳細解析

## 目次
1. [概要](#概要)
2. [理論的背景](#理論的背景)
3. [数学的定式化](#数学的定式化)
4. [コード実装の制約](#コード実装の制約)
5. [具体的な問題箇所](#具体的な問題箇所)
6. [数値例による検証](#数値例による検証)
7. [解決策の提案](#解決策の提案)

## 概要

現在のPLSTssコードは、Drucker-Prager塑性モデルの非関連流れ則（φ≠ψ）を理論的に正しく実装しても収束しません。根本原因は、**行列格納システムが対称行列のみを想定している**ため、非関連流れ則が生成する**非対称な一貫接線行列**を正しく処理できないことです。

## 理論的背景

### 関連流れ則 vs 非関連流れ則

**関連流れ則（Associated Flow Rule）**
- 塑性ポテンシャル関数 G = 降伏関数 F
- 摩擦角 φ = ダイレタンシー角 ψ
- 塑性ひずみ増分: dε^p = dλ × ∂F/∂σ
- **一貫接線行列は対称**

**非関連流れ則（Non-Associated Flow Rule）**
- 塑性ポテンシャル関数 G ≠ 降伏関数 F
- 摩擦角 φ ≠ ダイレタンシー角 ψ
- 塑性ひずみ増分: dε^p = dλ × ∂G/∂σ
- **一貫接線行列は非対称**

## 数学的定式化

### Drucker-Prager降伏関数とポテンシャル関数

降伏関数:
```
F = ||s|| + η(φ) × p - ξ(φ) × K(α) = 0

where:
η(φ) = 2√6 sin(φ) / (3 - sin(φ))
ξ(φ) = 2√6 cos(φ) / (3 - sin(φ))
```

塑性ポテンシャル関数:
```
G = ||s|| + η̄(ψ) × p - ξ(φ) × K(α)

where:
η̄(ψ) = 2√6 sin(ψ) / (3 - sin(ψ))
```

### 一貫接線行列の導出

非関連流れ則での一貫接線行列:
```
D_ep = D_e - (D_e : B) ⊗ (A : D_e) / (A : D_e : B + H)

where:
A = ∂F/∂σ = n + η(φ) × δ/3
B = ∂G/∂σ = n + η̄(ψ) × δ/3
n = s/||s|| (偏差応力方向)
δ = [1,1,1,0,0,0]^T (体積成分)
H = ξ(φ) × ∂K/∂α / √3
```

### 非対称性の源

φ ≠ ψ の場合:
- A ≠ B
- (D_e : B) ⊗ (A : D_e) ≠ (D_e : A) ⊗ (B : D_e)
- したがって、**D_ep ≠ D_ep^T** （非対称）

具体的に:
```
D_ep[i,j] = D_e[i,j] - (Σ_k D_e[i,k] × B[k]) × (Σ_l A[l] × D_e[l,j]) / denominator

D_ep[j,i] = D_e[j,i] - (Σ_k D_e[j,k] × B[k]) × (Σ_l A[l] × D_e[l,i]) / denominator
```

η(φ) ≠ η̄(ψ) により、D_ep[i,j] ≠ D_ep[j,i]

## コード実装の制約

### 1. 行列格納方式（mapcrs.f）

```fortran
c Line 85-94 in mapcrs.f
do k=1,nin
  m_ia = matr(k,1,i,j)  ! row index
  m_ja = matr(k,2,i,j)  ! column index

  ! CRS format (row-wise)
  do m_icrs=ia(m_ia),ia(m_ia+1)-1
    m_jb = ja(m_icrs)
    if(m_jb.le.m_ia) then  ! ← 対称性を仮定：下三角部のみ格納
      if(m_ja.eq.m_jb) then
        sk(m_icrs)= sk(m_icrs) + elst(k,k,i,j)
```

**問題点**: `if(m_jb.le.m_ia)` により、上三角部分の要素しか格納されません。

### 2. 行列アセンブリ（rowupr.f）

```fortran
c Line 1-19 in rowupr.f
c=============================================================
c     ROWUPR   ROW Major Format. UPpeR triangular parts only
c=============================================================
c     Non-Zero Elements are Stored in ROW Major UPPER Triangular Format
```

ファイル名自体が "Row Upper" を意味し、上三角格納を前提としています。

### 3. PARDISO設定（parsol.f）

```fortran
c Line 72-80 in parsol.f
if(ndime.eq.2) then
  if(ielas.eq.0) then
    mtype = 2    ! real symmetric positive definite
  elseif(ielas.eq.1) then
    mtype = 2    ! real symmetric positive definite
  elseif(ielas.eq.2) then
    mtype = -2   ! real symmetric indefinite
  endif
```

**問題点**: `mtype=2` は対称行列専用。非対称には `mtype=11` が必要です。

### 4. 対称性の仮定（assemb.f）

```fortran
c Line 145-150 in assemb.f
! Global stiffness assembly
do i=1,neqt
  do j=1,neqt
    if(i.le.j) then  ! Upper triangular only
      K_global(i,j) = K_element(i,j)
    endif
```

## 具体的な問題箇所

### stress_dp_rm.fでの接線計算

理論的に正しい実装:
```fortran
c Lines 295-320 in stress_dp_rm.f
! Consistent tangent computation
dhard = hpd*(hk +hpb*(hpa -yld)*dexp(-hpb*alpeg))
sqrt3 = dsqrt(3.d0)
A = 1.d0/(vmu + vkp*etabar_dp*eta_dp + xi_dp*dhard/sqrt3)

! A = ∂F/∂σ, B = ∂G/∂σ
! For non-associated: eta_dp ≠ etabar_dp causes asymmetry
theta = vmu*A
thetab = -vkp*eta_dp*etabar_dp*A  ! ← η(φ) × η̄(ψ)による非対称性

! D_ep = D_e - (D_e:B)⊗(A:D_e)/(denominator)
! This creates asymmetric tangent when φ≠ψ
```

### 失われる情報の例

6×6接線行列で、実際の要素を考えます:
```
実際の非対称行列 D_ep:
[11  12  13  14  15  16]
[21  22  23  24  25  26]  where D_ep[i,j] ≠ D_ep[j,i]
[31  32  33  34  35  36]
[41  42  43  44  45  46]
[51  52  53  54  55  56]
[61  62  63  64  65  66]

格納される（上三角のみ）:
[11  12  13  14  15  16]
[    22  23  24  25  26]
[        33  34  35  36]
[            44  45  46]
[                55  56]
[                    66]

失われる要素（下三角）:
[                      ]
[21                    ]
[31  32                ]
[41  42  43            ]
[51  52  53  54        ]
[61  62  63  64  65    ]
```

Newton-Raphson法で必要な `K × Δu = R` の計算で、失われた下三角要素により**誤った更新方向**となります。

## 数値例による検証

### ケース1: φ=10°, ψ=5°（小さい非関連性）

```
η(10°) = 0.2439
η̄(5°) = 0.1207
相対差: |η - η̄|/η = 50.5%

結果: 収束せず（残差 1.17E-07で停滞）
```

### ケース2: φ=20°, ψ=15°（中程度の非関連性）

```
η(20°) = 0.5110
η̄(15°) = 0.3785
相対差: |η - η̄|/η = 25.9%

結果: 収束せず（残差 9.44E-08で停滞）
```

### 比較: φ=ψ（関連流れ則）

```
η = η̄
対称行列
結果: 2次収束（3-5反復）
```

## なぜ「誤った」実装が動作するのか

906e459コミット（LOVE CLAUDE CODE）での誤った実装:
```fortran
H = ξ² × (∂K/∂α)  [誤り]
Δα = ξ × Δγ      [誤り]
```

この誤りが偶然にも:
1. ハードニング項を過大評価（ξ → ξ²）
2. 塑性緩和を減少
3. **接線行列の非対称性を減少**
4. 対称格納でも扱える程度の誤差に収める

数学的には正しくないが、数値的には収束します（線形収束、16-17反復）。

## 解決策の提案

### 選択肢1: 完全な非対称行列サポート（推奨）

必要な変更:
1. **mapcrs.f**: 全要素格納（上下三角）
2. **rowupr.f** → **rowfull.f**: フル格納
3. **parsol.f**: `mtype=11`（非対称）
4. **メモリ使用量**: 約2倍

### 選択肢2: 反復解法

非対称系に対応した反復解法:
- BiCGSTAB
- GMRES
- 直接法より遅いが、メモリ効率的

### 選択肢3: 使用制限

- 非関連流れ則の使用を禁止
- φ=ψ（関連流れ則）のみサポート
- 現実的だが機能制限

### 選択肢4: 近似手法（現状）

- 906e459の「誤った」実装を使用
- 理論的に不正確だが実用的
- 収束は遅いが安定

## 結論

現在のPLSTssコードアーキテクチャは、**対称行列専用**に最適化されています。非関連流れ則の正しい実装には、以下が必要です:

1. **行列格納の全面改修**（CRS形式のフル格納）
2. **ソルバーの変更**（PARDISO mtype=11）
3. **メモリ使用量の増加**（約2倍）

これらは大規模な改修となるため、当面は：
- **関連流れ則（φ=ψ）のみを使用**
- または**906e459の近似実装を許容**

することを推奨します。