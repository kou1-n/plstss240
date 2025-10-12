# Yamamoto論文準拠チェック

## 論文: Yamamoto et al. (2021) "Simultaneously iterative procedure based on block Newton method for elastoplastic problems"

### 1. BOX 1アルゴリズム（論文のFigure 3）との対応

#### STEP 1: Compute trial stress
**論文（式73-74）:**
```
σ^tr = C^e : (ε - ε^p_n)
```

**実装（lines 174-190）:** ✅ 準拠
```fortran
c   試行偏差応力: s^tr = 2μ{dev[ε] - ε^p}
    stry = 2.d0*vmu*(str -emean*DELTA -plstrg_old)
    stno = dsqrt(stno)    ! ||s^tr||
```

#### STEP 2: Update consistency parameter
**論文（式75-76）:**
- 塑性の場合: Δγを更新
- 弾性の場合: Δγ = 0

**実装の問題点:** ❌ 不完全
- 現在の実装では、`delta_gamma_inc`はphexa8.fから受け取る（line 80）
- BOX 1内でのΔγ更新ロジックが不明確

#### STEP 3: Compute tangent moduli
**論文（式77-79）:**
```
C^ep = C - (1/N) L ⊗ M
```

**実装（lines 347-408）:** ✅ 準拠
```fortran
ctens(ii,jj,kk,ll) = C_ijkl - N_inv * L_mat(ii,jj) * M_mat(kk,ll)
```

### 2. BOX 2アルゴリズム（論文のFigure 4）との対応

#### STEP 1: Initialize
**論文:**
```
u^(0) = u_n, γ^(0) = γ_n
Δγ^(0) = 0
```

**実装（lines 221-237）:** ✅ 準拠
```fortran
if(itr.eq.1) then
    deltag = 0.d0
    delta_gamma_inc = 0.d0
    alpeg_new = alpeg_old
```

#### STEP 2: Compute global system matrices
**論文（式80-82）:**
```
K_T = Σ B^T C^ep B
L = Σ B^T L
M = Σ M^T B
N = Σ N
```

**実装の問題点:** ⚠️ 部分的
- L_mat, M_mat, N_scalarは計算されている（lines 241-259）
- しかし、グローバル組み立ては`phexa8.f`と`analys.f`で行われる

#### STEP 3: Newton update
**論文（式83）:**
```
[K_T  L] [δu]   [R_u]
[M^T  N] [δγ] = [g  ]
```

**実装の問題点:** ❌ 分散
- δγの計算は`phexa8.f`で実施（lines 269-273）
```fortran
delta_gamma_inc = -(M_eps + g_val_prev) / N_scalar
```

### 3. 主要な相違点と問題点

#### ❌ **BOX構造の分離**
- 論文: BOX 1（積分点）とBOX 2（大域）が明確に分離
- 実装: stress_dp_bn.f、phexa8.f、analys.fに分散

#### ❌ **δγ計算の位置**
- 論文: BOX 2 STEP 3で大域的に計算
- 実装: phexa8.f内で要素ごとに計算

#### ⚠️ **応力補正項psig**
- 論文には明示的な記載なし
- 実装: lines 361-378でpsig計算（Block Newton収束改善のため？）

### 4. 推奨修正事項

1. **BOX 1とBOX 2の明確な分離**
   - stress_dp_bn.f: BOX 1のみ実装
   - phexa8.f/analys.f: BOX 2実装

2. **δγ計算の統一**
   - 大域システムで一括計算する方式に変更

3. **収束判定基準の明確化**
   - 論文の式（84-85）に従った判定

### 5. 現在の実装の評価

**準拠度: 60%**

✅ 良好な点:
- 基本的な数式（L, M, N行列）は正確
- 接線剛性の計算は適切
- 初期化処理は論文準拠

❌ 改善必要:
- BOX構造が論文と異なる
- δγ計算が分散している
- 収束ロジックが不明確

### 6. 結論

現在の実装は論文の基本概念を取り入れているが、**BOXアルゴリズムの構造**が論文と異なっている。特に：

1. BOX 1（積分点処理）とBOX 2（大域Newton反復）の分離が不明確
2. consistency parameter δγの計算位置が論文と異なる
3. 応力補正項psigは論文にない追加実装

これらの相違が収束問題の原因となっている可能性がある。