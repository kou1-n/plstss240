1  呼出: main.f
2    call flopen(ierror, iRCM);        // 出力ファイル/ログを開く（/iodev/ を使用）
3    call initt3(...); call initt4(...); call initt6(...);   // 節点・要素・境界条件・材料読込
4    call elastc(..., Cmat, ...);      // 材料毎の弾性Cを一度生成（平面応力/ひずみは既存流儀）
5    call analys();                     // 以後は analys.f が主制御
6
7  subroutine analys()
8    初期化: u(:)=0; γ_g(:)=0; 反復カウンタ it=0
9    do istep = 1, nstep                         // 荷重/擬似時間ステップ
10       収束フラグ = .false.
11       do while (.not. 収束フラグ)             // 外側ニュートン反復
12         R(:)=0;  K(:,:)=0                      // 全体残差・全体接線のリセット
13         call forces(istep, Fext)               // 外力ベクトル（荷重・境界条件）
14
15         do e = 1, nelem                        // === 要素ループ ===
16           取り出し: ue = gather(u, lm(e,:))    // 既存の LM で要素自由度抽出
17           Re(:)=0;  Ke(:,:)=0                  // 要素残差・剛性
18           do g = 1, ngaus(e)                   // --- ガウスポイント ---
19             B_g, wJ = 形状関数微分・ヤコビアン(pquad4/8 内で算出)
20             eps_g = B_g * ue                   // 微小ひずみ
21
22             ! === DP-BlockNewton 材料: stress_dp_bn.f ===
23             call stress_dp_bn( eps_g, state_n(e,g), matp(e),
24      &                         gamma(g), sigma_g, C_alg_g,
25      &                         rcorr_g, mCeB_g, f_g, D_g, m_g, n_g )
26             ! ※ gamma(g): 反復内の現在値（inout）。state_n: 前ステップの内部変数
27
28             ! --- BN 仕様の要素組立て（rank-1補正込み）---
29             Re += wJ * ( B_g^T * (sigma_g + rcorr_g) )
30             Ke += wJ * ( B_g^T *  C_alg_g * B_g )
31
32             ! γ更新に必要な量はワークに保存
33             save_f(e,g)   = f_g
34             save_D(e,g)   = D_g
35             save_mCeB(e,g,:)= mCeB_g(:)        // (1×ndof_e)
36           end do  ! g
37
38           call assemble(Ke, Re, K, R, lm(e,:))  // 既存のアセンブリ流儀
39         end do    ! e
40
41         ! --- 全体連立 K δu = Fext - R の解法（既存の Skyline）---
42         rhs(:) = Fext(:) - R(:)
43         call solut(K, rhs, du)                  // → 内部で skylin.f を使用
44         u(:) = u(:) + du(:)
45
46         ! --- 同時反復の局所更新：局所反復なし・式で一発更新 ---
47         do e = 1, nelem
48           du_e = gather(du, lm(e,:))
49           do g = 1, ngaus(e)
50             dgamma = ( save_f(e,g) + dot( save_mCeB(e,g,:), du_e ) )
51      &                                      / save_D(e,g)
52             gamma(g) = max( 0.0d0, gamma(g) + dgamma )  // 非負制約
53           end do
54         end do
55
56         ! --- 収束判定（釣合 + 降伏残差の同時判定）---
57         res_glo = norm2(Fext - R) / max(1.0d0, norm2(Fext))
58         max_f   = max_g | save_f(e,g) |
59         収束フラグ = (res_glo <= tol_R) .and. (max_f <= tol_f)
60
61         it = it + 1
62       end do  ! while
63
64       ! === ステップ確定：内部変数の更新と出力 ===
65       do e = 1, nelem; do g = 1, ngaus(e)
66         state_{n+1}(e,g)%epsp = state_n(e,g)%epsp + gamma(g)*n_g
67         state_{n+1}(e,g)%kappa= state_n(e,g)%kappa+ gamma(g)
68         state_{n+1}(e,g)%sig  = 最新の sigma_g                     // 保存
69       end do; end do
70       call postpr(istep, u, state_{n+1}, ...)    // STS/DIS/NOR/ENE など既存出力
71     end do  ! istep
72  end subroutine analys
