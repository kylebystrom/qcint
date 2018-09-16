#!/usr/bin/env clisp 
;;;; Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>

(load "gen-code.cl")

(gen-cint "intor1.c"
  ;'("int1e_ovlp"                ( \| ))
  ;'("int1e_nuc"                 ( \| nuc \| ))
  '("int1e_kin"                 (.5 \| p dot p))
  '("int1e_ia01p"               (#C(0 1) \| nabla-rinv \| cross p))
  '("int1e_giao_irjxp"          (#C(0 1) \| r cross p))
  '("int1e_cg_irxp"             (#C(0 1) \| rc cross p))
  '("int1e_giao_a11part"        (-.5 \| nabla-rinv \| r))
  '("int1e_cg_a11part"          (-.5 \| nabla-rinv \| rc))
  '("int1e_a01gp"               (g \| nabla-rinv cross p \|))
  '("int1e_igkin"               (#C(0 .5) g \| p dot p))
  '("int1e_igovlp"              (#C(0 1) g \|))
  '("int1e_ignuc"               (#C(0 1) g \| nuc \|))
  '("int1e_pnucp"               (p* \| nuc dot p \| ))
  '("int1e_z"                   ( \| zc \| ))
  '("int1e_zz"                  ( \| zc zc \| ))
  '("int1e_r"                   ( \| rc \| ))
  '("int1e_r2"                  ( \| rc dot rc \| ))
  '("int1e_rr"                  ( \| rc rc \| ))
  '("int1e_rrr"                 ( \| rc rc rc \| ))
  '("int1e_rrrr"                ( \| rc rc rc rc \| ))
  '("int1e_z_origj"             ( \| z \| ))
  '("int1e_zz_origj"            ( \| z z \| ))
  '("int1e_r_origj"             ( \| r \| ))
  '("int1e_rr_origj"            ( \| r r \| ))
  '("int1e_r2_origj"            ( \| r dot r \| ))
  '("int1e_r4_origj"            ( \| r dot r r dot r))
  '("int1e_p4"                  ( p dot p \| p dot p ))
  ; use p* instead of p, to ignore the operator after it, then it can
  ; cross to the next p
  '("int1e_prinvp"              (p* \| rinv dot p \| ))
  '("int1e_prinvxp"             (p* \| rinv cross p \| ))
  '("int1e_pnucxp"              (p* \| nuc cross p \| ))
  '("int1e_irp"                 ( \| rc nabla \| ))
  '("int1e_irrp"                ( \| rc rc nabla \| ))
  '("int1e_irpr"                ( \| rc nabla rc \| ))
  '("int1e_gg"                  ( \| g g \|))
  '("int1e_ggkin"               ( \| g g nabla dot nabla \|))
  '("int1e_ggnuc"               ( \| g g nuc \|))
  '("int1e_grjxp"               ( \| g r cross p \|))
)

(gen-cint "intor2.c"
  ;'("int2e"                     ( \, \| \, ))
  '("int2e_ig1"                 (#C(0 1) g \, \| \, ))
  '("int2e_gg1"                 (g g \, \| \, ))
  '("int2e_g1g2"                (-1 g \, \| g \, ))
  '("int2e_p1vxp1"              (p* \, cross p \| \, )) ; SSO
  '("int2e_ip1v_rc1"            ( \, rc \| nabla-r12 \| \, ))
  '("int2e_ip1v_r1"             ( \, r  \| nabla-r12 \| \, ))
  '("int2e_ipvg1_xp1"           (g \, \| nabla-r12 cross p \| \, ))
  '("int2e_ipvg2_xp1"           (  \, \| nabla-r12 cross p \| g \, ))
  '("int1e_inuc_rcxp"           (#C(0 1) \| nuc \| rc cross p ))
  '("int1e_inuc_rxp"            (#C(0 1) \| nuc \| r cross p ))
)

(gen-cint "intor3.c"
  '("int1e_sigma"               ( \| sigma ))
  '("int1e_spsigmasp"           (sigma dot p \| sigma sigma dot p))
  '("int1e_srsr"                (sigma dot r \| sigma dot r))
  '("int1e_sr"                  (sigma dot r \|))
  '("int1e_srsp"                (sigma dot r \| sigma dot p))
  '("int1e_spsp"                (sigma dot p \| sigma dot p))
  '("int1e_sp"                  (sigma dot p \|))
  '("int1e_spnucsp"             (sigma dot p \| nuc \| sigma dot p))
  '("int1e_sprinvsp"            (sigma dot p \| rinv \| sigma dot p))
  '("int1e_srnucsr"             (sigma dot r \| nuc \| sigma dot r))
  '("int1e_sprsp"               (sigma dot p \| rc \| sigma dot p))
  '("int1e_govlp"               (g \|))
  '("int1e_gnuc"                (g \| nuc \|))
  '("int1e_cg_sa10sa01"         (.5 sigma cross rc \| sigma cross nabla-rinv \|))
  '("int1e_cg_sa10sp"           (.5 rc cross sigma \| sigma dot p))
  '("int1e_cg_sa10nucsp"        (.5 rc cross sigma \| nuc \| sigma dot p))
  '("int1e_giao_sa10sa01"       (.5 sigma cross r \| sigma cross nabla-rinv \|))
  '("int1e_giao_sa10sp"         (.5 r cross sigma \| sigma dot p))
  '("int1e_giao_sa10nucsp"      (.5 r cross sigma \| nuc \| sigma dot p))
  '("int1e_sa01sp"              (\| nabla-rinv cross sigma \| sigma dot p))
  '("int1e_spgsp"               (g sigma dot p \| sigma dot p))
  '("int1e_spgnucsp"            (g sigma dot p \| nuc \| sigma dot p))
  '("int1e_spgsa01"             (g sigma dot p \| nabla-rinv cross sigma \|))
)

(gen-cint "intor4.c"
  '("int2e_spsp1"               (sigma dot p \, sigma dot p \| \, ))
  '("int2e_spsp1spsp2"          (sigma dot p \, sigma dot p \| sigma dot p \, sigma dot p ))
  '("int2e_srsr1"               (sigma dot r \, sigma dot r \| \,))
  '("int2e_srsr1srsr2"          (sigma dot r \, sigma dot r \| sigma dot r \, sigma dot r))
  '("int2e_cg_sa10sp1"          (.5 rc cross sigma \, sigma dot p \| \,))
  '("int2e_cg_sa10sp1spsp2"     (.5 rc cross sigma \, sigma dot p \| sigma dot p \, sigma dot p ))
  '("int2e_giao_sa10sp1"        (.5 r cross sigma \, sigma dot p \| \,))
  '("int2e_giao_sa10sp1spsp2"   (.5 r cross sigma \, sigma dot p \| sigma dot p \, sigma dot p ))
  '("int2e_g1"                  (g \, \| \,))
  '("int2e_spgsp1"              (g sigma dot p \, sigma dot p \| \,))
  '("int2e_g1spsp2"             (g \, \| sigma dot p \, sigma dot p))
  '("int2e_spgsp1spsp2"         (g sigma dot p \, sigma dot p \| sigma dot p \, sigma dot p))
  '("int2e_pp1"                 (p* \, dot p \| \,))
  '("int2e_pp2"                 (\, \| p* \, dot p))
  '("int2e_pp1pp2"              (p* \, dot p \| p* \, dot p))
)

(gen-cint "dkb.c"
  ; for DKB
  '("int1e_spspsp"              (sigma dot p \| sigma dot p sigma dot p))
  '("int1e_spnuc"               (sigma dot p \| nuc \|))
  '("int2e_spv1"                (sigma dot p \, \| \,))
  '("int2e_vsp1"                (\, sigma dot p \| \,))
  '("int2e_spsp2"               (\, \| sigma dot p \, sigma dot p))
  '("int2e_spv1spv2"            (sigma dot p \, \| sigma dot p \,))
  '("int2e_vsp1spv2"            (\, sigma dot p \| sigma dot p \,))
  '("int2e_spv1vsp2"            (sigma dot p \, \| \, sigma dot p))
  '("int2e_vsp1vsp2"            (\, sigma dot p \| \, sigma dot p))
  '("int2e_spv1spsp2"           (sigma dot p \, \| sigma dot p \, sigma dot p))
  '("int2e_vsp1spsp2"           (\, sigma dot p \| sigma dot p \, sigma dot p))
)

(gen-cint "grad1.c"
  '("int1e_ipovlp"              (nabla \|))
  '("int1e_ovlpip"              (\| nabla))
  '("int1e_ipkin"               (.5 nabla \| p dot p))
  '("int1e_kinip"               (.5 \| p dot p nabla))
  '("int1e_ipnuc"               (nabla \| nuc \|))
  '("int1e_iprinv"              (nabla \| rinv \|))
  '("int1e_rinv"                (\| rinv \|))
  '("int1e_ipspnucsp"           (nabla sigma dot p \| nuc \| sigma dot p))
  '("int1e_ipsprinvsp"          (nabla sigma dot p \| rinv \| sigma dot p))
  '("int1e_ippnucp"             (p* nabla \| nuc dot p \|))
  '("int1e_ipprinvp"            (p* nabla \| rinv dot p \|))
)

(gen-cint "grad2.c"
  '("int2e_ip1"                 (nabla \, \| r12 \| \,))
  '("int2e_ip2"                 ( \, \| nabla \,))
  '("int2e_ipspsp1"             (nabla sigma dot p \, sigma dot p \| \,))
  '("int2e_ip1spsp2"            (nabla \, \| sigma dot p \, sigma dot p))
  '("int2e_ipspsp1spsp2"        (nabla sigma dot p \, sigma dot p \| sigma dot p \, sigma dot p))
  '("int2e_ipsrsr1"             (nabla sigma dot r \, sigma dot r \| \,))
  '("int2e_ip1srsr2"            (nabla \, \| sigma dot r \, sigma dot r))
  '("int2e_ipsrsr1srsr2"        (nabla sigma dot r \, sigma dot r \| sigma dot r \, sigma dot r))
)

(gen-cint "gaunt1.c"
  '("int2e_ssp1ssp2"            ( \, sigma dot p \| gaunt \| \, sigma dot p))
  '("int2e_ssp1sps2"            ( \, sigma dot p \| gaunt \| sigma dot p \,))
  '("int2e_sps1ssp2"            ( sigma dot p \, \| gaunt \| \, sigma dot p))
  '("int2e_sps1sps2"            ( sigma dot p \, \| gaunt \| sigma dot p \,))
  '("int2e_cg_ssa10ssp2"        (rc cross sigma \, \| gaunt \| \, sigma dot p))
  '("int2e_giao_ssa10ssp2"      (r cross sigma  \, \| gaunt \| \, sigma dot p))
  '("int2e_gssp1ssp2"           (g \, sigma dot p  \| gaunt \| \, sigma dot p))
)

(gen-cint "breit1.c"
  '("int2e_gauge_r1_ssp1ssp2"   ( \, sigma dot p \| breit-r1 \| \, sigma dot p))
  '("int2e_gauge_r1_ssp1sps2"   ( \, sigma dot p \| breit-r1 \| sigma dot p \,))
  '("int2e_gauge_r1_sps1ssp2"   ( sigma dot p \, \| breit-r1 \| \, sigma dot p))
  '("int2e_gauge_r1_sps1sps2"   ( sigma dot p \, \| breit-r1 \| sigma dot p \,))
  '("int2e_gauge_r2_ssp1ssp2"   ( \, sigma dot p \| breit-r2 \| \, sigma dot p))
  '("int2e_gauge_r2_ssp1sps2"   ( \, sigma dot p \| breit-r2 \| sigma dot p \,))
  '("int2e_gauge_r2_sps1ssp2"   ( sigma dot p \, \| breit-r2 \| \, sigma dot p))
  '("int2e_gauge_r2_sps1sps2"   ( sigma dot p \, \| breit-r2 \| sigma dot p \,))
)

(gen-cint "hess.c"
  '("int1e_ipipovlp"            ( nabla nabla \| ))
  '("int1e_ipovlpip"            ( nabla \| nabla ))
  '("int1e_ipipkin"             ( .5 nabla nabla \| p dot p \| ))
  '("int1e_ipkinip"             ( .5 nabla \| p dot p \| nabla ))
  '("int1e_ipipnuc"             ( nabla nabla \| nuc \| ))
  '("int1e_ipnucip"             ( nabla \| nuc \| nabla ))
  '("int1e_ipiprinv"            ( nabla nabla \| rinv \| ))
  '("int1e_iprinvip"            ( nabla \| rinv \| nabla ))
  '("int2e_ipip1"               ( nabla nabla \, \| \, ))
  '("int2e_ipvip1"              ( nabla \, nabla \| \, ))
  '("int2e_ip1ip2"              ( nabla \, \| nabla \, ))
  '("int1e_ipippnucp"           ( p* nabla nabla \| nuc dot p \| ))
  '("int1e_ippnucpip"           ( p* nabla \| nuc dot p \| nabla ))
  '("int1e_ipipprinvp"          ( p* nabla nabla \| rinv dot p \| ))
  '("int1e_ipprinvpip"          ( p* nabla \| rinv dot p \| nabla ))
  '("int1e_ipipspnucsp"         ( nabla nabla sigma dot p \| nuc sigma dot p \| ))
  '("int1e_ipspnucspip"         ( nabla sigma dot p \| nuc sigma dot p \| nabla ))
  '("int1e_ipipsprinvsp"        ( nabla nabla sigma dot p \| rinv sigma dot p \| ))
  '("int1e_ipsprinvspip"        ( nabla sigma dot p \| rinv sigma dot p \| nabla ))
)

(gen-cint "int3c2e.c"
  '("int3c2e_ip1"               (nabla \, \| ))
  '("int3c2e_ip2"               ( \, \| nabla))
  '("int3c2e_pvp1"              (p* \, dot p \|))
  '("int3c2e_pvxp1"             (p* \, cross p \|))
  '("int2c2e_ip1"               (nabla \| r12 \| ))
  '("int2c2e_ip2"               ( \| r12 \| nabla))
  '("int3c2e_ig1"               (#C(0 1) g \, \| ))
  '("int3c2e_spsp1"             (sigma dot p \, sigma dot p \| ))
  '("int3c2e_ipspsp1"           (nabla sigma dot p \, sigma dot p \| ))
  '("int3c2e_spsp1ip2"          (sigma dot p \, sigma dot p \| nabla ))
;
  '("int3c2e_ipip1"             ( nabla nabla \, \| ))
  '("int3c2e_ipvip1"            ( nabla \, nabla \| ))
  '("int3c2e_ip1ip2"            ( nabla \, \| nabla ))
  '("int2c2e_ip1ip2"            ( nabla \| r12 \| nabla))
)

(gen-cint "int3c1e.c"
  '("int3c1e_p2"                ( \, \, p dot p))
  '("int3c1e_iprinv"            ( p \, \| rinv \| ))
)

